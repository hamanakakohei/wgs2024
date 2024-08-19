BamList=/betelgeuse07/analysis/hamanaka/wgs_sv/melt/bam60.list
MantaVcfList=/betelgeuse07/analysis/hamanaka/wgs_sv/manta/sample_MantaVcf3.txt
REF__gtexhg38=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
MELT="java -jar /usr/local/genome/melt-2.2.0/MELT.jar"
MELTGENEBED=/usr/local/genome/melt-2.2.0/add_bed_files/Hg38/Hg38.genes.bed
MEIREF=/mira05/analysis/hamanaka/resource/mei_list.hg38.txt
export PATH=$PATH:/usr/local/genome/bowtie2-2.3.3.1/bin
SKIPCHR=/mira05/analysis/hamanaka/wgs_sv2/SkipChr.txt


# MELT indiv
$MELT IndivAnalysis \
  -w aaa \
  -bamfile /mira03/tmp_hamanaka/DA0000007155.bam \
  -h $REF__gtexhg38 \
  -t $MEIREF \
  -b $SKIPCHR > log.melt.indiv.txt 2>&1


# link files
cd group_step
less ../sample1340.list|while read SAMPLE; do
  bash group_step_link.sh $SAMPLE
done


# MELT group
$MELT GroupAnalysis \
  -discoverydir group_step \
  -w group_step \
  -t $MEIREF \
  -h $REF__gtexhg38 \
  -n $MELTGENEBED > log.melt.group.txt 2>&1


# MELT genotype
less bam1336.list|while read BAM; do
  $MELT Genotype \
    -bamfile $BAM \
    -t $MEIREF \
    -h $REF__gtexhg38 \
    -w genotype_step \
    -p group_step 
done > log.melt.genotype.txt 2>&1


# MELT MakeVCF
$MELT MakeVCF \
  -genotypingdir genotype_step \
  -h $REF__gtexhg38 \
  -t $MEIREF \
  -w make_vcf_step \
  -p group_step > log.melt.make_vcf.txt 2>&1


# cat vcfs
cat make_vcf_step/ALU.final_comp.vcf \
  <(awk '$1!~/^#/' make_vcf_step/HERVK.final_comp.vcf) \ 
  <(awk '$1!~/^#/' make_vcf_step/LINE1.final_comp.vcf) \ 
  <(awk '$1!~/^#/' make_vcf_step/SVA.final_comp.vcf  ) \ 
  > all.final_comp.vcf


# sort vcf
cat \
  <(awk '$1 ~/^#/' all.final_comp.vcf) \
  <(awk '$1!~/^#/' all.final_comp.vcf|sort -k1,1 -k2,2n) \
  > all.final_comp.sort.vcf


# add SUPP & IDLIST
awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^##/{
  print $0
}$1 == "#CHROM"{
  for(i=10; i<=NF; i++){
    COL_SAMPLE[i] = $i
  }
  print $0
}$1 !~ /^#/{
  SUPP = 0
  IDLIST = "sampleXXX_1"
  for(i=10; i<=NF; i++){
    SAMPLE = COL_SAMPLE[i]
    if($i !~ /^0\/0/ && $i !~ /\.\/\./){
      SUPP = SUPP + 1
      IDLIST = IDLIST","SAMPLE"_1"
    }
  }
  $8 = $8";SUPP="SUPP";IDLIST="IDLIST
  print $0
}' all.final_comp.sort.vcf > all.final_comp.sort.supp.vcf 


# annotate
source /usr/local/genome/python3-venv/env.sh
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
export ANNOTSV=/usr/local/bio/src/AnnotSV
export PATH=/usr/local/genome/bedtools2-2.29.0/bin:$PATH
export PATH=/usr/local/genome/bcftools-1.8/bin:$PATH
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped
CONTROL=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.txt
cat $PEDREFORMAT | while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_sample0512 ../all.final_comp.sort.supp.vcf tmp.$Proband.txt > tmp.vcf
  excludae_no_carrier_variant tmp.vcf > tmp2.vcf
  mv tmp2.vcf tmp4.vcf
  filter_by_supp tmp4.vcf 50     > tmp5.vcf
  $ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 #-promoterSize 2000 -REselect1 0 -REselect2 1
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done
