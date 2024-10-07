REF__gtexhg38=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
MELT="java -jar MELT.jar"
MELTGENEBED=melt-2.2.0/add_bed_files/Hg38/Hg38.genes.bed
MEIREF=mei_list.hg38.txt
export PATH=$PATH:bowtie2-2.3.3.1/bin


# MELT indiv
$MELT IndivAnalysis \
  -w aaa \
  -bamfile $BAM \
  -h $REF__gtexhg38 \
  -t $MEIREF \
  -b $SKIPCHR > log.melt.indiv.txt 2>&1


# link files
cd group_step
while read SAMPLE; do
  bash group_step_link.sh $SAMPLE
done < sample.list


# MELT group
$MELT GroupAnalysis \
  -discoverydir group_step \
  -w group_step \
  -t $MEIREF \
  -h $REF__gtexhg38 \
  -n $MELTGENEBED > log.melt.group.txt 2>&1


# MELT genotype
cat bam.list | while read BAM; do
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
source python3-venv/env.sh
source VcfFunctions.sh

while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_samples ../all.final_comp.sort.supp.vcf tmp.$Proband.txt > tmp.vcf
  excludae_no_carrier_variant tmp.vcf > tmp2.vcf
  filter_by_supp tmp2.vcf 50     > tmp3.vcf
  AnnotSV -genomeBuild GRCh38 -SVinputFile tmp3.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done < family_members.txt
