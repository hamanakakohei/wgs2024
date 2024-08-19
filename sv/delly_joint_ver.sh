DELLY=/usr/local/genome/delly-0.8.3/delly 
REF=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
EXCL=/mira05/analysis/hamanaka/resource/delly/human.hg38.excl.tsv
MAP=/mira05/analysis/hamanaka/resource/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
BCFTOOLS=/usr/local/genome/bcftools-1.10/bin/bcftools


# step 1: SV calling 
cat ../../bam/cram1336.list | while read CRAM; do
  qsub qsub__delly_step1.sh $CRAM
  sleep 1h
done


# step 2: Merge SV sites into a unified site list
bash ./ss_delly_joint_ver_step2.sh


# step 3: Genotype this merged SV site list across all samples
cat ../../bam/cram1336.list | while read CRAM; do
  qsub qsub__delly_step3.sh $CRAM
  sleep 1m
done


# step 4: Merge all genotyped samples to get a single VCF/BCF using bcftools merge
bash ./ss_delly_joint_ver_step4.sh


# option
# skip to increase the sensitivity
# step 5: Apply the germline SV filter which requires at least 20 unrelated samples
$DELLY filter -f germline -o germline.sv.bcf merged.sv.bcf


# step 6: bcf to vcf
$BCFTOOLS view -O v -o germline.sv.geno.vcf germline.sv.geno.bcf


# step 7: add SUPP & IDLIST
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
}' germline.sv.geno.vcf > germline.sv.geno.supp.vcf


# step 8: annotate
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
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_sample0512 ../../germline.sv.geno.supp.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf > tmp2.vcf
  mv tmp2.vcf tmp4.vcf
  filter_by_supp tmp4.vcf 50     > tmp5.vcf
  $ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done

