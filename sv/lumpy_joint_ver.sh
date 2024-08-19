BED=/mira05/analysis/hamanaka/resource/smoove/exclude.cnvnator_100bp.GRCh38.20170403.bed
REF=/betelgeuse07/analysis/hamanaka/resource/resources_broad_hg38_v0_Homo_sapiens_assembly38
INPUT=/betelgeuse07/analysis/hamanaka/wgs_sv/lumpy
SMOOVE="docker run --rm -v $INPUT:/data brentp/smoove smoove"
BEDFILE=exclude.cnvnator_100bp.GRCh38.20170403.bed
REFFILE=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta


#1. genotype each sample
cat /betelgeuse07/analysis/hamanaka/bam/cram1336.list | while read CRAM; do
  DIR=$(dirname $CRAM)
  cd $DIR
  cp $BED      .
  cp ${REF}.*  .
  SMOOVE="docker run --rm -v $DIR:/data brentp/smoove smoove"
  SAMPLE=`echo $CRAM|awk -F"/" '{print $NF}'|cut -d"." -f1`
  echo $SMOOVE
  echo $SAMPLE
  echo $CRAM

  $SMOOVE call \
    --outdir /data/results-smoove/ \
    --exclude /data/$BEDFILE \
    --name $SAMPLE \
    --fasta /data/$REFFILE \
    -p 1 \
    --genotype /data/$SAMPLE.cram \
    > /betelgeuse07/analysis/hamanaka/wgs_sv/lumpy/log.lumpy.step1.$SAMPLE.txt 2>&1

  rm exclude.cnvnator_100bp.GRCh38.20170403.bed
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.dict
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.alt
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.amb
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.ann
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.bwt
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.pac
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.sa
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.index
  rm results-smoove/$SAMPLE.split.bam.orig.bam
  rm results-smoove/$SAMPLE.split.bam
  rm results-smoove/$SAMPLE.split.bam.bai
  rm results-smoove/$SAMPLE.disc.bam.orig.bam
  rm results-smoove/$SAMPLE.disc.bam
  rm results-smoove/$SAMPLE.disc.bam.bai
  rm results-smoove/$SAMPLE.histo
  rm results-smoove/$SAMPLE-lumpy-cmd.sh
done


#2. Get the union of sites across all samples --> merged.sites.vcf.gz
bash ./lumpy_joint_step2.sh 


#3. genotype each sample at all merged sites
cat /betelgeuse07/analysis/hamanaka/bam/cram1336.list | while read CRAM; do
  DIR=$(dirname $CRAM)
  cd $DIR
  cp /betelgeuse07/analysis/hamanaka/wgs_sv/lumpy/merged.sites.vcf.gz .
  cp ${REF}.*  .
  SMOOVE="docker run --rm -v $DIR:/data brentp/smoove smoove"
  SAMPLE=`echo $CRAM|awk -F"/" '{print $NF}'|cut -d"." -f1`
  echo $SMOOVE
  echo $SAMPLE
  echo $CRAM
  $SMOOVE genotype \
    -d -x -p 1 \
    --name $SAMPLE-joint \
    --outdir /data/results-genotyped/ \
    --fasta /data/$REFFILE \
    --vcf /data/merged.sites.vcf.gz \
    /data/$SAMPLE.cram
  
  rm merged.sites.vcf.gz
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.dict
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.alt
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.amb
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.ann
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.bwt
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.pac
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.sa
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai
  rm resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.index
  sleep 1s
done


#4. paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.
bash ./lumpy_joint_step4.sh 


#5. add SUPP & IDLIST
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
}' <(zcat merged.smoove.square.vcf.gz) > merged.smoove.square.supp.vcf


#6. annotate
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
  cd ${Proband}
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_sample0512 ../../merged.smoove.square.supp.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf  > tmp2.vcf
  mv tmp2.vcf tmp4.vcf
  filter_by_supp          tmp4.vcf 50 > tmp5.vcf
  $ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done


