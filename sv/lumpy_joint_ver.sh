BED=smoove/exclude.cnvnator_100bp.GRCh38.20170403.bed
REF=resources_broad_hg38_v0_Homo_sapiens_assembly38
SMOOVE="docker run --rm -v $INPUT:/data brentp/smoove smoove"
BEDFILE=exclude.cnvnator_100bp.GRCh38.20170403.bed
REFFILE=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta


# 1. genotype each sample
while read CRAM; do
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
    > log.lumpy.step1.$SAMPLE.txt 2>&1

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
done < cram.list


#2. Get the union of sites across all samples --> merged.sites.vcf.gz
bash ./lumpy_joint_step2.sh 


#3. genotype each sample at all merged sites
while read CRAM; do
  DIR=$(dirname $CRAM)
  cd $DIR
  cp ../merged.sites.vcf.gz .
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
done < cram.list


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
source python3-venv/env.sh
source VcfFunctions.sh

while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd ${Proband}
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_samples ../merged.smoove.square.supp.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf  > tmp2.vcf
  filter_by_supp          tmp2.vcf 50 > tmp3.vcf
  AnnotSV -genomeBuild GRCh38 -SVinputFile tmp3.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done < family_members.list


