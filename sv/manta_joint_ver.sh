source VcfFunctions.sh
REF__gtexhg38=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta


# genotype each sample
while read LINE; do
  MEMB1=`echo $LINE|awk '{print $2}'`
  MEMB2=`echo $LINE|awk '{print $3}'`
  MEMB3=`echo $LINE|awk '{print $4}'`
  #MEMB4=`echo $LINE|awk '{print $5}'`
  MEMB1CRAM=`grep $MEMB1 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  MEMB2CRAM=`grep $MEMB2 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  MEMB3CRAM=`grep $MEMB3 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  #MEMB4CRAM=`grep $MEMB4 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  PT=`echo $MEMB1CRAM|awk -F"/" '{print $NF}'|cut -d"." -f1`  
  
  bin/configManta.py \
    --bam $MEMB1CRAM \
    --bam $MEMB2CRAM \
    --bam $MEMB3CRAM \
    --referenceFasta $REF__gtexhg38 \
    --runDir ${PT}family

  ./${PT}family/runWorkflow.py -m local -g 60 --quiet
  sleep 2s
done < family_members.txt


# for INV
while read LINE; do
  SAMPLE=`echo $LINE|awk '{print $1}'`; VCF=`echo $LINE|awk '{print $2}'`
  qsub qsub__mantaOnlyInv.sh $SAMPLE $REF__gtexhg38 $VCF
  sleep 3s
done < $MantaVcfList


# prep for jasmine
while read VCF; do
  split_manta_family_vcf $VCF
done < manta_family_vcf.txt


# rename & merge
ls *.split.2.vcf |  while read VCF; do 
  SAMPLE=`echo $VCF|cut -d"." -f1`                                                               
  rename_id_and_sort $VCF > $SAMPLE.split.2.rename.vcf
done


# jasmine merge
source jasmine-20210811/env.sh
find `pwd` -maxdepth 1 -a -name '*split.2.rename.vcf' > AllVcfPath.list
jasmine file_list=AllVcfPath.list out_file=jasmine.vcf threads=8 --output_genotypes > jasmine.log 2>&1
rename_jasmine_vcf jasmine.vcf > jasmine.rename.vcf


# anntation
source python3-venv/env.sh

while read LINE; do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_samples ../../jasmine.rename.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf > tmp2.vcf
  filter_by_supp tmp2.vcf 50         > tmp3.vcf
  AnnotSV -genomeBuild GRCh38 -SVinputFile tmp3.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done < family_members.txt
