# cnvnator
source VcfFunctions.sh
source sve/env.sh
export PATH=$PATH:bowtie2-2.3.3.1/bin
REF__gtexhg38=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

while read BAM; do
  SAMPLE=`basename $BAM|cut -d"." -f1`
  mkdir -p $SAMPLE/$SAMPLE
  cd $SAMPLE
  sve call -t 2 -r $REF__gtexhg38 -g hg38 -a cnvnator -o $SAMPLE $BAM > log.cnvnator.txt 2>&1 &
  cd ..
  sleep 3m
done < bamlist


# prep for jasmine
while read VCF; do
  SAMPLE=`echo $VCF | awk -F"/" '{print $NF}' | sed 's/.vcf//'`
  reformat_cnvnator_vcf $VCF > $SAMPLE.reformat.vcf
  rename_id_and_sort $SAMPLE.reformat.vcf > $SAMPLE.reformat.rename.vcf
done < cnvnator_vcf.list


# jasmine 
source jasmine-20210811/env.sh
find `pwd` -maxdepth 1 -a -name '*.reformat.rename.vcf' > AllVcfPath.list
jasmine file_list=AllVcfPath.list out_file=jasmine.vcf threads=8 --output_genotypes > jasmine.log 2>&1
rename_jasmine_vcf jasmine.vcf > jasmine.rename.vcf


# exclude bad samples
source python3-venv/env.sh

while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}' > tmp.$Proband.txt
  select_samples ../jasmine.rename.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf     > tmp2.vcf
  filter_by_supp             tmp2.vcf 50 > tmp3.vcf
  AnnotSV -genomeBuild GRCh38 -SVinputFile tmp3.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done < family_members.txt

