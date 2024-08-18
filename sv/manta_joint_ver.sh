REF__gtexhg38=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped

#awk 'NF==5' $PEDREFORMAT|grep -v -f tmp.manta.sumi.txt|while read LINE; do
grep DA0000000406 $PEDREFORMAT | while read LINE; do
  MEMB1=`echo $LINE|awk '{print $2}'`
  MEMB2=`echo $LINE|awk '{print $3}'`
  MEMB3=`echo $LINE|awk '{print $4}'`
  #MEMB4=`echo $LINE|awk '{print $5}'`
  MEMB1CRAM=`grep $MEMB1 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  MEMB2CRAM=`grep $MEMB2 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  MEMB3CRAM=`grep $MEMB3 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  #MEMB4CRAM=`grep $MEMB4 /betelgeuse07/analysis/ncgm/PathCram__ncgm.txt`
  PT=`echo $MEMB1CRAM|awk -F"/" '{print $NF}'|cut -d"." -f1`  
  echo $PT
  echo $MEMB1
  echo $MEMB2
  echo $MEMB3
  #echo $MEMB4
  echo $MEMB1CRAM
  echo $MEMB2CRAM
  echo $MEMB3CRAM
  #echo $MEMB4CRAM
  
  /usr/local/genome/manta-1.5.0/bin/configManta.py \
    --bam $MEMB1CRAM \
    --bam $MEMB2CRAM \
    --bam $MEMB3CRAM \
    --referenceFasta $REF__gtexhg38 \
    --runDir ${PT}family
    #--bam $MEMB4CRAM \

  ./${PT}family/runWorkflow.py -m local -g 60 --quiet
  sleep 2s
done

# for INV
#MantaVcfList=/betelgeuse07/analysis/hamanaka/wgs_sv/manta/sample_mantaVcf.txt
MantaVcfList=/betelgeuse07/analysis/hamanaka/wgs_sv/manta/sample_mantaVcf2.txt
REF__gtexhg38=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
less $MantaVcfList | while read LINE; do
  SAMPLE=`echo $LINE|awk '{print $1}'`; VCF=`echo $LINE|awk '{print $2}'`
  #qsub /mira05/analysis/hamanaka/wgs_sv_ycu/manta/qsub__mantaOnlyInv.sh             $SAMPLE $REF__gtexhg38 $VCF
  qsub /mira05/analysis/hamanaka/wgs_sv_ycu/manta/qsub__mantaOnlyInv_ForErrorVcf.sh $SAMPLE $REF__gtexhg38 $VCF
  sleep 3s
done

# prep for jasmine
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
less manta_family_vcf.txt | tail -n+3 | while read VCF; do
  split_manta_family_vcf $VCF
done

split_manta_family_vcf(){
  local VCF=$1
  NCOL=`less $VCF|head -n5000|awk -F"\t" '$1=="#CHROM"{print NF}'`
  for i in `seq 10 $NCOL`; do
    SAMPLE=`less $VCF|head -n5000|awk -F"\t" 'BEGIN{OFS="\t"}$1=="#CHROM"'|cut -f $i`
    echo $SAMPLE
    less $VCF|cut -f1-9,$i | awk -v SAMPLE="${SAMPLE}" 'BEGIN{OFS="\t"}$1 ~ /^#/{print $0}$1 !~ /^#/{$3 = SAMPLE"_"NR; print $0}' > $SAMPLE.split.vcf
    exclude_no_carrier_variant $SAMPLE.split.vcf > $SAMPLE.split.2.vcf
  done
}

# rename & merge
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
export PATH=/antares01/analysis/hamanaka/function:$PATH
ls DA*.split.2.vcf |  while read VCF; do 
  SAMPLE=`echo $VCF|cut -d"." -f1`                                                               
  rename_id_and_sort $VCF > $SAMPLE.split.2.rename.vcf
done

# jasmine merge
GNOMADSV=/mira05/analysis/hamanaka/resource/nstd166.GRCh38.variant_call.AddGt.vcf
source /usr/local/genome/jasmine-20210811/env.sh
# DA000002105/6/7 ha nozoku, izenno jasmine.rename.vcf tokaga mazattenaiyouni chuui
find `pwd` -maxdepth 1 -a -name '*rename.vcf'|grep -v -e 2105 -e 2106 -e 2107 > AllVcfPath.list
echo $GNOMADSV >> AllVcfPath.list
jasmine file_list=AllVcfPath.list out_file=jasmine.vcf threads=8 --output_genotypes > jasmine.log 2>&1
rename_jasmine_vcf jasmine.vcf > jasmine.rename.vcf

source /antares01/analysis/hamanaka/function/VcfFunctions.sh
source /usr/local/genome/python3-venv/env.sh
export ANNOTSV=/usr/local/bio/src/AnnotSV
export PATH=/usr/local/genome/bedtools2-2.29.0/bin:$PATH
export PATH=/usr/local/genome/bcftools-1.8/bin:$PATH
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped
CONTROL=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.txt
less $PEDREFORMAT | tail -n+3 | while read LINE; do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_sample0512 ../../jasmine.rename.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf > tmp2.vcf
  #exclude_too_large_cnv tmp2.vcf > tmp3.vcf
  mv tmp2.vcf tmp3.vcf
  exclude_gnomadsv_only tmp3.vcf > tmp4.vcf
  filter_by_supp tmp4.vcf 50     > tmp5.vcf
  $ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 #-promoterSize 2000 -REselect1 0 -REselect2 1
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done
