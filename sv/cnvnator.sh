# cnvnator
source /usr/local/genome/sve/env.sh
export PATH=$PATH:/usr/local/genome/bowtie2-2.3.3.1/bin
#REF__gtexhg38=/betelgeuse07/analysis/hamanaka/resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
REF__gtexhg38=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

less /betelgeuse07/analysis/hamanaka/wgs_sv/bam120.list | while read BAM; do
  SAMPLE=`basename $BAM|cut -d"." -f1`
  mkdir -p $SAMPLE/$SAMPLE
  cd $SAMPLE
  sve call -t 2 -r $REF__gtexhg38 -g hg38 -a cnvnator -o $SAMPLE $BAM > log.cnvnator.txt 2>&1 &
  cd ..
  sleep 3m
done


# prep for jasmine
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
less cnvnator_vcfs.txt | tail -n+4 | while read VCF; do
  SAMPLE=`echo $VCF | awk -F"/" '{print $NF}' | sed 's/.vcf//'`
  #reformat_cnvnator_vcf $VCF > $SAMPLE.reformat.vcf
  rename_id_and_sort $SAMPLE.reformat.vcf > $SAMPLE.reformat.rename.vcf
done

reformat_cnvnator_vcf(){
  local VCF=$1
  awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^##/{
    print $0
  }$1 == "#CHROM"{
    N = split($10, A, "\/")
    $10 = A[(N-2)]
    print $0
  }$1 !~ /^#/ && $1 !~ /^HLA/ && $1 !~ /_/ && $1 !~ /EBV/{
    for(i=10;i<=NF;i++){
      if($i ~ /^\.\/1/){
        sub("^\./1", "0/1", $i)
      }
    }
    print $0
  }' <(less $VCF)
}


# jasmine 
GNOMADSV=/mira05/analysis/hamanaka/resource/nstd166.GRCh38.variant_call.AddGt.vcf
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
source /usr/local/genome/jasmine-20210811/env.sh
find `pwd` -maxdepth 1 -a -name 'DA*_S10.reformat.rename.vcf'|grep -v -e 2105 -e 2106 -e 2107 > AllVcfPath.list
echo $GNOMADSV >> AllVcfPath.list
jasmine file_list=AllVcfPath.list out_file=jasmine.vcf threads=8 --output_genotypes > jasmine.log 2>&1
rename_jasmine_vcf jasmine.vcf > jasmine.rename.vcf


# exclude bad samples
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
source /usr/local/genome/python3-venv/env.sh
export ANNOTSV=/usr/local/bio/src/AnnotSV
export PATH=/usr/local/genome/bedtools2-2.29.0/bin:$PATH
export PATH=/usr/local/genome/bcftools-1.8/bin:$PATH
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped
CONTROL=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.txt
less $PEDREFORMAT | while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  #echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  #select_sample0512 ../../jasmine.rename.vcf tmp.$Proband.txt > tmp.vcf
  #exclude_no_carrier_variant tmp.vcf     > tmp2.vcf
  ##exclude_too_large_cnv      tmp2.vcf    > tmp3.vcf
  #mv tmp2.vcf tmp3.vcf
  #exclude_gnomadsv_only      tmp3.vcf    > tmp4.vcf
  #filter_by_supp             tmp4.vcf 50 > tmp5.vcf
  #$ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 #-promoterSize 2000 -REselect1 0 -REselect2 1
  #AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  #python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  awk -F"\t" 'NR == 1 || $(NF-85) == "split"{split($(NF-74), A, "-"); if(A[1] != A[2] || A[1] !~ /intron/) print $0}' annotsv.$Proband.final.Maf0Filtered.tsv    > annotsv.$Proband.final.Maf0ExonicFiltered.tsv
  awk -F"\t" 'NR == 1 || $(NF-85) == "split"{split($(NF-74), A, "-"); if(A[1] != A[2] || A[1] !~ /intron/) print $0}' annotsv.$Proband.final.Maf0pLIFiltered.tsv > annotsv.$Proband.final.Maf0pLIExonicFiltered.tsv
  cd ..
done

