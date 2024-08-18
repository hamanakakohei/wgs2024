DELLY=/usr/local/genome/delly-0.8.3/delly 
REF=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
EXCL=/mira05/analysis/hamanaka/resource/delly/human.hg38.excl.tsv
MAP=/mira05/analysis/hamanaka/resource/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
BCFTOOLS=/usr/local/genome/bcftools-1.10/bin/bcftools

# step 1: SV calling 
less ../../bam/cram1336.list | grep -f tmp.remainsample.txt | while read CRAM; do
  qsub ../qsub__delly_step1.sh $CRAM
  sleep 1h
done

# step 2: Merge SV sites into a unified site list
# see ss_delly_joint_ver_step2.sh
$DELLY merge -o sites.sv.bcf \
  DA0000002148.sv.bcf \
  DA0000002149.sv.bcf \
  DA0000002150.sv.bcf 

# step 3: Genotype this merged SV site list across all samples
less ../../bam/cram1336.list | grep 007226 | while read CRAM; do
  qsub ../qsub__delly_step3.sh $CRAM
  sleep 1m
done
#$DELLY call -g $REF -v sites.sv.bcf -o DA0000002148.sv.geno.bcf -x $EXCL /betelgeuse07/analysis/ncgm/data05/G011/DA0000002148/parabricks/DA0000002148.cram > log.2148.step2.txt 2>&1 &
#$DELLY call -g $REF -v sites.sv.bcf -o DA0000002149.sv.geno.bcf -x $EXCL /betelgeuse07/analysis/ncgm/data05/G011/DA0000002149/parabricks/DA0000002149.cram > log.2149.step2.txt 2>&1 &
#$DELLY call -g $REF -v sites.sv.bcf -o DA0000002150.sv.geno.bcf -x $EXCL /betelgeuse07/analysis/ncgm/data05/G011/DA0000002150/parabricks/DA0000002150.cram > log.2150.step2.txt 2>&1 &

# step 4: Merge all genotyped samples to get a single VCF/BCF using bcftools merge
# see ss_delly_joint_ver_step4.sh
$BCFTOOLS merge -m id -O b -o merged.sv.bcf \
  DA0000002148.sv.geno.bcf \
  DA0000002149.sv.geno.bcf \
  DA0000002150.sv.geno.bcf

# step 5: Apply the germline SV filter which requires at least 20 unrelated samples
$DELLY filter -f germline -o germline.sv.bcf merged.sv.bcf

# bcf to vcf
$BCFTOOLS view -O v -o germline.sv.geno.vcf germline.sv.geno.bcf

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
}' germline.sv.geno.vcf > germline.sv.geno.supp.vcf


# exclude bad samples
source /usr/local/genome/python3-venv/env.sh
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
export ANNOTSV=/usr/local/bio/src/AnnotSV
export PATH=/usr/local/genome/bedtools2-2.29.0/bin:$PATH
export PATH=/usr/local/genome/bcftools-1.8/bin:$PATH
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped
CONTROL=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.txt
less $PEDREFORMAT | tail -n+3 | while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_sample0512 ../../germline.sv.geno.supp.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf > tmp2.vcf
  #exclude_too_large_cnv tmp2.vcf > tmp3.vcf
  #exclude_gnomadsv_only tmp3.vcf > tmp4.vcf
  mv tmp2.vcf tmp4.vcf
  filter_by_supp tmp4.vcf 50     > tmp5.vcf
  $ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 #-promoterSize 2000 -REselect1 0 -REselect2 1
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done



##Delly uses GC and mappability fragment correction to call CNVs. This requires a mappability map.
##Call CNVs for each sample and optionally refine breakpoints using delly SV calls
#$DELLY cnv -o DA0000002148.cn.bcf -g $REF -m $MAP -l DA0000002148.sv.bcf $CRAM
#$DELLY cnv -o DA0000002149.cn.bcf -g $REF -m $MAP -l DA0000002149.sv.bcf $CRAM
#$DELLY cnv -o DA0000002150.cn.bcf -g $REF -m $MAP -l DA0000002150.sv.bcf $CRAM
#
##Merge CNVs into a unified site list
#$DELLY merge -e -p -o sites.cn.bcf -m 1000 -n 100000 \
#  DA0000002148.cn.bcf \
#  DA0000002149.cn.bcf \
#  DA0000002150.cn.bcf
#
##Genotype CNVs for each sample
#$DELLY cnv -u -v sites.cn.bcf -g $REF -m $MAP -o DA0000002148.cn.geno.bcf /betelgeuse07/analysis/ncgm/data05/G011/DA0000002148/parabricks/DA0000002148.cram
#$DELLY cnv -u -v sites.cn.bcf -g $REF -m $MAP -o DA0000002149.cn.geno.bcf /betelgeuse07/analysis/ncgm/data05/G011/DA0000002149/parabricks/DA0000002149.cram
#$DELLY cnv -u -v sites.cn.bcf -g $REF -m $MAP -o DA0000002150.cn.geno.bcf /betelgeuse07/analysis/ncgm/data05/G011/DA0000002150/parabricks/DA0000002150.cram
#
##Merge genotypes using bcftools
#$BCFTOOLS merge -m id -O b -o merged.cn.bcf \
#  DA0000002148.cn.geno.bcf \
#  DA0000002149.cn.geno.bcf \
#  DA0000002150.cn.geno.bcf 
#
##Filter for germline CNVs
#$DELLY classify -f germline -o filtered.cn.bcf merged.cn.bcf
#
##Optional: Plot copy-number distribution for large number of samples (>>100)
#bcftools query -f "%ID[\t%RDCN]\n" filtered.bcf > plot.tsv
#
#Rscript R/cnv.R plot.tsv
#
#
##$DELLY call -g $REF -o $SAMPLE.bcf -x $EXCL $CRAM
##$DELLY call -g $REF -o DA0000002148.sv.bcf -x $EXCL /betelgeuse07/analysis/ncgm/data05/G011/DA0000002148/parabricks/DA0000002148.cram
##$DELLY call -g $REF -o DA0000002149.sv.bcf -x $EXCL /betelgeuse07/analysis/ncgm/data05/G011/DA0000002149/parabricks/DA0000002149.cram > log.2149.txt 2>&1 &
##$DELLY call -g $REF -o DA0000002150.sv.bcf -x $EXCL /betelgeuse07/analysis/ncgm/data05/G011/DA0000002150/parabricks/DA0000002150.cram > log.2150.txt 2>&1 &
