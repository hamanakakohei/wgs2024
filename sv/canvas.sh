# step 0: prep HQ site vcf @ hamanaka/wgs/
REF=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
#for CHR in `seq 22 -1 1`; do
for CHR in X; do
  #VCF=/betelgeuse07/analysis/ncgm/FY2020-2021_G011/G011_jointcall.VQSR.chr$CHR.vcf.gz
  VCF=/betelgeuse07/analysis/ncgm/FY2020-2021_G011/G011_jointcall.chr$CHR.vcf.gz
  OUT=G011_jointcall.VQSR.chr$CHR.ForCanvas.vcf.gz
  java -jar /usr/local/bio/src/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants \
    -V $VCF \
    -o $OUT \
    -R $REF \
    --sites_only \
    --selectTypeToInclude SNP \
    --restrictAllelesTo BIALLELIC \
    --selectexpressions "AF > 0.1" 
    #--excludeFiltered \
done

cat \
  <(less  G011_jointcall.VQSR.chr1.ForCanvas.vcf.gz) \
  <(less  G011_jointcall.VQSR.chr2.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr3.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr4.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr5.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr6.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr7.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr8.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chr9.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr10.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr11.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr12.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr13.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr14.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr15.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr16.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr17.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr18.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr19.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr20.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr21.ForCanvas.vcf.gz|tail -n+3416) \
  <(less G011_jointcall.VQSR.chr22.ForCanvas.vcf.gz|tail -n+3416) \
  <(less  G011_jointcall.VQSR.chrX.ForCanvas.vcf.gz|tail -n+3399) \
  > G011_jointcall.VQSR.chrall.ForCanvas.vcf

# step 0-2: bam copy
# trio
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped | awk '$1==3' |sort| while read LINE; do
PT=`echo $LINE|awk '{print $2}'`
FA=`echo $LINE|awk '{print $3}'`
MO=`echo $LINE|awk '{print $4}'`
PTBAM=`grep $PT bam1216.list`
FABAM=`grep $FA bam1216.list`
MOBAM=`grep $MO bam1216.list`
echo $PTBAM
echo $FABAM
echo $MOBAM
cp $PTBAM .
cp $PTBAM.bai .
cp $FABAM .
cp $FABAM.bai .
cp $MOBAM .
cp $MOBAM.bai .
sleep 1s
done

# quad
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped | awk '$1==4' |grep -e 002729 -e 002734 -e 002741 -e 002758 -e 002953|sort| while read LINE; do
PT=`echo $LINE|awk '{print $2}'`
OT=`echo $LINE|awk '{print $3}'`
FA=`echo $LINE|awk '{print $4}'`
MO=`echo $LINE|awk '{print $5}'`
PTBAM=`grep $PT bam1216.list`
OTBAM=`grep $OT bam1216.list`
FABAM=`grep $FA bam1216.list`
MOBAM=`grep $MO bam1216.list`
echo $PTBAM
echo $OTBAM
echo $FABAM
echo $MOBAM
cp $PTBAM .
cp $PTBAM.bai .
cp $OTBAM .
cp $OTBAM.bai .
cp $FABAM .
cp $FABAM.bai .
cp $MOBAM .
cp $MOBAM.bai .
sleep 1s
done


# step 1
KMER=/mira05/analysis/hamanaka/resource/canvas/kmer.fa
FASTA=/mira05/analysis/hamanaka/resource/canvas/genome.fa
DBSNP=/mira05/analysis/hamanaka/resource/canvas/dbsnp.vcf
FILTER=/mira05/analysis/hamanaka/resource/canvas/filter13.bed
GENOMESIZE=/mira05/analysis/hamanaka/resource/canvas/GenomeSize.xml
GENOMEFOLDER=/mira05/analysis/hamanaka/resource/canvas/
VCF=/betelgeuse07/analysis/hamanaka/wgs/G011_jointcall.VQSR.chrall.ForCanvas.vcf.gz
export PATH=$PATH:/usr/local/bio/src/Canvas-1.40.0.1613+master_x64/

# trio
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped|awk '$1==3'|grep -e 000550 -e 000591|sort|while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  FA=`echo $LINE|awk '{print $3}'`
  MO=`echo $LINE|awk '{print $4}'`
  PTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  FABAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$FA.bam
  MOBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$MO.bam
  echo $PT
  echo $FA
  echo $MO
  echo $PTBAM
  echo $FABAM
  echo $MOBAM

  /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/dotnet \
    /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT \
    --bam=$FABAM father  $FA \
    --bam=$MOBAM mother  $MO \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$FA.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$MO.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam.bai
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$FA.bam.bai
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$MO.bam.bai
done

# quad
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped|awk '$1==4'|grep -e 002729 -e 002734 -e 002741 -e 002758 -e 002953|sort|while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  OT=`echo $LINE|awk '{print $3}'`
  FA=`echo $LINE|awk '{print $4}'`
  MO=`echo $LINE|awk '{print $5}'`
  PTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  OTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT.bam
  FABAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$FA.bam
  MOBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$MO.bam
  echo $PT
  echo $OT
  echo $FA
  echo $MO
  echo $PTBAM
  echo $OTBAM
  echo $FABAM
  echo $MOBAM

  /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/dotnet \
    /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT \
    --bam=$OTBAM other   $OT \
    --bam=$FABAM father  $FA \
    --bam=$MOBAM mother  $MO \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$FA.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$MO.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam.bai
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT.bam.bai
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$FA.bam.bai
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$MO.bam.bai
done

# duo
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped|awk '$1==2'|grep -v -e 002729 -e 002734 -e 002741 -e 002758 -e 002953|sort|while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  OT=`echo $LINE|awk '{print $3}'`
  PTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  OTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT.bam
  echo $PT
  echo $OT
  echo $PTBAM
  echo $OTBAM

  /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/dotnet \
    /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT \
    --bam=$OTBAM other   $OT \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam.bai
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT.bam.bai
done

# single
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped|awk '$1==1'|grep -v -e 002729 -e 002734 -e 002741 -e 002758 -e 002953|sort|while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  PTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  echo $PT
  echo $PTBAM

  /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/dotnet \
    /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam.bai
done

# six
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped|awk '$1==6'|grep -v -e 002729 -e 002734 -e 002741 -e 002758 -e 002953|sort|while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  OT1=`echo $LINE|awk '{print $3}'`
  OT2=`echo $LINE|awk '{print $4}'`
  OT3=`echo $LINE|awk '{print $5}'`
  OT4=`echo $LINE|awk '{print $6}'`
  OT5=`echo $LINE|awk '{print $7}'`
  PTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  OT1BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT1.bam
  OT2BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT2.bam
  OT3BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT3.bam
  OT4BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT4.bam
  OT5BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT5.bam
  echo $PT
  echo $OT1
  echo $OT2
  echo $OT3
  echo $OT4
  echo $OT5
  echo $PTBAM
  echo $OT1BAM
  echo $OT2BAM
  echo $OT3BAM
  echo $OT4BAM
  echo $OT5BAM

  /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/dotnet \
    /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT \
    --bam=$OT1BAM other  $OT1 \
    --bam=$OT2BAM other  $OT2 \
    --bam=$OT3BAM other  $OT3 \
    --bam=$OT4BAM other  $OT4 \
    --bam=$OT5BAM other  $OT5 \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT1.bam
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT2.bam
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT3.bam
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT4.bam
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT5.bam
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam.bai
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT1.bam.bai
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT2.bam.bai
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT3.bam.bai
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT4.bam.bai
 # rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT5.bam.bai
done

# five
less /betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped|awk '$1==5'|grep -v -e 002729 -e 002734 -e 002741 -e 002758 -e 002953|sort|while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  OT1=`echo $LINE|awk '{print $3}'`
  OT2=`echo $LINE|awk '{print $4}'`
  OT3=`echo $LINE|awk '{print $5}'`
  OT4=`echo $LINE|awk '{print $6}'`
  PTBAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  OT1BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT1.bam
  OT2BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT2.bam
  OT3BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT3.bam
  OT4BAM=/betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT4.bam
  echo $PT
  echo $OT1
  echo $OT2
  echo $OT3
  echo $OT4
  echo $PTBAM
  echo $OT1BAM
  echo $OT2BAM
  echo $OT3BAM
  echo $OT4BAM

  /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/dotnet \
    /usr/local/bio/src/Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT \
    --bam=$OT1BAM other  $OT1 \
    --bam=$OT2BAM other  $OT2 \
    --bam=$OT3BAM other  $OT3 \
    --bam=$OT4BAM other  $OT4 \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT1.bam
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT2.bam
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT3.bam
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT4.bam
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$PT.bam.bai
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT1.bam.bai
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT2.bam.bai
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT3.bam.bai
  #rm /betelgeuse07/analysis/hamanaka/wgs_sv/canvas/$OT4.bam.bai
done


# prep for jasmine
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
less canvas_vcfs.txt|while read VCF; do
  FAMILY=`echo $VCF|cut -d"/" -f1`
  reformat_canvas_vcf $VCF > $FAMILY.reformat.vcf
  split_canvas_family_vcf $FAMILY.reformat.vcf
done

reformat_canvas_vcf(){
  #1: CNVLEN -> SVLEN
  #2: SVTYPE no CNV -> DEL or DUP
  #3: ALT col ga "," wo fukundehainai
  #4: ID col ga "COMPLEXCNV" wo fukundehainai
  local VCF=$1
  awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^#/{
    print $0
  }$1 !~ /^#/ && $3 ~ /LOSS/ && $5 !~ /,/{
    sub("CNVLEN", "SVLEN", $8)
    sub("SVTYPE=CNV", "SVTYPE=DEL", $8)
    for(i=10;i<=NF;i++){
      if($i ~ /^\.\/1/){
        sub("^\./1", "0/1", $i)
      }
    }
    print $0
  }$1 !~ /^#/ && $3 ~ /GAIN/ && $5 !~ /,/{
    sub("CNVLEN", "SVLEN", $8)
    sub("SVTYPE=CNV", "SVTYPE=DUP", $8)
    for(i=10;i<=NF;i++){
      if($i ~ /^\.\/1/){
        sub("^\./1", "0/1", $i)
      }
    }
    print $0
  }' <(less $VCF)
}

split_canvas_family_vcf(){
  local VCF=$1
  NCOL=`less $VCF|head -n5000|awk -F"\t" '$1=="#CHROM"{print NF}'`
  for i in `seq 10 $NCOL`; do
    SAMPLE=`less $VCF|head -n5000|awk -F"\t" 'BEGIN{OFS="\t"}$1=="#CHROM"'|cut -f $i`
    echo $SAMPLE
    less $VCF|cut -f1-9,$i | awk -v SAMPLE="${SAMPLE}" 'BEGIN{OFS="\t"}$1 ~ /^#/{print $0}$1 !~ /^#/{$3 = SAMPLE"_"NR; print $0}' > $SAMPLE.vcf
    exclude_no_carrier_variant3 $SAMPLE.vcf > $SAMPLE.2.vcf
  done
}

# jasmine 
GNOMADSV=/mira05/analysis/hamanaka/resource/nstd166.GRCh38.variant_call.AddGt.vcf
source /antares01/analysis/hamanaka/function/VcfFunctions.sh
source /usr/local/genome/jasmine-20210811/env.sh
find `pwd` -maxdepth 1 -a -name 'DA*.2.vcf'|grep -v -e 2105 -e 2106 -e 2107 > AllVcfPath.list
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
less $PEDREFORMAT | tail -n+3 | while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}END{print "gnomad"}' > tmp.$Proband.txt
  select_sample0512 ../../jasmine.rename.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf     > tmp2.vcf
  #exclude_too_large_cnv      tmp2.vcf    > tmp3.vcf
  mv tmp2.vcf tmp3.vcf
  exclude_gnomadsv_only      tmp3.vcf    > tmp4.vcf
  filter_by_supp             tmp4.vcf 50 > tmp5.vcf
  $ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile tmp5.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4 #-promoterSize 2000 -REselect1 0 -REselect2 1
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python /antares01/analysis/hamanaka/function/annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  gzip tmp.vcf
  gzip tmp3.vcf
  gzip tmp4.vcf
  gzip tmp5.vcf
  cd ..
done

