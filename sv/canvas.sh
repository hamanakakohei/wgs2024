source VcfFunctions.sh
REF=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

# step 0: prep HQ site vcf
for CHR in `seq 1 22` X; do
  VCF=G011_jointcall.VQSR.chr$CHR.vcf.gz
  OUT=G011_jointcall.VQSR.chr$CHR.ForCanvas.vcf.gz
  java -jar GenomeAnalysisTK.jar -T SelectVariants \
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
  <(zcat  G011_jointcall.VQSR.chr1.ForCanvas.vcf.gz) \
  <(zcat  G011_jointcall.VQSR.chr2.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr3.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr4.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr5.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr6.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr7.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr8.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chr9.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr10.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr11.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr12.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr13.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr14.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr15.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr16.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr17.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr18.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr19.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr20.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr21.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat G011_jointcall.VQSR.chr22.ForCanvas.vcf.gz|tail -n+3416) \
  <(zcat  G011_jointcall.VQSR.chrX.ForCanvas.vcf.gz|tail -n+3399) \
  > G011_jointcall.VQSR.chrall.ForCanvas.vcf


# step 1
KMER=canvas/kmer.fa
FASTA=canvas/genome.fa
DBSNP=canvas/dbsnp.vcf
FILTER=canvas/filter13.bed
GENOMESIZE=canvas/GenomeSize.xml
GENOMEFOLDER=canvas/
VCF=G011_jointcall.VQSR.chrall.ForCanvas.vcf.gz
export PATH=$PATH:Canvas-1.40.0.1613+master_x64/

# trio
while read LINE; do
  PT=`echo $LINE|awk '{print $2}'`
  FA=`echo $LINE|awk '{print $3}'`
  MO=`echo $LINE|awk '{print $4}'`

  Canvas-1.40.0.1613+master_x64/dotnet \
    Canvas-1.40.0.1613+master_x64/Canvas.dll \
    SmallPedigree-WGS \
    --bam=$PTBAM proband $PT_BAM_PATH \
    --bam=$FABAM father  $FA_BAM_PATH \
    --bam=$MOBAM mother  $MO_BAM_PATH \
    -r $KMER \
    -g $GENOMEFOLDER \
    --population-b-allele-vcf $VCF \
    -f $FILTER \
    -o ${PT}family 

done < family_members.txt


# prep for jasmine
while read VCF; do
  FAMILY=`echo $VCF|cut -d"/" -f1`
  reformat_canvas_vcf $VCF > $FAMILY.reformat.vcf
  split_canvas_family_vcf $FAMILY.reformat.vcf
done < canvas_vcfs.txt


# jasmine 
source jasmine-20210811/env.sh
find `pwd` -maxdepth 1 -a -name 'DA*.2.vcf' > AllVcfPath.list
jasmine file_list=AllVcfPath.list out_file=jasmine.vcf threads=8 --output_genotypes > jasmine.log 2>&1
rename_jasmine_vcf jasmine.vcf > jasmine.rename.vcf


# exclude bad samples
source python3-venv/env.sh

while read LINE;do
  Proband=`echo $LINE|awk '{print $2}'`
  mkdir $Proband
  cd $Proband
  echo $LINE|awk '{for(i=2;i<=NF;i++)print $i}' > tmp.$Proband.txt
  select_samples ../../jasmine.rename.vcf tmp.$Proband.txt > tmp.vcf
  exclude_no_carrier_variant tmp.vcf     > tmp2.vcf
  filter_by_supp             tmp2.vcf 50 > tmp3.vcf
  AnnotSV -genomeBuild GRCh38 -SVinputFile tmp3.vcf -outputFile annotsv.$Proband.tsv -svtBEDcol 4
  AnnotsvFile=`ls 2023*_AnnotSV/annotsv.DA0*.tsv|grep -v unanno`
  python annotsv_filter.py -annotsv_file $AnnotsvFile -sample $Proband -control $CONTROL
  cd ..
done < family_members.txt

