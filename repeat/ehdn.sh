EHDN=/usr/local/genome/ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo 
REF=/antares01/analysis/hamanaka/resource/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta 


# step 1: ehdn for each sample
do_ehdn_profile(){
  BAMPATH=$1
  REF=$2
  WESWGS=$3
  SAMPLE=`basename $BAMPATH`
  
  $EHDN profile --log-reads \
    --reads $BAMPATH \
    --reference $REF \
    --output-prefix $SAMPLE.$WESWGS.ehdn \
    --min-anchor-mapq 50 \
    --max-irr-mapq 40 \
    > log.$SAMPLE.$WESWGS.ehdn.txt 2>&1
  
  mkdir $SAMPLE.$WESWGS.ehdn
  
  mv \
    $SAMPLE.$WESWGS.ehdn.str_profile.json \
    $SAMPLE.$WESWGS.ehdn.locus.tsv \
    $SAMPLE.$WESWGS.ehdn.reads.tsv \
    $SAMPLE.$WESWGS.ehdn.motif.tsv \
    log.$SAMPLE.$WESWGS.ehdn.txt \
    $SAMPLE.$WESWGS.ehdn/
}

cat mira03_cram.txt | while read BAMPATH; do 
  do_ehdn_profile $BAMPATH  $REF wgs
done 


# step 2: merge -> maf filtering
source /usr/local/genome/python3-venv/env.sh
OUTLIER=/betelgeuse04/analysis/hamanaka/resource/ExpansionHunterDenovo/scripts/outlier.py
ANNOTATE=/betelgeuse04/analysis/hamanaka/resource/ExpansionHunterDenovo/scripts/annotate_ehdn.sh
HELPER=/betelgeuse04/analysis/hamanaka/resource/Fazal2020Scripts/EHDn-v0.6.2_HelperScripts
MANIFEST=/betelgeuse07/analysis/hamanaka/ehdn/manifest1333.noDA0000002674.txt
CONTROL=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.noDA0000002674.txt
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.noDA0000002674.ped 
PREFIX=sample1333
HumanDb=/betelgeuse01/system/db/annovar/201612/hg38
ANNOVAR=/usr/local/bio/src/annovar2016oct/annotate_variation.pl

$EHDN merge --reference $REF --manifest $MANIFEST --output-prefix $PREFIX
$OUTLIER locus --manifest $MANIFEST --multisample-profile $PREFIX.multisample_profile.json --output $PREFIX.outlier_locus.tsv
$OUTLIER motif --manifest $MANIFEST --multisample-profile $PREFIX.multisample_profile.json --output $PREFIX.outlier_motif.tsv
$ANNOTATE --ehdn-results $PREFIX.outlier_locus.tsv --ehdn-annotated-results $PREFIX.outlier_locus.annotated.tsv --annovar-annotate-variation $ANNOVAR --annovar-humandb $HumanDb --annovar-buildver hg38

