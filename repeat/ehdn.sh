EHDN=ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo 
REF=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta 


# step 1: ehdn for each sample
do_ehdn_profile(){
  BAMPATH=$1
  REF=$2
  SAMPLE=`basename $BAMPATH`
  
  $EHDN profile --log-reads \
    --reads $BAMPATH \
    --reference $REF \
    --output-prefix $SAMPLE.ehdn \
    --min-anchor-mapq 50 \
    --max-irr-mapq 40 \
    > log.$SAMPLE.ehdn.txt 2>&1
  
  mkdir $SAMPLE.ehdn
  
  mv \
    $SAMPLE.ehdn.str_profile.json \
    $SAMPLE.ehdn.locus.tsv \
    $SAMPLE.ehdn.reads.tsv \
    $SAMPLE.ehdn.motif.tsv \
    log.$SAMPLE.ehdn.txt \
    $SAMPLE.ehdn/
}

while read CRAMPATH; do 
  do_ehdn_profile $CRAMPATH $REF
done < cram.list


# step 2: merge -> maf filtering
source python3-venv/env.sh
OUTLIER=ExpansionHunterDenovo/scripts/outlier.py
ANNOTATE=ExpansionHunterDenovo/scripts/annotate_ehdn.sh
MANIFEST=manifest.txt
HumanDb=annovar/201612/hg38
ANNOVAR=annovar2016oct/annotate_variation.pl

$EHDN merge --reference $REF --manifest $MANIFEST --output-prefix $PREFIX
$OUTLIER locus --manifest $MANIFEST --multisample-profile $PREFIX.multisample_profile.json --output $PREFIX.outlier_locus.tsv
$OUTLIER motif --manifest $MANIFEST --multisample-profile $PREFIX.multisample_profile.json --output $PREFIX.outlier_motif.tsv
$ANNOTATE --ehdn-results $PREFIX.outlier_locus.tsv --ehdn-annotated-results $PREFIX.outlier_locus.annotated.tsv --annovar-annotate-variation $ANNOVAR --annovar-humandb $HumanDb --annovar-buildver hg38

