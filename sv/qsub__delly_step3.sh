#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V


CRAM=$1
SAMPLE=`basename $CRAM | cut -d"." -f1`

DELLY=delly
REF=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
EXCL=human.hg38.excl.tsv

# step 3: Genotype this merged SV site list across all samples
$DELLY call -g $REF -v sites.sv.bcf -o $SAMPLE.sv.geno.bcf -x $EXCL $CRAM > log.delly.$SAMPLE.step3.txt 2>&1

