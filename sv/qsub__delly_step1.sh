#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V


CRAM=$1
SAMPLE=`basename $CRAM | cut -d"." -f1`

DELLY=delly
REF=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
EXCL=human.hg38.excl.tsv

# step 1: SV calling 
$DELLY call -g $REF -o $SAMPLE.bcf -x $EXCL $CRAM > log.delly.$SAMPLE.step1.txt 2>&1

