#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V

SAMPLE=$1
REF=$2
VCF=$3
ConvInv=convertInversion.py
SAMTOOLS=samtools

mkdir $SAMPLE
cd $SAMPLE
$ConvInv $SAMTOOLS $REF__ycu $VCF >  diploidSV.inv.vcf
