#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V

VCF=$1
REF=$2
OUT=$3
D=$4 #default 50; max=4999

source python2-venv/env.sh
spliceai -I $VCF -O $OUT -R $REF -A grch38 -D $D

