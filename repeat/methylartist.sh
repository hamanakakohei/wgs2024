REF=~/resource/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF=~/resource/gencode.v44.annotation__CHR.sort.gtf.gz
INTERVAL=chr17:7883000-7889000
RATIO=10,35,10,30,30 

methylartist locus -r $REF -g $GTF -n CG --svg --genes CHD3 \
  -b $BAM \
  -i $INTERVAL \
  -p $RATIO \
  --readopenmarkeredgecolor grey \
  --samplepalette grey \
  --width 16 --height 5.75 


