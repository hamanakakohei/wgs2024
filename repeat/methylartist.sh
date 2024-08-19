


# DNAme

#REF=~/resource/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
REF=~/resource/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
#REF=~/resource/ref/GRCh38.primary_assembly.genome.fa
GTF=~/resource/gencode.v44.annotation__CHR.sort.gtf.gz
#BAM=/home/kohei.hamanaka/wgs/dname/ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.bam
#BAM=/home/kohei.hamanaka/wgs/dname/ONT_adp_218_5mc_sup_v14.allele1.sorted.bam
#BAM=/home/kohei.hamanaka/wgs/dname/ONT_adp_218_5mc_sup_v14.allele2.sorted.bam
#BAM=/home/kohei.hamanaka/wgs/dname/ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.bam
#BAM=/home/kohei.hamanaka/wgs/dname/ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.bam
#BAM=/home/kohei.hamanaka/wgs/dname/ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.bam
#BAM=ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.bam,ONT_adp_218_5mc_sup_v14.allele1.sorted.bam,ONT_adp_218_5mc_sup_v14.allele2.sorted.bam,ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.bam,ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.bam,ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.bam
BAM=ONT_adp_218_5mc_sup_v14.allele1.sorted.bam,ONT_adp_218_5mc_sup_v14.allele2.sorted.bam,ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.bam,ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.bam
BAM=ONT_adp_218_5mc_sup_v14.allele1.sorted.bam,ONT_adp_218_5mc_sup_v14.allele2.sorted.bam
BAM=ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.bam,ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.bam
INTERVAL=chr17:7880000-7890000
INTERVAL=chr17:7883000-7889000
#HIGHLIGHT=7884000-7886000
HIGHLIGHT=7885308-7885345
RATIO=10,35,10,30,30 # hight 5.75, BAM mo; 8.9, BAM pt
#RATIO=10,57,10,30,30 

methylartist locus -r $REF -g $GTF -n CG --svg --genes CHD3 \
  -b $BAM \
  -i $INTERVAL \
  -p $RATIO \
  --readopenmarkeredgecolor grey \
  --samplepalette grey \
  --width 16 --height 5.75 # 6.85
  #-l $HIGHLIGHT \


# DNAme within repeat expansion
samtools  fastq  ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.bam  > ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.fq                   
samtools  fastq  ONT_adp_218_5mc_sup_v14.allele1.sorted.bam         > ONT_adp_218_5mc_sup_v14.allele1.sorted.fq                   
samtools  fastq  ONT_adp_218_5mc_sup_v14.allele2.sorted.bam         > ONT_adp_218_5mc_sup_v14.allele2.sorted.fq                   
samtools  fastq  ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.bam    > ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.fq                   
samtools  fastq  ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.bam    > ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.fq                   
samtools  fastq  ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.bam  > ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.fq                   

echo -e 'chr17\t7880000\t7890000' > region_oi.bed
bedtools getfasta -fi $REF -bed region_oi.bed -fo chd3.fasta

GGCCTCCTCCCCTCCCCCGCCGCACCCCTCCCCC no atoni (CCG)250/160 sounyuu

REF250=chd3.ccg250.fasta
REF160=chd3.ccg160.fasta
minimap2 -a $REF250 ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.fq > ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.sam
minimap2 -a $REF250 ONT_adp_218_5mc_sup_v14.allele1.sorted.fq        > ONT_adp_218_5mc_sup_v14.allele1.sorted.minimap2.sam
minimap2 -a $REF250 ONT_adp_218_5mc_sup_v14.allele2.sorted.fq        > ONT_adp_218_5mc_sup_v14.allele2.sorted.minimap2.sam
minimap2 -a $REF160 ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.fq   > ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.minimap2.sam
minimap2 -a $REF160 ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.fq   > ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.minimap2.sam
minimap2 -a $REF160 ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.fq > ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.sam

samtools sort -O bam -o ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.bam ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.sam
samtools sort -O bam -o ONT_adp_218_5mc_sup_v14.allele1.sorted.minimap2.bam        ONT_adp_218_5mc_sup_v14.allele1.sorted.minimap2.sam
samtools sort -O bam -o ONT_adp_218_5mc_sup_v14.allele2.sorted.minimap2.bam        ONT_adp_218_5mc_sup_v14.allele2.sorted.minimap2.sam
samtools sort -O bam -o ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.minimap2.bam   ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.minimap2.sam
samtools sort -O bam -o ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.minimap2.bam   ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.minimap2.sam
samtools sort -O bam -o ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.bam ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.sam

samtools index ONT_adp_218_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.bam 
samtools index ONT_adp_218_5mc_sup_v14.allele1.sorted.minimap2.bam        
samtools index ONT_adp_218_5mc_sup_v14.allele2.sorted.minimap2.bam        
samtools index ONT_adp_219_5khz_5mc_sup_v14.allele1.sorted.minimap2.bam   
samtools index ONT_adp_219_5khz_5mc_sup_v14.allele2.sorted.minimap2.bam   
samtools index ONT_adp_219_5khz_5mc_sup_v14.sorted.oneRegion.minimap2.bam 
