#! /usr/bin/Rscript

library(tidyverse)

FILE1 = commandArgs(trailingOnly=TRUE)[1]
FILE2 = commandArgs(trailingOnly=TRUE)[2]
FILE3 = commandArgs(trailingOnly=TRUE)[3]
FILE4 = commandArgs(trailingOnly=TRUE)[4]
FILE5 = commandArgs(trailingOnly=TRUE)[5]
FILE6 = commandArgs(trailingOnly=TRUE)[6]
FILE7 = commandArgs(trailingOnly=TRUE)[7]
OUT = commandArgs(trailingOnly=TRUE)[8]


split_genes_and_join = function(GENES, ANNO__df){
  tmp = tibble(gene=str_split(GENES, ";")[[1]]) %>% left_join(ANNO__df) %>% select(-gene)
  for(i in colnames(tmp)){
    eval(parse(text=paste0(i,' = paste(tmp$',i,',collapse=";")')))
  }
  return( eval(parse(text=paste0('tibble(',paste0(colnames(tmp),collapse=","),')'))) )
}

varid_proband            = read_delim(FILE1,delim=' ',col_names=c('varid','sample')) %>% mutate(varid=str_replace_all(varid,'-','_'))
fam_sample		 = read_tsv(FILE2,col_names=c('fam','sample'))
sample_varid_gt          = read_tsv(FILE3, col_names=c('sample','varid','gt'))
varid_AnnovarAnnos	 = read_tsv(FILE4,guess_max=Inf) %>% unite(varid, c(chr__vcf,pos__vcf,ref__vcf,alt__vcf), sep="_") %>% select(-DA0000000348)
gene_annos__ddd          = read_tsv(FILE5) %>% rename(gene=gene_symbol) %>% select(c(gene,disease_name,DDD_category,allelic_requirement,mutation_consequence))
gene_annos__omim         = read_tsv(FILE6, col_names=c('gene','omim'))
gene_annos__gnomad       = read_tsv(FILE7) %>% filter(canonical==TRUE) %>% select(c(gene,pLI,oe_lof,oe_lof_upper,oe_mis,mis_z))

sample_membs = fam_sample %>% nest(-fam) %>% mutate(membs=map_chr(data,~paste(.x %>% arrange(sample) %>% pull(sample), collapse=';'))) %>% unnest %>% select(fam,membs) %>% rename(sample=fam) %>% distinct
sample_varid_gts = sample_varid_gt %>% inner_join(fam_sample) %>% nest(-c(fam,varid)) %>% mutate(gts=map_chr(data,~paste(.x %>% arrange(sample) %>% pull(gt), collapse=';'))) %>% unnest %>% select(fam,varid,gts) %>% rename(sample=fam) %>% distinct

varid_proband %>%
  left_join(sample_membs) %>% 
  left_join(sample_varid_gts) %>% 
  left_join(varid_AnnovarAnnos) %>% 
  nest(-Gene.refGene) %>%
  mutate(ddd   =map(Gene.refGene,~split_genes_and_join(.x,gene_annos__ddd))) %>%
  mutate(gnomad=map(Gene.refGene,~split_genes_and_join(.x,gene_annos__gnomad))) %>%
  mutate(omim  =map(Gene.refGene,~split_genes_and_join(.x,gene_annos__omim))) %>%
  unnest %>%
  write_tsv(OUT)
