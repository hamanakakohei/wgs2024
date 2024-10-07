#! /usr/bin/Rscript

library(tidyverse)
FILE1 = commandArgs(trailingOnly=TRUE)[1]
FILE2 = commandArgs(trailingOnly=TRUE)[2]
FILE3 = commandArgs(trailingOnly=TRUE)[3]
FILE4 = commandArgs(trailingOnly=TRUE)[4]
FILE5 = commandArgs(trailingOnly=TRUE)[5]
OUT   = commandArgs(trailingOnly=TRUE)[6]

split_genes_and_join = function(GENES, ANNO__df){
  tmp = tibble(gene=str_split(GENES, ";")[[1]]) %>% left_join(ANNO__df) %>% select(-gene)
  for(i in colnames(tmp)){
    eval(parse(text=paste0(i,' = paste(tmp$',i,',collapse=";")')))
  }
  return( eval(parse(text=paste0('tibble(',paste0(colnames(tmp),collapse=","),')'))) )
}

varid_sample_fa_mo_annos = read_tsv(FILE1,guess_max=Inf) %>% rename(varid=varid__vcf) 
sample_varid_gt	 	 = read_tsv(FILE2,guess_max=Inf, col_names=c('sample','varid','gt'))
gene_annos__ddd   	 = read_tsv(FILE3,guess_max=Inf) %>% rename(gene=gene_symbol) %>% select(c(gene,disease_name,DDD_category,allelic_requirement,mutation_consequence))
gene_annos__omim   	 = read_tsv(FILE4,guess_max=Inf, col_names=c('gene','omim'))
gene_annos__gnomad 	 = read_tsv(FILE5,guess_max=Inf) %>% filter(canonical==TRUE) %>% select(c(gene,pLI,oe_lof,oe_lof_upper,oe_mis,mis_z)) 
left_join(varid_sample_fa_mo_annos, sample_varid_gt) %>% 
  left_join(rename(sample_varid_gt, fa=sample, gt.fa=gt)) %>% 
  left_join(rename(sample_varid_gt, mo=sample, gt.mo=gt)) %>% 
  nest(-Gene.refGene) %>% 
  mutate(ddd   =map(Gene.refGene,~split_genes_and_join(.x,gene_annos__ddd))) %>% 
  mutate(gnomad=map(Gene.refGene,~split_genes_and_join(.x,gene_annos__gnomad))) %>% 
  mutate(omim  =map(Gene.refGene,~split_genes_and_join(.x,gene_annos__omim))) %>% 
  unnest %>%
  write_tsv(OUT)


