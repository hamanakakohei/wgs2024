#! /usr/bin/Rscript

library(tidyverse)
FILE1 = commandArgs(trailingOnly=TRUE)[1]
FILE2 = commandArgs(trailingOnly=TRUE)[2]
OUT   = commandArgs(trailingOnly=TRUE)[3]

varid__vcf_DenovoAnnos    = read_tsv(FILE1,guess_max=Inf) %>% unite(varid__vcf, c(chr,pos,ref,alt), sep="_")
varid__vcf_AnnovarAnnos	  = read_tsv(FILE2,guess_max=Inf) %>% unite(varid__vcf, c(chr__vcf,pos__vcf,ref__vcf,alt__vcf), sep="_") %>% select(-DA0000000348)
inner_join(varid__vcf_DenovoAnnos, varid__vcf_AnnovarAnnos) %>% write_tsv(OUT)

