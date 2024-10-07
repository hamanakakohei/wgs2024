#! /usr/bin/Rscript

library(tidyverse)
FILE1 = commandArgs(trailingOnly=TRUE)[1]
FILE2 = commandArgs(trailingOnly=TRUE)[2]
OUTPUT = commandArgs(trailingOnly=TRUE)[3]

chr_pos_ref_alt_sample_fa_mo = read_tsv(FILE1,col_types=cols(.default="c"),col_names=c("id","sample","fa","mo")) %>% separate(id,c("chr","pos","ref","alt"),sep="-")
chr_pos_ref_alt_sp    	     = read_tsv(FILE2,col_types=cols(.default="c"),col_names=c("chr","pos","ref","alt","sp"))

chr_pos_ref_alt_sample_fa_mo %>%
    mutate(type=if_else(str_length(ref)==1 & str_length(alt)==1,"snv","indel")) %>%
    left_join(chr_pos_ref_alt_sp) %>%
    rename(spliceai_manual=sp) %>%
    write_tsv(OUTPUT)
