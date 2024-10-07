#!/usr/bin/env python

import argparse
import vcf
import pandas as pd
import numpy as np

def parse_spai(df):
  df['spaiscore'] = df['spliceai'].map(lambda x: spai_to_score(x))
  return df

def spai_to_score(spais):
  scores = []
  if pd.isna(spais):
    return np.nan
  else:
    for spai in spais.split(';')[:-1]:
      scores = scores + list(map(float,spai.split('|')[2:6]))
    return max(scores)

def main():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument("-annovar", type=str, help="") 
  parser.add_argument("-out",     type=str, help="") 
  args = parser.parse_args()
  out     = args.out          
  snv    = pd.read_table(args.annovar,low_memory=False) 
  snv = parse_spai(snv)
  snv.replace({'CADD_phred'      :{'.':'0'}}).astype({'CADD_phred'      :float}).\
    query(" \
      ( \
        `ExonicFunc.refGene`=='frameshift deletion' or \
        `ExonicFunc.refGene`=='frameshift insertion' or \
        `ExonicFunc.refGene`=='nonframeshift deletion' or \
        `ExonicFunc.refGene`=='nonframeshift insertion' or \
        `ExonicFunc.refGene`=='stopgain' or \
        `ExonicFunc.refGene`=='stoploss' or \
        'athogenic' in CLNSIG) or \
      ( \
        `ExonicFunc.refGene`=='nonsynonymous SNV' and \
        CADD_phred > 20 \
      ) or \
      spaiscore > 0.15 \
    ")['VarId'].drop_duplicates().to_csv(out,index=False,header=False)

if __name__ == "__main__":
  main()
