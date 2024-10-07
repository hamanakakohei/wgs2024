#!/usr/bin/env python
import argparse
import pandas as pd

def main():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument("-annovar" , type = str, help="")
  parser.add_argument("-spliceai", type = str, help="")
  parser.add_argument("-topmed"  , type = str, help="")
  parser.add_argument("-tommo"   , type = str, help="")
  parser.add_argument("-wbbc"    , type = str, help="")
  parser.add_argument("-ncbn"    , type = str, help="")
  parser.add_argument("-ncgm"    , type = str, help="")
  parser.add_argument("-ycu"     , type = str, help="")
  parser.add_argument("-out"     , type = str, help="") 
  args = parser.parse_args()
  annovar  = pd.read_table(args.annovar ,                                          low_memory=False) 
  spliceai = pd.read_table(args.spliceai,  names=('VarId','spliceai')             ,low_memory=False) 
  topmed   = pd.read_table(args.topmed  ,  names=('VarId','topmed_q','topmed_af') ,low_memory=False) 
  tommo    = pd.read_table(args.tommo   ,  names=('VarId','tommo_q','tommo_af')   ,low_memory=False) 
  wbbc     = pd.read_table(args.wbbc    ,  names=('VarId','wbbc_q','wbbc_af')     ,low_memory=False) 
  ncbn     = pd.read_table(args.ncbn    ,  names=('VarId','ncbn_q','ncbn_all_af','ncbn_hondo_af','ncbn_ryukyu_af')     ,low_memory=False) 
  ncgm     = pd.read_table(args.ncgm    ,  names=('VarId','ncgm_af')              ,low_memory=False) 
  ycu      = pd.read_table(args.ycu     ,  names=('VarId','ycu_af')               ,low_memory=False) 
  out      = args.out
  annovar['VarId'] = annovar['chr__vcf']+'-'+annovar['pos__vcf'].astype(str)+'-'+annovar['ref__vcf']+'-'+annovar['alt__vcf']
  annovar = pd.merge(annovar,spliceai,how='left')
  annovar = pd.merge(annovar,topmed  ,how='left')
  annovar = pd.merge(annovar,tommo   ,how='left')
  annovar = pd.merge(annovar,wbbc    ,how='left')
  annovar = pd.merge(annovar,ncbn    ,how='left')
  annovar = pd.merge(annovar,ncgm    ,how='left')
  annovar = pd.merge(annovar,ycu     ,how='left')
  annovar = annovar.fillna({'topmed_af':0,'tommo_af':0,'wbbc_af':0,'ncbn_all_af':0,'ncbn_hondo_af':0,'ncbn_ryukyu_af':0,'ncgm_af':0,'ycu_af':0}) 
  annovar.to_csv(out,sep='\t',index=False)

if __name__ == "__main__":
  main()
