#!/usr/bin/env python
import argparse
import vcf
import pandas as pd


def main():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument("-annovar"         , type=str, help="") 
  parser.add_argument('-esp6500siv2_all' , type=float, default=0.0002,help="") 
  parser.add_argument('-gnomAD_exome_ALL', type=float, default=0.00001,help="") 
  parser.add_argument('-gnomAD_exome_AFR', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_AMR', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_ASJ', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_EAS', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_FIN', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_NFE', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_OTH', type=float, default=1,help="") 
  parser.add_argument('-gnomAD_exome_SAS', type=float, default=1,help="") 
  parser.add_argument('-AF'              , type=float, default=0.00002,help="") 
  parser.add_argument('-AF_afr'          , type=float, default=1,help="") 
  parser.add_argument('-AF_ami'          , type=float, default=1,help="") 
  parser.add_argument('-AF_amr'          , type=float, default=1,help="") 
  parser.add_argument('-AF_asj'          , type=float, default=1,help="") 
  parser.add_argument('-AF_eas'          , type=float, default=1,help="") 
  parser.add_argument('-AF_fin'          , type=float, default=1,help="") 
  parser.add_argument('-AF_nfe'          , type=float, default=1,help="") 
  parser.add_argument('-AF_oth'          , type=float, default=1,help="") 
  parser.add_argument('-AF_sas'          , type=float, default=1,help="") 
  parser.add_argument("-topmed"          , type=float, default=0.00002,help="")
  parser.add_argument("-tommo"           , type=float, default=0.0001,help="")
  parser.add_argument("-wbbc"            , type=float, default=0.0003,help="")
  parser.add_argument("-ncbn_all"        , type=float, default=0.0001,help="")
  parser.add_argument("-ncbn_hondo"      , type=float, default=1,help="")
  parser.add_argument("-ncbn_ryukyu"     , type=float, default=1,help="")
  parser.add_argument("-ncgm"            , type=float, default=0.001,help="")
  parser.add_argument("-out"             , type=str  , help="") 
  args = parser.parse_args()
  annovar = pd.read_table(args.annovar ,low_memory=False)
  out     = args.out          
  esp6500siv2_all  = args.esp6500siv2_all 
  gnomAD_exome_ALL = args.gnomAD_exome_ALL
  gnomAD_exome_AFR = args.gnomAD_exome_AFR
  gnomAD_exome_AMR = args.gnomAD_exome_AMR
  gnomAD_exome_ASJ = args.gnomAD_exome_ASJ
  gnomAD_exome_EAS = args.gnomAD_exome_EAS
  gnomAD_exome_FIN = args.gnomAD_exome_FIN
  gnomAD_exome_NFE = args.gnomAD_exome_NFE
  gnomAD_exome_OTH = args.gnomAD_exome_OTH
  gnomAD_exome_SAS = args.gnomAD_exome_SAS
  AF               = args.AF              
  AF_afr           = args.AF_afr          
  AF_ami           = args.AF_ami          
  AF_amr           = args.AF_amr          
  AF_asj           = args.AF_asj          
  AF_eas           = args.AF_eas          
  AF_fin           = args.AF_fin          
  AF_nfe           = args.AF_nfe          
  AF_oth           = args.AF_oth          
  AF_sas           = args.AF_sas          
  topmed           = args.topmed          
  tommo            = args.tommo           
  wbbc             = args.wbbc            
  ncbn_all         = args.ncbn_all        
  ncbn_hondo       = args.ncbn_hondo      
  ncbn_ryukyu      = args.ncbn_ryukyu     
  ncgm             = args.ncgm            
  annovar = annovar.\
    astype({'esp6500siv2_all' :str}).replace({'esp6500siv2_all' :{'.':'0'}}).astype({'esp6500siv2_all' :float}).\
    astype({'gnomAD_exome_ALL':str}).replace({'gnomAD_exome_ALL':{'.':'0'}}).astype({'gnomAD_exome_ALL':float}).\
    astype({'gnomAD_exome_AFR':str}).replace({'gnomAD_exome_AFR':{'.':'0'}}).astype({'gnomAD_exome_AFR':float}).\
    astype({'gnomAD_exome_AMR':str}).replace({'gnomAD_exome_AMR':{'.':'0'}}).astype({'gnomAD_exome_AMR':float}).\
    astype({'gnomAD_exome_ASJ':str}).replace({'gnomAD_exome_ASJ':{'.':'0'}}).astype({'gnomAD_exome_ASJ':float}).\
    astype({'gnomAD_exome_EAS':str}).replace({'gnomAD_exome_EAS':{'.':'0'}}).astype({'gnomAD_exome_EAS':float}).\
    astype({'gnomAD_exome_FIN':str}).replace({'gnomAD_exome_FIN':{'.':'0'}}).astype({'gnomAD_exome_FIN':float}).\
    astype({'gnomAD_exome_NFE':str}).replace({'gnomAD_exome_NFE':{'.':'0'}}).astype({'gnomAD_exome_NFE':float}).\
    astype({'gnomAD_exome_OTH':str}).replace({'gnomAD_exome_OTH':{'.':'0'}}).astype({'gnomAD_exome_OTH':float}).\
    astype({'gnomAD_exome_SAS':str}).replace({'gnomAD_exome_SAS':{'.':'0'}}).astype({'gnomAD_exome_SAS':float}).\
    astype({'AF'              :str}).replace({'AF'              :{'.':'0'}}).astype({'AF'              :float}).\
    astype({'AF_afr'          :str}).replace({'AF_afr'          :{'.':'0'}}).astype({'AF_afr'          :float}).\
    astype({'AF_ami'          :str}).replace({'AF_ami'          :{'.':'0'}}).astype({'AF_ami'          :float}).\
    astype({'AF_amr'          :str}).replace({'AF_amr'          :{'.':'0'}}).astype({'AF_amr'          :float}).\
    astype({'AF_asj'          :str}).replace({'AF_asj'          :{'.':'0'}}).astype({'AF_asj'          :float}).\
    astype({'AF_eas'          :str}).replace({'AF_eas'          :{'.':'0'}}).astype({'AF_eas'          :float}).\
    astype({'AF_fin'          :str}).replace({'AF_fin'          :{'.':'0'}}).astype({'AF_fin'          :float}).\
    astype({'AF_nfe'          :str}).replace({'AF_nfe'          :{'.':'0'}}).astype({'AF_nfe'          :float}).\
    astype({'AF_oth'          :str}).replace({'AF_oth'          :{'.':'0'}}).astype({'AF_oth'          :float}).\
    astype({'AF_sas'          :str}).replace({'AF_sas'          :{'.':'0'}}).astype({'AF_sas'          :float}).\
    astype({'topmed_af'       :str}).replace({'topmed_af'       :{'.':'0'}}).astype({'topmed_af'       :float}).\
    astype({'tommo_af'        :str}).replace({'tommo_af'        :{'.':'0'}}).astype({'tommo_af'        :float}).\
    astype({'wbbc_af'         :str}).replace({'wbbc_af'         :{'.':'0'}}).astype({'wbbc_af'         :float}).\
    astype({'ncbn_all_af'     :str}).replace({'ncbn_all_af'     :{'.':'0'}}).astype({'ncbn_all_af'     :float}).\
    astype({'ncbn_hondo_af'   :str}).replace({'ncbn_hondo_af'   :{'.':'0'}}).astype({'ncbn_hondo_af'   :float}).\
    astype({'ncbn_ryukyu_af'  :str}).replace({'ncbn_ryukyu_af'  :{'.':'0'}}).astype({'ncbn_ryukyu_af'  :float}).\
    astype({'ncgm_af'         :str}).replace({'ncgm_af'         :{'.':'0'}}).astype({'ncgm_af'         :float}).\
    query(" \
    esp6500siv2_all  <  @esp6500siv2_all   and \
    gnomAD_exome_ALL <  @gnomAD_exome_ALL  and \
    gnomAD_exome_AFR <  @gnomAD_exome_AFR  and \
    gnomAD_exome_AMR <  @gnomAD_exome_AMR  and \
    gnomAD_exome_ASJ <  @gnomAD_exome_ASJ  and \
    gnomAD_exome_EAS <  @gnomAD_exome_EAS  and \
    gnomAD_exome_FIN <  @gnomAD_exome_FIN  and \
    gnomAD_exome_NFE <  @gnomAD_exome_NFE  and \
    gnomAD_exome_OTH <  @gnomAD_exome_OTH  and \
    gnomAD_exome_SAS <  @gnomAD_exome_SAS  and \
    AF               <  @AF                and \
    AF_afr           <  @AF_afr            and \
    AF_ami           <  @AF_ami            and \
    AF_amr           <  @AF_amr            and \
    AF_asj           <  @AF_asj            and \
    AF_eas           <  @AF_eas            and \
    AF_fin           <  @AF_fin            and \
    AF_nfe           <  @AF_nfe            and \
    AF_oth           <  @AF_oth            and \
    AF_sas           <  @AF_sas            and \
    (topmed_q!='PASS' or topmed_q!=topmed_q or topmed_af <  @topmed ) and \
    ( tommo_q!='PASS' or  tommo_q!=tommo_q  or tommo_af  <  @tommo  ) and \
    (  wbbc_q!='PASS' or   wbbc_q!=wbbc_q   or wbbc_af   <  @wbbc   ) and \
    (  ncbn_q!='PASS' or   ncbn_q!=ncbn_q   or (ncbn_all_af < @ncbn_all and ncbn_hondo_af < @ncbn_hondo and ncbn_ryukyu_af < @ncbn_ryukyu )) and \
    ncgm_af          <  @ncgm              ")
  annovar['VarId'].drop_duplicates().to_csv(out,index=False,header=False)

if __name__ == "__main__":
  main()
