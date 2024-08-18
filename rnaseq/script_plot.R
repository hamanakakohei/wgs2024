
source('C:/Users/hamanakakohei/Dropbox/backup20230402/_antares01_analysis_hamanaka_function_ggplot_setting.R')
#source('/Users/hamanaka/Dropbox/backup20230402/_antares01_analysis_hamanaka_function_ggplot_setting.R')

library(tidyverse)
#DATA = 'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/sample27.tmm-norm_count.txt'
DATA = 'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/sample23.tmm-norm_count.expression.txt.gz'
#DATA = '/Users/hamanaka/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/sample23.tmm-norm_count.expression.txt.gz'
BAD_SAMPLES = c(
  #'100222'
   '140698'
  ,'180767'
  ,'182848'
  ,'200132'
  ,'201535'
  #,'210488'
  #,'221122'
)
BAD_SAMPLES = c( 'U173381' ) #'H120109'
HEALTHYS = c( 'H110283', 'H120109', 'H131696', 'H133133', 'H133326', 'H141189', 'H144750' )


#ENSG = 'ENSG00000173575.20' # CHD2
#ENSG = 'ENSG00000272888.5' # CHASERR
#ENSG = 'ENSG00000177565.16' # TBL1XR1
#ENSG = 'ENSG00000140992.18' # PDPK1
#ENSG = 'ENSG00000185883.11' # ATP6V0C
# pdpk1周辺
ENSGS = c(
  # "ENSG00000103024.7"   # NME3"
  #,"ENSG00000074071.14"  # MRPS34"
  #,"ENSG00000197774.12"  # EME2"
  #,"ENSG00000162032.15"  # SPSB3"
  #,"ENSG00000095906.16"  # NUBP2"
  #,"ENSG00000099769.5"   # IGFALS"
  #,"ENSG00000063854.12"  # HAGH"
  #,"ENSG00000180185.11"  # FAHD1"
  #,"ENSG00000162039.14"  # MEIOB"
  #,"ENSG00000260541.1"   # LA16c-429E7.1"
  #,"ENSG00000281219.2"   # LINC00254"
  #,"ENSG00000259947.1"   # LINC02124"
  #,"ENSG00000162040.6"   # HS3ST6"
  #,"ENSG00000198736.11"  # MSRB1"
  #,"ENSG00000140986.7"   # RPL3L"
  #,"ENSG00000140990.14"  # NDUFB10"
  #,"ENSG00000140988.15"  # RPS2"
  #,"ENSG00000206811.1"   # SNORA10"
  #,"ENSG00000207405.1"   # SNORA64"
  #,"ENSG00000255513.1"   # AC005363.9"
  #,"ENSG00000255198.4"   # SNHG9"
  #,"ENSG00000273587.1"   # SNORA78"
  #,"ENSG00000179580.9"   # RNF151"
  #,"ENSG00000277602.1"   # AC005363.11"
  #,"ENSG00000183751.14"  # TBL3"
  #,"ENSG00000196408.11"  # NOXO1"
  #,"ENSG00000127554.12"  # GFER"
  #,"ENSG00000261790.1"   # AC005606.14"
  #,"ENSG00000127561.14"  # SYNGR3"
  #,"ENSG00000260107.1"   # AC005606.15"
  #,"ENSG00000167962.13"  # ZNF598"
  #,'ENSG00000183971.6' # NPW
  #,'ENSG00000065054.13' # SLC9A3R2
  #,'ENSG00000065057.7' # NTHL1
  #,'ENSG00000103197.16' # TSC2
  #,'ENSG00000008710.19' # PKD1
  #,'ENSG00000167964.12' # RAB26
  #,'ENSG00000131653.12' # TRAF7
  #,'ENSG00000167971.15' # CASKIN1
   'ENSG00000167965.17' # MLST8
  ,'ENSG00000182685.7' # BRICD5
  ,'ENSG00000184207.8' # PGP
  ,'ENSG00000167967.15' # E4F1
  #,'ENSG00000167968.12' # DNASE1L2
  ,'ENSG00000167969.12' # ECI1
  ,'ENSG00000205937.11' # RNPS1
  ,'ENSG00000167972.13' # ABCA3       1.5? aaaa
  ,'ENSG00000162063.12' # CCNF        aaaaa
  ,'ENSG00000162062.14' # C16orf59    aaaa
  #,'ENSG00000162068.1' # NTN3         NA
  ,'ENSG00000162065.12' # TBC1D24      aaa
  ,'ENSG00000185883.11' # ATP6V0C     aaaa
  ,'ENSG00000162066.14' # AMDHD2    aaaa
  ,'ENSG00000140992.18' # PDPK1    aaaa
  ,'ENSG00000167977.8' # KCTD5     1.5? aaaaa
  #,'ENSG00000172382.9' # PRSS27    NA
  ,'ENSG00000167978.16' # SRRM2
  ,'ENSG00000103363.14' # ELOB
  #,'ENSG00000103355.13' # PRSS33
  #,'ENSG00000215148.8' # PRSS41
  #,'ENSG00000007038.10' # PRSS21
  #,'ENSG00000162078.11' # ZG16B
  #,'ENSG00000005001.9' # PRSS22
  #,'ENSG00000162076.12' # FLYWCH2
  #,'ENSG00000059122.16' # FLYWCH1
  #,'ENSG00000131650.13' # KREMEN2
  #,'ENSG00000127564.16' # PKMYT1
  #,'ENSG00000162073.13' # PAQR4
  #,"ENSG00000274367.1"  # LA16c-380H5.6"
  #,"ENSG00000270168.2"  # LA16c-380H5.3"
  #,"ENSG00000262152.6"  # LINC00514"
  #,"ENSG00000262362.1"  # LA16c-380H5.2"
  #,"ENSG00000272079.2"  # LA16c-380H5.5"
  #,"ENSG00000213937.3"  # CLDN9"
  #,"ENSG00000184697.6"  # CLDN6"
  #,"ENSG00000006327.13" # TNFRSF12A"
  #,"ENSG00000103145.10" # HCFC1R1"
  #,"ENSG00000131652.13" # THOC6"
  #,"ENSG00000162069.14" # BICDL2"
  #,"ENSG00000205890.3"  # RP11-473M20.5"
  #,"ENSG00000008516.16" # MMP25"
  #,"ENSG00000261971.6"  # MMP25-AS1"
  #,"ENSG00000008517.16" # IL32"
  #,"ENSG00000252561.1"  # RNU1-125P"
  #,"ENSG00000262370.5"  # RP11-473M20.9"
  #,"ENSG00000200204.1"  # RNU1-22P"
  #,"ENSG00000130182.7"  # ZSCAN10"
  #,"ENSG00000263011.1"  # RP11-473M20.11"
  #,"ENSG00000263072.6"  # ZNF213-AS1"
  #,"ENSG00000122386.10" # ZNF205"
  #,"ENSG00000085644.13" # ZNF213"
  #,"ENSG00000228146.7"  # CASP16P"
  #,"ENSG00000261889.1"  # RP11-473M20.16"
  #,"ENSG00000262521.1"  # AJ003147.8"
  #,"ENSG00000262668.1"  # AJ003147.9"
  #,"ENSG00000168124.2"  # OR1F1"
  #,"ENSG00000203581.7"  # OR1F2P"
  #,"ENSG00000010539.11" # ZNF200"
  #,"ENSG00000103313.12" # MEFV"
  #,"ENSG00000281005.1"  # LINC00921"
  #,"ENSG00000006194.9"  # ZNF263"
  #,"ENSG00000279330.1"  # AJ003147.11"
  #,"ENSG00000279031.1"  # LA16c-360H6.1"
  #,"ENSG00000140993.10" # TIGD7"
  #,"ENSG00000162086.14" # ZNF75A"
  #,"ENSG00000262899.1"  # LA16c-360H6.3"
  #,"ENSG00000262554.1"  # LA16c-360H6.2"
  #,"ENSG00000262316.1"  # RP11-433P17.1"
  #,"ENSG00000168158.3"  # OR2C1"
  #,"ENSG00000262621.4"  # LA16c-306E5.2"
  #,"ENSG00000262118.1"  # MTCO1P28"
  #,"ENSG00000263177.1"  # MTND1P8"
  #,"ENSG00000232196.3"  # MTRNR2L4"
  #,"ENSG00000140987.19" # ZSCAN32"
  #,"ENSG00000103343.12" # ZNF174"
  #,"ENSG00000167981.5"  # ZNF597"
  #,"ENSG00000122390.17" # NAA60"
  #,"ENSG00000263212.2"  # LA16c-306E5.3"
)

# tbl1xr1周辺
ENSGS = c(
   "ENSG00000201648.1"  # "RNU4-91P"
  ,"ENSG00000236241.1"  # "RP11-744O11.2"
  ,"ENSG00000214192.3"  # "UBE2V1P2"
  ,"ENSG00000201452.1"  # "RNU6-1317P"
  ,"ENSG00000225552.1"  # "NAALADL2-AS1"
  ,"ENSG00000260217.1"  # "RP11-809F4.3"
  ,"ENSG00000230362.1"  # "ACTG1P23"
  ,"ENSG00000236357.1"  # "EI24P1"
  ,"ENSG00000271413.1"  # "RP11-809F4.4"
  ,"ENSG00000238026.1"  # "RP11-78E6.1"
  ,"ENSG00000283333.1"  # "MIR7977"
  ,"ENSG00000223715.1"  # "LINC01208"
  ,"ENSG00000232461.1"  # "RP11-644C3.1"
  ,"ENSG00000252302.1"  # "RNA5SP147"
  ,"ENSG00000228308.1"  # "LINC01209"
  ,"ENSG00000231888.1"  # "MTND5P15"
  ,"ENSG00000177565.16" # "TBL1XR1"
  ,"ENSG00000231310.3"  # "TBL1XR1-AS1"
  ,"ENSG00000201343.1"  # "Y_RNA"
  ,"ENSG00000200882.1"  # "RNU6-681P"
  ,"ENSG00000203645.2"  # "LINC00501"
  ,"ENSG00000236893.1"  # "ASS1P7"
  ,"ENSG00000226782.2"  # "RP11-706D8.3"
  ,"ENSG00000228221.5"  # "LINC00578"
  ,"ENSG00000252028.1"  # "RN7SKP52"
  ,"ENSG00000236385.1"  # "RP11-114M1.2"
  ,"ENSG00000200288.1"  # "SNORA18"
  ,"ENSG00000228561.2"  # "RP11-114M1.1"
  ,"ENSG00000277241.1"  # "RP11-114M1.3"
  ,"ENSG00000231574.5"  # "LINC02015"
  ,"ENSG00000225790.1"  # "RP11-2L8.1"
  ,"ENSG00000270321.1"  # "RP11-2L8.2"
  ,"ENSG00000199858.1"  # "RNU6-1120P"
  ,"ENSG00000223930.5"  # "RP11-33A14.1"
  ,"ENSG00000197584.11" # "KCNMB2"
  ,"ENSG00000223941.5"  # "LINC01014"
  ,"ENSG00000252336.1"  # "RNA5SP148"
  ,"ENSG00000271131.1"  # "RP11-33A14.3"
  ,"ENSG00000237978.5"  # "KCNMB2-AS1"
  ,"ENSG00000172667.10" # "ZMAT3"
  ,"ENSG00000200616.1"  # "Y_RNA"
  ,"ENSG00000229102.1"  # "RP11-360P21.2"
  ,"ENSG00000121879.3"  # "PIK3CA"
  ,"ENSG00000201957.1"  # "SNORA25"
  ,"ENSG00000171121.16" # "KCNMB3"
  ,"ENSG00000240429.1"  # "LRRFIP1P1"
  ,"ENSG00000121864.9"  # "ZNF639"
  ,"ENSG00000260743.1"  # "RP11-255C15.3"
  ,"ENSG00000171109.18" # "MFN1"
  ,"ENSG00000114450.9"  # "GNB4"
  ,"ENSG00000272699.1"  # "RP11-145M9.5"
  ,"ENSG00000239255.1"  # "RP11-145M9.2"
  ,"ENSG00000181260.8"  # "MTHFD2P7"
)

# ogdhl周辺
ENSGS = c(
   "ENSG00000204172.12" # "AGAP9"
  ,"ENSG00000252877.1"  # "RNA5SP312"
  ,"ENSG00000251079.6"  # "BMS1P2"
  ,"ENSG00000265630.2"  # "GLUD1P8"
  ,"ENSG00000215065.3"  # "DUSP8P4"
  ,"ENSG00000266217.2"  # "CTSLP2"
  ,"ENSG00000264404.2"  # "RP11-342C24.8"
  ,"ENSG00000189014.7"  # "FAM35DP"
  ,"ENSG00000276544.1"  # "RN7SL453P"
  ,"ENSG00000226964.1"  # "RHEBP2"
  ,"ENSG00000277758.4"  # "ABC7-42404400C24.1"
  ,"ENSG00000276850.4"  # "CH17-360D5.2"
  ,"ENSG00000264717.5"  # "CH17-360D5.1"
  ,"ENSG00000273760.1"  # "CH17-360D5.3"
  ,"ENSG00000276430.2"  # "FAM25C"
  ,"ENSG00000265018.6"  # "AGAP12P"
  ,"ENSG00000252149.1"  # "RNA5SP315"
  ,"ENSG00000270025.2"  # "BMS1P7"
  ,"ENSG00000278561.1"  # "PTPN20CP"
  ,"ENSG00000170324.20" # "FRMPD2"
  ,"ENSG00000225482.1"  # "RPS6P14"
  ,"ENSG00000107643.15" # "MAPK8"
  ,"ENSG00000279822.1"  # "RP11-541M12.6"
  ,"ENSG00000128805.14" # "ARHGAP22"
  ,"ENSG00000248682.1"  # "ARHGAP22-IT1"
  ,"ENSG00000231906.1"  # "RP11-541M12.3"
  ,"ENSG00000251413.1"  # "RP11-534L6.5"
  ,"ENSG00000236800.1"  # "RP11-534L6.2"
  ,"ENSG00000128815.19" # "WDFY4"
  ,"ENSG00000228754.1"  # "RPL13AP19"
  ,"ENSG00000228403.1"  # "RP11-563N6.6"
  ,"ENSG00000241577.1"  # "RP11-523O18.7"
  ,"ENSG00000165383.11" # "LRRC18"
  ,"ENSG00000233665.8"  # "RP11-523O18.5"
  ,"ENSG00000226576.1"  # "RP11-523O18.1"
  ,"ENSG00000264800.1"  # "MIR4294"
  ,"ENSG00000165633.12" # "VSTM4"
  ,"ENSG00000234736.5"  # "FAM170B-AS1"
  ,"ENSG00000172538.6"  # "FAM170B"
  ,"ENSG00000204161.13" # "C10orf128"
  ,"ENSG00000236208.1"  # "C10orf71-AS1"
  ,"ENSG00000177354.11" # "C10orf71"
  ,"ENSG00000165606.8"  # "DRGX"
  ,"ENSG00000235939.1"  # "RP11-123B3.2"
  ,"ENSG00000225830.12" # "ERCC6"
  ,"ENSG00000271237.1"  # "HMGB1P50"
  ,"ENSG00000070748.18" # "CHAT"
  ,"ENSG00000187714.6"  # "SLC18A3"
  ,"ENSG00000178645.12" # "C10orf53"
  ,"ENSG00000197444.9"  # "OGDHL"
  ,"ENSG00000226389.3"  # "MAPK6PS6"
  ,"ENSG00000229870.1"  # "RPL21P89"
  ,"ENSG00000227345.8"  # "PARG"
  ,"ENSG00000230166.1"  # "RPL35AP24"
  ,"ENSG00000204152.10" # "TIMM23B"
  ,"ENSG00000223111.1"  # "SNORA74"
  ,"ENSG00000178440.5"  # "LINC00843"
  ,"ENSG00000222108.1"  # "RNA5SP317"
  ,"ENSG00000204149.10" # "AGAP6"
  ,"ENSG00000235618.7"  # "FAM21EP"
  ,"ENSG00000226631.1"  # "SLC9A3P3"
  ,"ENSG00000099290.15" # "WASHC2A"
  ,"ENSG00000233011.1"  # "SLC9A3P1"
  ,"ENSG00000188611.14" # "ASAH2"
  ,"ENSG00000225137.1"  # "DYNC1I2P1"
  ,"ENSG00000198964.13" # "SGMS1"
  ,"ENSG00000279863.1"  # "RP11-521C22.2"
  ,"ENSG00000225303.2"  # "RP11-512N4.2"
  ,"ENSG00000238523.1"  # "RNU7-107P"
  ,"ENSG00000226200.6"  # "SGMS1-AS1"
  ,"ENSG00000231588.1"  # "SHQ1P1"
  ,"ENSG00000231345.3"  # "BEND3P1"
  ,"ENSG00000232706.3"  # "NUTM2HP"
  ,"ENSG00000226168.1"  # "RP11-564C4.4"
  ,"ENSG00000230011.2"  # "CTSLP4"
  ,"ENSG00000213667.3"  # "PGGT1BP1"
  ,"ENSG00000204147.9"  # "ASAH2B"
  ,"ENSG00000148584.14" # "A1CF"
  ,"ENSG00000229711.1"  # "RP11-449O16.2"
  ,"ENSG00000236944.1"  # "CCDC58P2"
  ,"ENSG00000261368.1"  # "RP11-96B5.4"
  ,"ENSG00000223502.1"  # "RP11-96B5.3"
  ,"ENSG00000231132.1"  # "RP11-40C11.2"
  ,"ENSG00000207813.1"  # "MIR605"
  ,"ENSG00000235279.1"  # "RP11-539E19.2"
  ,"ENSG00000213659.4"  # "RSU1P3"
)

ENSGS = c(
  'ENSG00000185883.11' # ATP6V0C
  ,'ENSG00000140992.18' # PDPK1
)

ENSGS = c(
   'ENSG00000272888.5' # CHASERR
  ,'ENSG00000173575.20' # CHD2
)
  
sample_group = tibble(
  sample = c(
    'H110283',
    'H120109',
    'H131696',
    'H133133',
    'H133326',
    'H141189',
    'H144750',
    'H170738_1',
    'H170738_2',
    'H173468_1',
    'H173468_2',
    'H181707_1',
    'H181707_2',
    'M130158',
    'M170869',
    'M171873',
    'U151997',
    'U170082',
    'U170530',
    'U173381',
    'U181173',
    'U200375_1',
    'U200375_2'
  ),
  group = c(
    'ct_y',
    'ct_y',
    'ct_a',
    'ct_a',
    'ct_a',
    'ct_y',
    'ct_y',
    'pdpk',
    'pdpk',
    'tbl',
    'tbl',
    'chas',
    'chas',
    'dis_m',
    'dis_m',
    'dis_m',
    'dis_u',
    'dis_u',
    'dis_u',
    'id',
    'id',
    'id',
    'id'
  )
)


sample_group = tibble(
  sample = c(
    'H110283',
    'H120109',
    'H131696',
    'H133133',
    'H133326',
    'H141189',
    'H144750',
    'H170738',
    'H170738_1',
    'H170738_2',
    'H173468',
    'H173468_1',
    'H173468_2',
    'H181707',
    'H181707_1',
    'H181707_2',
    'M130158',
    'M170869',
    'M171873',
    'U151997',
    'U170082',
    'U170530',
    #'U173381',
    'U181173',
    'U200375',
    'U200375_1',
    'U200375_2'
  ),
  group = c(
    'ct_y',
    'ct_y',
    'ct_z', # ''ct_a',
    'ct_z', # ''ct_a',
    'ct_z', # ''ct_a',
    'ct_y',
    'ct_y',
    'dis_a', #'pdpk',
    'dis_a', #'pdpk',
    'dis_a', #'pdpk',
    'zzz', # 'dis_a', #'tbl',
    'zzz', # 'dis_a', #'tbl',
    'zzz', # 'dis_a', #'tbl',
    'dis_a', # 'zzz', # 'chas',
    'dis_a', # 'zzz', # 'chas',
    'dis_a', # 'zzz', # 'chas',
    'dis_a', #'dis_m',
    'dis_a', #'dis_m',
    'dis_a', #'dis_m',
    'dis_b', #'dis_u',
    'dis_b', #'dis_u',
    'dis_b', #'dis_u',
    #'dis_N', #'id',
    'dis_a', #'id',
    'dis_a', #'id',
    'dis_a', #'id',
    'dis_a' #'id'
  )
)



sample_group = tibble(
  sample = c(
    # '100222'
    '140093'
    ,'140430'
    #,'140698'
    ,'140761'
    ,'141097'
    ,'142033'
    ,'142083' 
    ,'142784' 
    ,'142947'
    ,'143620'
    ,'151767'
    ,'152149'
    ,'152281'
    ,'152761'
    ,'161937'
    ,'162981'
    ,'163464'
    ,'171863'
    ,'172478'
    ,'172844'
    ,'153552'
    ,'173502'
    #,'180767'
    ,'181707'
    #,'182848'
    #,'200132'
    #,'201535'
    #,'210488'
    #,'221122'
  ),
  group = c(
    # 'dis_n' #'miyatake'
    'ct_y' # ,'control'
    ,'ct_y' # ,'control'
    #,'ct_y' # ,'control'
    ,'dis_n' #,'ComplexSV1'
    ,'dis_n' #,'HH'
    ,'dis_n' #,'HH'
    ,'ct_y' # ,'control'
    ,'dis_n' #,'MECP2'
    ,'dis_n' #,'HH'
    ,'dis_n' #,'HH'
    ,'dis_n' #,'HH'
    ,'dis_n' #,'HH'
    ,'dis_n' #,'HH'
    ,'dis_n' #,'MEF2C'
    ,'dis_n' #,'ATP7A'
    ,'dis_n' #,'TAD'
    ,'dis_n' #,'SMC1A'
    ,'dis_n' #,'PURA'
    ,'dis_n' #,'AEBP2'
    ,'dis_n' #,'TAD2'
    ,'dis_n' #,'DLG4'
    ,'dis_n' #,'GATA4_NLGN3'
    #,'dis_n' #,'HH'
    ,'CHASERR'
    #,'dis_n' #,'ComplexSV2'
    #,'dis_n' #,'Consanguinity'
    #,'dis_n' #,'LAMB1'
    #,'blood'
    #,'blood'
  )
)


dt = read_tsv(DATA)

# PCA
res_pca = dt %>%
  select( -gene_id ) %>%
  mutate(across(everything(), ~. + sort(unique(.))[2] )) %>%
  mutate_all(~log2(.)) %>%
  t %>%
  prcomp( scale = T )

pcs = res_pca$x %>% 
  as.data.frame %>% 
  rownames_to_column(var = "sample") %>% 
  as_tibble %>%
  inner_join( sample_group ) 


# residualize
dt %>%
  filter( gene_id == ENSG ) %>%
  select( -gene_id ) %>%
  gather( key = sample, value = exp ) %>%
  inner_join( pcs ) %>%
  lm( exp ~ PC1 + PC2 + PC3, data = . ) %>%
  .$residuals
  

# pca plot
pcs %>%  
  #ggplot( aes( x = PC1, y = PC2, col = group ) ) +
  ggplot( aes( x = PC5, y = PC6, label = sample ) ) +
  geom_point() + 
  geom_text_repel()


# expression plot
dt2 = dt %>%
  #filter(gene_id == ENSG) %>%
  filter(gene_id %in% ENSGS) %>%
  select(-c(BAD_SAMPLES)) %>%
  # replicateマージオプション
  mutate(
    H181707 = ( H181707_1 + H181707_2 ) / 2, # chas
    #H170738 = ( H170738_1 + H170738_2 ) / 2, # pdpk
    H173468 = ( H173468_1 + H173468_2 ) / 2, # tbl
    U200375 = ( U200375_1 + U200375_2 ) / 2  # uchiya
  ) %>%
  #
  select( -c( 
    H181707_1, H181707_2, # chas
    #H170738_1, H170738_2, # pdpk
    H173468_1, H173468_2, # tbl
    U200375_1, U200375_2  # uchiya
  ) ) %>%
  gather(key = 'sample', value = 'exp', -gene_id) %>%
  # どのサンプルを選ぶか
  #inner_join( sample_group ) %>%
  #mutate(group = ifelse(sample == '181707', 'pt', 'ct')) %>% # chas
  #mutate(group = ifelse(sample == 'H181707' | sample == 'H181707_1' | sample == 'H181707_2', 'pt', 'ct')) %>% # chas
  mutate(group = ifelse(sample == 'H170738' | sample == 'H170738_1' | sample == 'H170738_2', 'pt', 'ct')) %>% # pdpk
  #mutate(group = ifelse(sample == 'H173468' | sample == 'H173468_1' | sample == 'H173468_2', 'pt', 'ct')) %>% # tbl1xr1
  #mutate(disease = ifelse(sample %in% HEALTHYS, 'healthy', 'disease')) %>%
  mutate(replicate = ifelse(grepl('_', sample), 'replicate', 'sample'))

# option: divide with ct mean
gene_MeanExpInCt = dt2 %>%
  group_by(gene_id, group) %>% 
  summarise(meanInCt = mean(exp)) %>% 
  #filter(group=='ct_z') %>% 
  #filter(group=='ct_y') %>% 
  filter(group=='ct') %>% 
  select(gene_id, meanInCt) %>% 
  ungroup %>% 
  # と発現量低いのを削る
  filter( meanInCt > 1)

dt2 = inner_join( dt2, gene_MeanExpInCt ) %>%
  mutate( exp = exp / meanInCt ) %>%
  select( -meanInCt ) # としないとこのコマンドを何回もしておかしくなる
#####

dt2 %>%
  #ggplot(aes(x = group, y = exp)) +
  #stat_summary(fun.y = 'mean', geom = 'bar', width = 0.7) + #, fill = 'white', color = 'black') +
  #geom_jitter(size = 1, width = 0.15, height = 0) + #, aes(color = disease)) +
  #stat_summary(fun.y = 'mean', fun.ymin = function(x)mean(x) - sd(x), fun.ymax = function(x)mean(x) + sd(x), geom = 'errorbar', width = 0.4) +
  #scale_color_manual(values = c('red', 'blue')) +
  #
  ggplot(aes(x = factor(gene_id), y = exp, fill = factor(group))) +
  stat_summary(fun.y = 'mean', geom = 'bar', width = 0.7, position = position_dodge(width = 0.8)) + 
  geom_point( position = position_jitterdodge(dodge.width = .8)  ) +
  #stat_summary(fun.y = 'mean', fun.ymin = function(x)mean(x) - sd(x), fun.ymax = function(x)mean(x) + sd(x), geom = 'errorbar', width = 0.4, position = position_dodge(width = 0.8)) +
  #scale_color_manual(values = c('red', 'blue')) +
  scale_x_discrete( limits = ENSGS ) + 
  #scale_shape_manual(values=c(4,16)) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4))) +
  #scale_fill_manual(limits = c('ct_y', 'dis_n', 'zzzz'), values = c(rgb(0,0,1,1), rgb(0.6,0.3,0.6,1), rgb(1,0,0,1))) +
  #scale_fill_manual(limits = c('ct_y', 'ct_z', 'dis_a', 'dis_b', 'zzz'), values = c( rgb(0,0,1,alpha=1), rgb(0,0,1,alpha=0.3), rgb(0.6,0.3,0.6,alpha=1), rgb(0.6,0.3,0.6,alpha=0.3), rgb(1,0,0,alpha=1) )) +
  #ylim(c(0,1.5)) +
  g -> gg

gg

OUT = paste0( ENSG, '.png' )
ggsave('tmp.png', plot = gg, width = 10, height = 5, units = 'cm', dpi = 600)

