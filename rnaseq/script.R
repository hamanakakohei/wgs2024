
source('C:/Users/hamanakakohei/Dropbox/backup20230402/_antares01_analysis_hamanaka_function_ggplot_setting.R')

library(tidyverse)

#ENSG = 'ENSG00000173575.20' # CHD2
#ENSG = 'ENSG00000272888.5' # CHASERR
#ENSG = 'ENSG00000177565.16' # TBL1XR1
#ENSG = 'ENSG00000140992.18' # PDPK1
#ENSG = 'ENSG00000185883.11' # ATP6V0C
#ENSG = 'ENSG00000183971.6' # NPW
ENSG = 'ENSG00000065054.13' # SLC9A3R2
ENSG = 'ENSG00000065057.7' # NTHL1
ENSG = 'ENSG00000103197.16' # TSC2
ENSG = 'ENSG00000008710.19' # PKD1
ENSG = 'ENSG00000167964.12' # RAB26
ENSG = 'ENSG00000131653.12' # TRAF7
ENSG = 'ENSG00000167971.15' # CASKIN1
ENSG = 'ENSG00000167965.17' # MLST8
ENSG = 'ENSG00000182685.7' # BRICD5
ENSG = 'ENSG00000184207.8' # PGP
ENSG = 'ENSG00000167967.15' # E4F1
ENSG = 'ENSG00000167968.12' # DNASE1L2
ENSG = 'ENSG00000167969.12' # ECI1
ENSG = 'ENSG00000205937.11' # RNPS1
ENSG = 'ENSG00000167972.13' # ABCA3
ENSG = 'ENSG00000162063.12' # CCNF
ENSG = 'ENSG00000162062.14' # C16orf59
ENSG = 'ENSG00000162068.1' # NTN3
ENSG = 'ENSG00000162065.12' # TBC1D24
ENSG = 'ENSG00000185883.11' # ATP6V0C
ENSG = 'ENSG00000162066.14' # AMDHD2
ENSG = 'ENSG00000140992.18' # PDPK1
ENSG = 'ENSG00000167977.8' # KCTD5
ENSG = 'ENSG00000172382.9' # PRSS27
ENSG = 'ENSG00000167978.16' # SRRM2
ENSG = 'ENSG00000103363.14' # ELOB
ENSG = 'ENSG00000103355.13' # PRSS33
ENSG = 'ENSG00000215148.8' # PRSS41
ENSG = 'ENSG00000007038.10' # PRSS21
ENSG = 'ENSG00000162078.11' # ZG16B
ENSG = 'ENSG00000005001.9' # PRSS22
ENSG = 'ENSG00000162076.12' # FLYWCH2
ENSG = 'ENSG00000059122.16' # FLYWCH1
ENSG = 'ENSG00000131650.13' # KREMEN2
ENSG = 'ENSG00000127564.16' # PKMYT1
ENSG = 'ENSG00000162073.13' # PAQR4

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
    'dis_n', #'dis_m',
    'dis_n', #'dis_m',
    'dis_n', #'dis_m',
    'dis_b', #'dis_u',
    'dis_b', #'dis_u',
    'dis_b', #'dis_u',
    'dis_n', #'id',
    'dis_n', #'id',
    'dis_n', #'id',
    'dis_n' #'id'
  )
)



sample_group = tibble(
  sample = c(
    # '100222'
     '140093'
    ,'140430'
    ,'140698'
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
    ,'180767'
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
    ,'ct_y' # ,'control'
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
    ,'dis_n' #,'HH'
    ,'CHASERR'
    #,'dis_n' #,'ComplexSV2'
    #,'dis_n' #,'Consanguinity'
    #,'dis_n' #,'LAMB1'
    #,'blood'
    #,'blood'
  )
)


DATA = 'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/sample27.tmm-norm_count.txt'
DATA = 'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/sample23.tmm-norm_count.expression.txt.gz'
dt = read_tsv(DATA)


# PCA
res_pca = dt %>%
  select( -gene_id ) %>%
  mutate(across(everything(), ~. + sort(unique(.))[2] )) %>%
  mutate_all(~log2(.)) %>%
  t %>%
  prcomp( scale = T )

## それか
#res_pca = readRDS( "C:/Users/hamanakakohei/Dropbox/wgs/難病ゲノム送付/Rproj_RNAseq/RNAseq_20230216.gtex.blindTrue.20240306.rds" )
#res_pca = readRDS( "C:/Users/hamanakakohei/Dropbox/wgs/難病ゲノム送付/Rproj_RNAseq/RNAseq_20240105.gtex.blindTrue.20240306.rds" )

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
  lm( exp ~ PC1, data = . ) %>%
  .$residuals
  

# pca plot
library("ggrepel")
pcs %>%  
  #ggplot( aes( x = PC1, y = PC2, col = group ) ) +
  ggplot( aes( x = PC5, y = PC6, label = sample ) ) +
  geom_point( shape = 16 ) + 
  geom_text_repel()
  scale_color_manual(values = c(
    rgb(1,0,0,alpha=1),
    #rgb(0,0,1,alpha=0.3),
    rgb(0,0,1,alpha=1),
    #rgb(0.6,0.3,0.6,alpha=0.3),
    rgb(0.6,0.3,0.6,alpha=1)
    #rgb(1,0,0,alpha=1),
    #rgb(1,0,0,alpha=1)
  )) + 
  theme( legend.position = "none" ) +
  g
  ggsave("tmp.png", width=4, height=4, units="cm", dpi=600)


z# expression plot
dt %>%
  filter(gene_id == ENSG) %>%
  select(-c(BAD_SAMPLES)) %>%
  mutate(
    H170738_1 = ( H170738_1 + H170738_2 ),
    H173468_1 = ( H173468_1 + H173468_2 ),
    H181707_1 = ( H181707_1 + H181707_2 ),
    U200375_1 = ( U200375_1 + U200375_2 )
  ) %>%
  select( -c( H170738_2, H173468_2, H181707_2, U200375_2 ) ) %>%
  gather(key = 'sample', value = 'exp', -gene_id) %>%
  #mutate(group = ifelse(sample == 'H181707_1' | sample == 'H181707_2', 'pt', 'ct')) %>% # chas
  #mutate(group = ifelse(sample == 'H170738_1' | sample == 'H170738_2', 'pt', 'ct')) %>% # pdpk
  mutate(group = ifelse(sample == 'H173468_1' | sample == 'H173468_2', 'pt', 'ct')) %>% # tbl1xr1
  mutate(disease = ifelse(sample %in% HEALTHYS, 'healthy', 'disease')) %>%
  #ggplot(aes(x = group, y = exp, fill = group)) +
  ggplot(aes(x = group, y = exp)) +
  #stat_summary(fun.y = 'mean', geom = 'bar', width = 0.7) +
  stat_summary(fun.y = 'mean', geom = 'bar', width = 0.7) + #, fill = 'white', color = 'black') +
  geom_jitter(size = 1, width = 0.15, height = 0) + #, aes(color = disease)) +
  stat_summary(fun.y = 'mean', fun.ymin = function(x)mean(x) - sd(x), fun.ymax = function(x)mean(x) + sd(x), geom = 'errorbar', width = 0.4) +
  scale_color_manual(values = c('red', 'blue')) +
  #theme(legend.position = 'none') +
  g -> gg

OUT = paste0( ENSG, '.png' )
ggsave(OUT, plot = gg, width = 6, height = 4.5, units = 'cm', dpi = 600)
