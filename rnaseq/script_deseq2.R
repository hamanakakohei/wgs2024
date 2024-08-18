
source('C:/Users/hamanakakohei/Dropbox/backup20230402/_antares01_analysis_hamanaka_function_DESeq2.R')

## 2024

FILES = c(
   'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H110283/H110283.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H120109/H120109.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H131696/H131696.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H133133/H133133.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H133326/H133326.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H141189/H141189.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H144750/H144750.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H170738_1/H170738_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H170738_2/H170738_2.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H173468_1/H173468_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H173468_2/H173468_2.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H181707_1/H181707_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/H181707_2/H181707_2.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/M130158/M130158.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/M170869/M170869.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/M171873/M171873.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U151997/U151997.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U170082/U170082.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U170530/U170530.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U173381/U173381.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U181173/U181173.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U200375_1/U200375_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/U200375_2/U200375_2.ReadsPerGene.out.tab'
)

FILES = c(
   'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H110283/H110283.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H120109/H120109.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H131696/H131696.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H133133/H133133.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H133326/H133326.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H141189/H141189.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H144750/H144750.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H170738_1/H170738_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H170738_2/H170738_2.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H173468_1/H173468_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H173468_2/H173468_2.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H181707_1/H181707_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/H181707_2/H181707_2.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/M130158/M130158.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/M170869/M170869.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/M171873/M171873.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U151997/U151997.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U170082/U170082.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U170530/U170530.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U173381/U173381.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U181173/U181173.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U200375_1/U200375_1.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20240105/gtex_method/U200375_2/U200375_2.ReadsPerGene.out.tab'
)

FILES = c(
   '/home/kohei.hamanaka/wgs/RNAseq_20240105/H110283/H110283.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H120109/H120109.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H131696/H131696.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H133133/H133133.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H133326/H133326.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H141189/H141189.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H144750/H144750.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H170738_1/H170738_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H170738_2/H170738_2.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H173468_1/H173468_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H173468_2/H173468_2.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H181707_1/H181707_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/H181707_2/H181707_2.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/M130158/M130158.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/M170869/M170869.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/M171873/M171873.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U151997/U151997.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U170082/U170082.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U170530/U170530.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U173381/U173381.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U181173/U181173.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U200375_1/U200375_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/U200375_2/U200375_2.ReadsPerGene.out.tab'
)

FILES = c(
   '/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H110283/H110283.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H120109/H120109.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H131696/H131696.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H133133/H133133.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H133326/H133326.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H141189/H141189.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H144750/H144750.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H170738_1/H170738_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H170738_2/H170738_2.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H173468_1/H173468_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H173468_2/H173468_2.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H181707_1/H181707_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/H181707_2/H181707_2.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/M130158/M130158.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/M170869/M170869.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/M171873/M171873.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U151997/U151997.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U170082/U170082.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U170530/U170530.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U173381/U173381.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U181173/U181173.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U200375_1/U200375_1.ReadsPerGene.out.tab'
  ,'/home/kohei.hamanaka/wgs/RNAseq_20240105/gtex_method/U200375_2/U200375_2.ReadsPerGene.out.tab'
)

WELLS = c(
  'H110283'
  ,'H120109'
  ,'H131696'
  ,'H133133'
  ,'H133326'
  ,'H141189'
  ,'H144750'
  ,'H170738_1'
  ,'H170738_2'
  ,'H173468_1'
  ,'H173468_2'
  ,'H181707_1'
  ,'H181707_2'
  ,'M130158'
  ,'M170869'
  ,'M171873'
  ,'U151997'
  ,'U170082'
  ,'U170530'
  ,'U173381'
  ,'U181173'
  ,'U200375_1'
  ,'U200375_2'
)

SAMPLES = c(
  'H110283'
  ,'H120109'
  ,'H131696'
  ,'H133133'
  ,'H133326'
  ,'H141189'
  ,'H144750'
  ,'H170738_1'
  ,'H170738_2'
  ,'H173468_1'
  ,'H173468_2'
  ,'H181707_1'
  ,'H181707_2'
  ,'M130158'
  ,'M170869'
  ,'M171873'
  ,'U151997'
  ,'U170082'
  ,'U170530'
  ,'U173381'
  ,'U181173'
  ,'U200375_1'
  ,'U200375_2'
)


CONDITIONS = c(
   'ct_y'
  ,'ct_y'
  ,'ct_a'
  ,'ct_a'
  ,'ct_a'
  ,'ct_y'
  ,'ct_y'
  ,'pdpk'
  ,'pdpk'
  ,'tbl'
  ,'tbl'
  ,'chas'
  ,'chas'
  ,'dis_m'
  ,'dis_m'
  ,'dis_m'
  ,'dis_u'
  ,'dis_u'
  ,'dis_u'
  ,'id'
  ,'id'
  ,'id'
  ,'id'
)

##########################

## 2023

# gtex

FILES = c(
# '/home/kohei.hamanaka/wgs/RNAseq_20230216/100222/100222.ReadsPerGene.out.tab'
 '/home/kohei.hamanaka/wgs/RNAseq_20230216/140093/140093.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/140430/140430.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/140698/140698.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/140761/140761.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/141097/141097.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/142033/142033.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/142083/142083.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/142784/142784.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/142947/142947.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/143620/143620.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/151767/151767.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/152149/152149.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/152281/152281.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/152761/152761.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/161937/161937.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/162981/162981.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/163464/163464.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/171863/171863.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/172478/172478.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/172844/172844.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/153552/153552.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/173502/173502.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/180767/180767.ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216/181707/181707.ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216/182848/182848.ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216/200132/200132.ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216/201535/201535.ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216/210488/210488.ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216/221122/221122.ReadsPerGene.out.tab'
)

FILES = c(
  # 'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/100222/100222.ReadsPerGene.out.tab'
   'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/140093/140093.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/140430/140430.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/140698/140698.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/140761/140761.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/141097/141097.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/142033/142033.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/142083/142083.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/142784/142784.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/142947/142947.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/143620/143620.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/151767/151767.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/152149/152149.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/152281/152281.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/152761/152761.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/161937/161937.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/162981/162981.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/163464/163464.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/171863/171863.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/172478/172478.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/172844/172844.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/153552/153552.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/173502/173502.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/180767/180767.ReadsPerGene.out.tab'
  ,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/181707/181707.ReadsPerGene.out.tab'
  #,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/182848/182848.ReadsPerGene.out.tab'
  #,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/200132/200132.ReadsPerGene.out.tab'
  #,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/201535/201535.ReadsPerGene.out.tab'
  #,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/210488/210488.ReadsPerGene.out.tab'
  #,'C:/Users/hamanakakohei/Dropbox/wgs/tmpcopy0114-2/RNAseq_20230216/gtex_method/221122/221122.ReadsPerGene.out.tab'
)

# default
 
FILES = c(
# '/home/kohei.hamanaka/wgs/RNAseq_20230216_default/100222/100222ReadsPerGene.out.tab'
 '/home/kohei.hamanaka/wgs/RNAseq_20230216_default/140093/140093ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/140430/140430ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/140698/140698ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/140761/140761ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/141097/141097ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/142033/142033ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/142083/142083ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/142784/142784ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/142947/142947ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/143620/143620ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/151767/151767ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/152149/152149ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/152281/152281ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/152761/152761ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/161937/161937ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/162981/162981ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/163464/163464ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/171863/171863ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/172478/172478ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/172844/172844ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/153552/153552ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/173502/173502ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/180767/180767ReadsPerGene.out.tab'
,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/181707/181707ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/182848/182848ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/200132/200132ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/201535/201535ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/210488/210488ReadsPerGene.out.tab'
#,'/home/kohei.hamanaka/wgs/RNAseq_20230216_default/221122/221122ReadsPerGene.out.tab'

)

WELLS = c(
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
)

SAMPLES = c(
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
)

CONDITIONS = c(
# 'miyatake'
 'control'
,'control'
,'control'
,'ComplexSV1'
,'HH'
,'HH'
,'control'
,'MECP2'
,'HH'
,'HH'
,'HH'
,'HH'
,'HH'
,'MEF2C'
,'ATP7A'
,'TAD'
,'SMC1A'
,'PURA'
,'AEBP2'
,'TAD2'
,'DLG4'
,'GATA4_NLGN3'
,'HH'
,'CHASERR'
#,'ComplexSV2'
#,'Consanguinity'
#,'LAMB1'
#,'blood'
#,'blood'
)

###########################

# output PCA plot
META <- function(FILES=...,WELLS=...,SAMPLES=...,CONDITIONS=...,TYPE="DESeq2"){
    # make coutData,colData objects
    rm(countData)
    for(FILE in FILES){
        if("countData" %in% ls()){
            if(TYPE=="DESeq2"){
                DATA = read.table(FILE,sep="\t",skip=4,header=F)
                countData = cbind(countData,DATA[,4])
            }else{
                DATA = read.table(FILE)
                countData = cbind(countData,round(DATA[,1]))
            }
        } else {
            if(TYPE=="DESeq2"){
                DATA = read.table(FILE,sep="\t",skip=4,header=F)
                countData = DATA[,4]
                names(countData) = DATA[,1]
            }else{
                DATA = read.table(FILE)
                countData = round(DATA[,1])
                names(countData) = rownames(DATA)
            }
        }
    }
    colnames(countData) = WELLS
    
    colData = data.frame(
        row.names = WELLS,
        condition = CONDITIONS,
        sample = SAMPLES,
        run = WELLS
    )
    
    return(list(countData,colData))
}


QC <- function(FILES=...,WELLS=...,SAMPLES=...,CONDITIONS=...,OUTPUT="output.png",OPTION="PCA",TYPE="DESeq2"){
  LIST = META(FILES=FILES,WELLS=WELLS,SAMPLES=SAMPLES,CONDITIONS=CONDITIONS,TYPE=TYPE)
  
  # DEG analysis
  library("DESeq2")
  DDS = DESeqDataSetFromMatrix(
    countData = LIST[[1]],
    colData = LIST[[2]],
    design = ~ condition
  )
  DDS = DESeq(DDS)
  RES = results(DDS)
  VSD = varianceStabilizingTransformation(DDS,blind=TRUE)
  
  if (OPTION == "PCA"){
    # PCA
    library(stats)
    MAT = data.matrix(assay(VSD))
    high_genes     = names( sort( apply(MAT,1,sum), decreasing = T ) )[1:10000]
    MAT = MAT[high_genes, ]
    variable_genes = names( sort( apply(MAT,1,var), decreasing = T ) )[1:1000]
    MAT = MAT[variable_genes, ]
    PCA = prcomp(t(MAT), scale = T)
    return( PCA )
  }
}

res_pca = QC( 
  FILES=FILES,
  WELLS=WELLS,
  SAMPLES=SAMPLES,
  CONDITIONS=CONDITIONS,
  OUTPUT="tmp.png",
  OPTION="PCA"
)

saveRDS( res_pca, "RNAseq_20230216.gtex.blindTrue.20240306.rds" )
saveRDS( res_pca, "RNAseq_20230216.gtex.blindFalse.20240306.rds" )

#saveRDS( res_pca, "RNAseq_20230216_default.blindTrue.rds" )
#saveRDS( res_pca, "RNAseq_20230216_default.blindFalse.rds" )

#saveRDS( res_pca, "RNAseq_20240105.default.blindTrue.20240306.rds" )
#saveRDS( res_pca, "RNAseq_20240105.default.blindFalse.20240306.rds" )

saveRDS( res_pca, "RNAseq_20240105.gtex.blindTrue.20240306.rds" )
saveRDS( res_pca, "RNAseq_20240105.gtex.blindFalse.20240306.rds" )
