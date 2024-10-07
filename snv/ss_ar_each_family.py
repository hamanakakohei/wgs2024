import pandas as pd
import copy


# 0. functions
def add_utr_anno(exonic_func, varid, utr_varid_list):
  if varid in utr_varid_list:
    exonic_func = 'uorf_or_polya; ' + exonic_func 
  #
  return exonic_func

def judge_inhe(gt_pt,gt_fa,gt_mo):
  if   gt_pt=='1/1' and gt_fa=='0/1' and gt_mo=='0/1':
    return 'homo'
  elif gt_pt=='1/1' and gt_fa=='0/1' and gt_mo=='0/0':
    return 'fa_upd_or_mo_cnv'
  elif gt_pt=='1/1' and gt_fa=='0/0' and gt_mo=='0/1':
    return 'mo_upd_or_fa_cnv'
  elif gt_pt=='1/1':
    return 'homo'
  elif gt_pt=='0/1' and gt_fa=='0/1' and gt_mo=='0/0':
    return 'fa'
  elif gt_pt=='0/1' and gt_fa=='0/0' and gt_mo=='0/1':
    return 'mo'
  elif gt_pt=='0/1':
    return 'het'
  else:
    return 'others'

def judge_ar(df):
  inhes = list(df['inhe'])
  if 'homo' in inhes or ('fa' in inhes and 'mo' in inhes):
    return 'ok'
  elif len([i for i in inhes if i=='het'])>1:
    return 'ok'
  elif 'fa_upd_or_mo_cnv' in inhes or 'mo_upd_or_fa_cnv' in inhes:
    return 'ok'
  elif ('fa' in inhes or 'mo' in inhes) and 'het' in inhes:
    return 'ok'
  else:
    return 'no'


# 1. data prep
AR_ANNO = 'ar_anno.txt'
TOOLS = ['manta', 'canvas', 'lumpy', 'melt', 'cnvnator', 'delly']
TRIO = 'sample1211.plus40trio.txt'
DDD = 'DDG2P_20_3_2023.edit.tsv'
OMIM = 'genemap2.edit.txt'
UTR3 = 'varid__3utr.ncgm.txt'
UTR5 = 'varid__GeneUtrannotator.txt'


ar = pd.read_table(AR_ANNO, low_memory = False).rename({'sample': 'proband'}, axis = 1)
ar['gts'] = ar['gts'].map(lambda x: x.replace('1/0', '0/1', 10))
ddd = pd.read_table(DDD)
omim = pd.read_table(OMIM, names = ['Gene.refGene', 'omim'])
utr3 = [i.replace('-', '_') for i in list(pd.read_table(UTR3, names = ['varid'])['varid'])]
utr5 = [i.replace('-', '_') for i in list(pd.read_table(UTR5, names = ['varid'])['varid'])]

trio = pd.read_table(TRIO)
trio_probands = list(trio['pt'])

interest_cols1 = [
'Gene.refGene', 
'disease_name', 
'DDD_category', 
'allelic_requirement', 
'mutation_consequence', 
'omim'] 

interest_cols2 = [
'ar',
'inhe',
'ExonicFunc.refGene',
'spliceai',
'Func.refGene',
'GeneDetail.refGene',
'sv_caller', 
'AnnotSV_ID', 
'Location', 
'Dist_nearest_SS',
'controlSUPP', 
'AllSUPP', 
'gnomadSV_AF'] 


# 2. split into each sample
for SAMPLE in trio_probands:
  print(SAMPLE)
  out = 'ar_anno.plus_sv.' + SAMPLE + '.txt'
  # prep snv
  ar_one_fam = ar.query('proband == @SAMPLE')
  membs = ar_one_fam.membs.iloc[0].split(';')
  for i in range(len(membs)):
    memb = membs[i]
    ar_one_fam[memb] = ar_one_fam['gts'].map(lambda x: x.split(';')[i])
  # prep sv
  for TOOL in TOOLS:
    ANNOTSV = TOOL + '/annotsv/' + SAMPLE + '/annotsv.' + SAMPLE + '.final.MafFiltered.tsv'
    dt = pd.read_table(ANNOTSV, low_memory=False)
    samples = [i for i in dt.columns if i.startswith('DA000')]
    dt['proband'] = SAMPLE
    dt['sv_caller'] = TOOL
    dt['SV_chrom'] = 'chr' + dt['SV_chrom'].astype(str)
    dt['omim'] = dt['OMIM_phenotype'] + ';' + dt['OMIM_inheritance']
    dt = dt.\
      query('Annotation_mode == "split"')\
      [samples + ['omim', 'DDD_status', 'DDD_mode', 'DDD_consequence', 'DDD_disease', 'sv_caller', 'SV_chrom', 'SV_start', 'SV_end', 'proband', 'AnnotSV_ID', 'Gene_name', 'Location', 'Dist_nearest_SS', 'controlSUPP', 'AllSUPP', 'gnomadSV_AF']].\
      rename({'DDD_status': 'DDD_category', 'DDD_mode': 'allelic_requirement', 'DDD_consequence': 'mutation_consequence', 'DDD_disease': 'disease_name', 'SV_chrom': 'Chr', 'SV_start': 'Start', 'SV_end': 'End', 'Gene_name': 'Gene.refGene'}, axis = 1)
    #
    exec('dt_' + TOOL + ' = dt.copy()')
  # 
  dt = pd.concat([ar_one_fam, dt_manta, dt_canvas, dt_lumpy, dt_melt, dt_cnvnator, dt_delly]).\
    sort_values(by = ['Gene.refGene', 'Chr', 'Start'])
  #
  dt['gt_pt'] = dt[SAMPLE].map(lambda x: x.split(':')[0].replace('./.', '0/0'))
  fa = trio.query('pt == @SAMPLE')['fa'].values[0]
  if fa == '0':
    dt['gt_fa'] = '0/0'
  else:
    dt['gt_fa'] = dt[fa].map(lambda x: x.split(':')[0].replace('./.', '0/0'))
  #
  mo = trio.query('pt == @SAMPLE')['mo'].values[0]
  if mo == '0':
    dt['gt_mo'] = '0/0'
  else:
    dt['gt_mo'] = dt[mo].map(lambda x: x.split(':')[0].replace('./.', '0/0'))
  #
  dt['inhe'] = dt.apply(lambda df: judge_inhe(df['gt_pt'],df['gt_fa'],df['gt_mo']       ),axis=1)
  ##
  #gene_gt = dt[['Gene.refGene', SAMPLE]] 
  #gene_gt['allele_count'] = gene_gt[SAMPLE].map(lambda gt: sum([int(i) if i != '.' else 0 for i in gt.split(':')[0].split('/')]))  
  #two_allele_genes = list(gene_gt.groupby('Gene.refGene').sum().query('allele_count >= 2').index)
  #
  gt_cols = [i for i in dt.columns if i.startswith('DA000')]
  #
  other_cols = [i for i in dt.columns if i not in gt_cols and i not in interest_cols1 and i not in interest_cols2]
  #
  dt_snv = dt.query('Location != Location')
  dt_sv = dt.query('Location == Location')
  dt_sv['Location1'] = dt_sv['Location'].map(lambda x: x.split('-')[0])
  dt_sv['Location2'] = dt_sv['Location'].map(lambda x: x.split('-')[1])
  #
  gene_ar = dt.groupby(['Gene.refGene']).apply(lambda df: judge_ar(df)).reset_index().rename({0: 'ar'} ,axis=1)
  dt = pd.merge(dt, gene_ar)
  dt = pd.merge(dt.drop(['disease_name', 'DDD_category', 'allelic_requirement', 'mutation_consequence'], axis = 1), ddd, how = 'left')
  print(len(dt))
  print(len(dt.columns))
  dt = pd.merge(dt.drop(['omim'], axis = 1), omim, how = 'left')
  print(len(dt))
  print(len(dt.columns))
  dt['ExonicFunc.refGene'] = dt.apply(lambda df: add_utr_anno(df['ExonicFunc.refGene'], df['varid'], (utr5 + utr3)), axis=1)
  dt[interest_cols1 + gt_cols + interest_cols2 + other_cols].\
    query('Chr != "chrX" and Chr != "chrY"').\
    query('gt_pt != "0/0"').\
    to_csv(out, sep = '\t', index = False)
    #query('`Gene.refGene` in @two_allele_genes').\


