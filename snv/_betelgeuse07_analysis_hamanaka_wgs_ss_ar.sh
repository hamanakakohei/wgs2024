import pandas as pd
import numpy as np

#PED='sample1192.ped'
PED='sample1211.plus40trio.ped'
DDD='/betelgeuse04/analysis/hamanaka/resource/DDG2P_14_10_2021.tsv'
ARDATA='ar_anno.txt'

# make gene_inhe (1 gene, mt 1 rows,,,,,
ddd = pd.read_table(DDD)
tmp1 = ddd[['gene_symbol','allelic_requirement']]
tmp2 = ddd[['prev_symbols','allelic_requirement']].query('prev_symbols==prev_symbols')
tmp3 = pd.DataFrame(columns=tmp2.columns)
for idx,row in tmp2.iterrows():
  for one_gene in row['prev_symbols'].split(';'):
    tmp3 = tmp3.append({'prev_symbols':one_gene,'allelic_requirement':row['allelic_requirement']},ignore_index=True)

gene_inhe = pd.concat([tmp1.rename({'gene_symbol':'Gene.refGene'},axis=1),tmp3.rename({'prev_symbols':'Gene.refGene'},axis=1)])
#################

ped = pd.read_table(PED, names=('fam','pt','fa','mo','sex','status')) #.query('fam==pt')
dt0 = pd.read_table(ARDATA, low_memory=False)
varid_annos = dt0.drop(['sample', 'membs', 'gts'], axis=1).drop_duplicates(keep = 'first')
dt = dt0.copy()

# extract varid OI
varid_oi_5utr = pd.read_table('varid__GeneUtrannotator.txt', names=['a'])['a'].tolist()
varid_oi_3utr = pd.read_table('varid__3utr.ncgm.txt',        names=['a'])['a'].tolist()
varid_spliceai = dt0[['varid','Func.refGene','spliceai']].copy() #.query('`Func.refGene`.str.contains("UTR") or `Func.refGene`.str.contains("intronic")').copy()
varid_spliceai['spliceai'] = varid_spliceai['spliceai'].map(lambda x: extract_spliceai(x))
varid_oi_spliceai = set(varid_spliceai.query('spliceai>0.15')['varid'].to_list())
#varid_oi1 = set(dt['varid'].to_list())
##################

dt = dt.fillna({'allelic_requirement':'aaa'}).query('allelic_requirement.str.contains("biallelic")')[['Gene.refGene','varid','sample','membs','gts']] 
dt['gt_pt'] = dt.apply(lambda df: split_gts(df['sample'],df['membs'],df['gts'],ped,'pt'),axis=1)
dt['gt_fa'] = dt.apply(lambda df: split_gts(df['sample'],df['membs'],df['gts'],ped,'fa'),axis=1)
dt['gt_mo'] = dt.apply(lambda df: split_gts(df['sample'],df['membs'],df['gts'],ped,'mo'),axis=1)
#dt = pd.concat([dt,annotsv_all])
dt['gt_pt'] = dt['gt_pt'].map(lambda x: correct_gt(x))
dt['gt_fa'] = dt['gt_fa'].map(lambda x: correct_gt(x))
dt['gt_mo'] = dt['gt_mo'].map(lambda x: correct_gt(x))
dt['inhe']  = dt.apply(lambda df: judge_inhe(df['gt_pt'],df['gt_fa'],df['gt_mo']       ),axis=1)
#dt.apply(lambda df: split_genes(df),axis=1)

dt_SinglGene = dt.query('not `Gene.refGene`.str.contains(";")')
dt_MultiGene = dt.query('`Gene.refGene`.str.contains(";")')
dt_MultiGene2 = pd.DataFrame(columns=dt_MultiGene.columns)
for idx,row in dt_MultiGene.iterrows():
  aa = split_one_col(row.copy(),'Gene.refGene')
  dt_MultiGene2 = pd.concat([dt_MultiGene2,aa])

varid_oi = set(varid_oi_spliceai) | set(varid_oi_sv) | set(varid_oi_3utr) | set(varid_oi_5utr)
dt_all = pd.concat([dt_SinglGene, dt_MultiGene2, annotsv_all[['Gene.refGene', 'varid', 'sample', 'gt_pt', 'gt_fa', 'gt_mo']]]      )
sample_gene_ArJudg = dt_all.groupby(['Gene.refGene','sample']).apply(lambda df: judge_ar(   df         )).reset_index().rename({0:'ar'   },axis=1)
sample_gene_varoi  = dt_all.groupby(['Gene.refGene','sample']).apply(lambda df: judge_varoi(df,varid_oi)).reset_index().rename({0:'varoi'},axis=1)
sample_gene_varid__oi = pd.merge(pd.merge(pd.merge(dt_all, sample_gene_ArJudg), sample_gene_varoi).query('ar=="ok" & varoi=="ok"'),gene_inhe.query('allelic_requirement=="biallelic"').drop_duplicates())[['sample','Gene.refGene','varid']].drop_duplicates()

res = pd.concat([
  pd.merge(sample_gene_varid__oi, annotsv_all),
  pd.merge(sample_gene_varid__oi, dt0.rename({'Gene.refGene':'Gene.refGene_multi'}, axis=1), on=['sample', 'varid'])
])
res = pd.merge(varid_vqslod, res, how='right')
res = pd.merge(varid_polya,  res, how='right')
res = pd.merge(varid_5utr,   res, how='right')
res.to_csv('ar_anno__oi.txt',index=False,sep='\t')


# snv additional annotation
VQSLOD = '/betelgeuse07/analysis/hamanaka/wgs/varid__VqslodPass.txt'
UTR5 = '/betelgeuse07/analysis/hamanaka/wgs/varid_gene_utrannotator.txt'
varid_vqslod = pd.read_table(VQSLOD, names=['varid'])
varid_vqslod['vqslod'] = 'PASS'
varid_polya = pd.DataFrame({'varid': varid_oi_3utr, 'polya': 'yes'})
varid_5utr = pd.read_table(UTR5, names=['varid', 'Gene.refene', '5utr'])[['varid', '5utr']]

# prep annotsv
annotsv_all = pd.DataFrame(columns=['Gene.refGene','varid','sample','gt_pt','gt_fa','gt_mo','Location','Location2','controlSUPP','AllSUPP'])
for idx,row in ped.query('fam==pt').iterrows():
  pt = row['pt']
  fa = row['fa']
  mo = row['mo']
  print(pt)
  annotsv_path = '/betelgeuse07/analysis/hamanaka/wgs_sv/manta/annotsv/'+pt+'/annotsv.'+pt+'.final.MafFiltered.tsv'
  annotsv = pd.read_table(annotsv_path,index_col=0, low_memory=False).query('Annotation_mode=="split"').rename({pt:'gt_pt','AnnotSV_ID':'varid','Gene_name':'Gene.refGene'},axis=1)
  annotsv['sample'] = pt
  if fa=='0' or fa!=fa or fa not in annotsv.columns:
    annotsv['gt_fa'] = np.nan
  else:
    annotsv = annotsv.rename({fa:'gt_fa'},axis=1)
  if mo=='0' or mo!=mo or mo not in annotsv.columns:
    annotsv['gt_mo'] = np.nan
  else:
    annotsv = annotsv.rename({mo:'gt_mo'},axis=1)
  #
  annotsv_all = pd.concat([annotsv_all,annotsv[['Gene.refGene','varid','sample','gt_pt','gt_fa','gt_mo','Location','Location2','controlSUPP','AllSUPP']]])

annotsv_all['gt_pt'] = annotsv_all['gt_pt'].map(lambda x: correct_gt_annotsv(x))
annotsv_all['gt_fa'] = annotsv_all['gt_fa'].map(lambda x: correct_gt_annotsv(x))
annotsv_all['gt_mo'] = annotsv_all['gt_mo'].map(lambda x: correct_gt_annotsv(x))
annotsv_all['inhe' ] = annotsv_all.apply(lambda df: judge_inhe(df['gt_pt'],df['gt_fa'],df['gt_mo']       ),axis=1)
annotsv_all = pd.merge(annotsv_all, gene_varid)
###############

# select damaging sv
for idx,row in ped.iterrows():
  pt = row['pt']
  print(pt)
  annotsv_path = '/betelgeuse07/analysis/hamanaka/wgs_sv/manta/annotsv/'+pt+'/annotsv.'+pt+'.final.MafFiltered.tsv'
  try:
    annotsv = pd.read_table(annotsv_path, index_col=0, low_memory=False).query('Annotation_mode=="split"')
    gene_varid__OneSample = annotsv.apply(lambda se: extract_deleterious_sv(se), axis=1)
    #print(gene_varid__OneSample)
    if 'gene_varid' in locals():
      gene_varid = pd.concat([gene_varid, gene_varid__OneSample])
    else:
      gene_varid = gene_varid__OneSample
  except:
    pass

gene_varid = gene_varid.query('AnnotSV_ID==AnnotSV_ID').rename({'AnnotSV_ID':'varid', 'Gene_name':'Gene.refGene'}, axis=1)[['varid', 'Gene.refGene']].drop_duplicates()
varid_oi_sv = set(gene_varid['varid'])
#####################

def judge_varoi(df,varid_oi):
  if len(set(varid_oi) & set(df['varid'].to_list()))>0:
    return 'ok'
  else:
    return 'no'

def extract_deleterious_sv(se):
  SvType = se['SV_type']
  SvSize = abs(se['SV_length'])
  Loc1 = se['Location'].split('-')[0]
  Loc2 = se['Location'].split('-')[1]
  if (Loc1!=Loc2 or 'intron' not in Loc1) and (SvSize<1000000 or SvType=='INV') and SvType!='BND':
    return se[['AnnotSV_ID', 'Gene_name']]
    #return se['AnnotSV_ID']

def extract_spliceai(spais):
  if spais==spais:
    maxs = []
    for spai in spais.split(';')[0:-1]:
      maxOne = max([float(ds) for ds in spai.split('|')[2:6]])
      maxs = maxs+[maxOne]
    return max(maxs)
  else:
    return np.nan

def split_one_col(se,COLUMN):
  genes = se[COLUMN] # gene toha kagiranai
  if len(genes.split(';'))==1:
    return pd.DataFrame([se])
  else:
    tmp = pd.DataFrame(columns=se.index)
    for gene in genes.split(';'):
      se[COLUMN] = gene
      tmp = tmp.append(se,ignore_index=True)
    return tmp

def judge_ar(df):
  inhes = df['inhe']
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

def judge_inhe(gt_pt,gt_fa,gt_mo):
  if   gt_pt=='1/1' and gt_fa=='0/1' and gt_mo=='0/1':
    return 'homo'
  elif gt_pt=='1/1' and gt_fa=='0/1' and gt_mo=='0/0':
    return 'fa_upd_or_mo_cnv'
  elif gt_pt=='1/1' and gt_fa=='0/0' and gt_mo=='0/1':
    return 'mo_upd_or_fa_cnv'
  elif gt_pt=='0/1' and gt_fa=='0/1' and gt_mo=='0/0':
    return 'fa'
  elif gt_pt=='0/1' and gt_fa=='0/0' and gt_mo=='0/1':
    return 'mo'
  elif gt_pt=='1/1' and gt_fa!=gt_fa and gt_mo!=gt_mo:
    return 'homo'
  elif gt_pt=='1/1' and gt_fa=='0/1' and gt_mo!=gt_mo:
    return 'homo'
  elif gt_pt=='1/1' and gt_fa!=gt_fa and gt_mo=='0/1':
    return 'homo'
  elif gt_pt=='0/1' and gt_fa=='0/1' and gt_mo!=gt_mo:
    return 'fa'
  elif gt_pt=='0/1' and gt_fa!=gt_fa and gt_mo=='0/1':
    return 'mo'
  elif gt_pt=='0/1' and gt_fa!=gt_fa and gt_mo!=gt_mo:
    return 'het'
  elif gt_pt=='0/1' and gt_fa!=gt_fa and gt_mo=="0/0":
    return 'het'
  elif gt_pt=='0/1' and gt_fa=="0/0" and gt_mo!=gt_mo:
    return 'het'
  else:
    return 'others'

def correct_gt(gt):
  if gt=='1/0':
    return '0/1'
  elif gt=='./.':
    return np.nan
  else:
    return gt

def correct_gt_annotsv(gt):
  try:
    gt = gt.split(':')[0]
    if gt=='./.':
      return np.nan
    else:
      return gt
  except:
    return np.nan

def split_gts(SAMPLE,membs,gts,ped,who='fa'):
  DA_of_who = ped.query('pt==@SAMPLE')[who].iloc[0]
  try:
    index_of_who = membs.split(';').index(DA_of_who)
    return gts.split(';')[index_of_who].split(':')[0]
  except:
    return np.nan 


# recessive
  
  JointVcf_to_EachFam3 jointcall.VQSR.chr$CHR.n2.vcf.gz $PEDREFORMAT
  cut -f2 PEDREFORMAT|while read PT; do
    exclude_no_carrier_variant2 $PED.fam.vcf $PT|awk '$1=="#CHROM"{for(i=10;i<=NF;i++){MEMB=$i"\t"};sub("\t$","",MEMB)}$1!~/^#/{print $1"-"$2"-"$3"-"$4"\t"MEMB}' > id.var.$PT_members.$CHR.txt
    join <(sort -u id.maf0.01.$CHR.txt) <(cut -f1 id.var.$PT_members.$CHR.txt|sort -u) > id.maf0.01.var.$PT.$CHR.txt
    join -t$'\t' id.maf0.01.var.$PT.$CHR.txt id.var.$PT_members.$CHR.txt > id.maf0.01.var.$PT_members.$CHR.txt
    source /usr/local/genome/gcc-8.5.0/env.sh
    merge_denovoqc4.R id.maf0.01.var.$PT_members.$CHR.txt $PT.maf0.01.$CHR.txt 
    annotate_denovo.R $PT.maf0.01.$CHR.txt exome_test_refGene.hg38_multianno.edit.txt $PT.maf0.01.$CHR.anno.txt
    awk 'BEGIN{OFS="\t"}NR>1{split($1,A,"_"); print $2,A[1],A[2],A[3],A[4]; print $3,A[1],A[2],A[3],A[4]; print $4,A[1],A[2],A[3],A[4]}' $PT.maf0.01.$CHR.anno.txt > $PT_chr_pos_ref_alt.txt
    extract_genotype_from_variant_sample_combi2 $PT_chr_pos_ref_alt.txt jointcall.VQSR.chr$CHR.n2.vcf.gz > $PT_varid_gt.txt
    annotate_denovo3.R $PT.maf0.01.$CHR.anno.txt $PT_varid_gt.txt $DDD $OMIM $PLI $PT.maf0.01.$CHR.anno2.txt
  done

  /antares01/analysis/hamanaka/function/get_ArVarVcf_ForManualSpai.py --annovar --VcfHead --out




#${ANNOVAR}convert2annovar.pl    -format vcf4 -allsample -withfreq ../jointcall.VQSR.chr$CHR.n2.vcf.gz > $CHR.avinput
#${ANNOVAR}annotate_variation.pl -filter -dbtype generic -buildver hg38 -outfile $CHR.gnomad $CHR.avinput $DB       -genericdbfile hg38_gnomad30_genome.txt -otherinfo
#${ANNOVAR}annotate_variation.pl -filter -dbtype vcf     -buildver hg38 -outfile $CHR.tommo  $CHR.avinput $RESOURCE -vcfdbfile     $TOMMOVCF                -otherinfo
#${ANNOVAR}table_annovar.pl      -filter -dbtype vcf     -buildver hg38 -outfile $CHR.tommo  ../jointcall.VQSR.chr$CHR.n2.vcf.gz $RESOURCE -vcfdbfile     $TOMMOVCF                -otherinfo
#make -f Makefile.annovar_whole_genome refGene INPUT=jointcall.VQSR.chr$CHR.n2.vcf.gz > log.annovar.$CHR.txt 2>&1

select_sample_and_calc_af vqsr_filter/n.vcf.gz $SAMPLEHEALTHY               > id_maf.ycu.ctrl.txt
awk '$2<0.0001{print $1}' id_maf.ycu.ctrl.txt                               > id.rare.ycu.txt
extract_rarevariantid_from_exomsummary $EXOMESUMMARY 0.00005                > id.rare.exac.txt



VCF="../vqsr_filter/n.vcf"
FUNCTION="/betelgeuse01/analysis/hamanaka/function/"
PED="/archive3/hamanaka/resource/2671trio_20191021.sexedit.txt"
PED2="/betelgeuse04/analysis/hamanaka/jointgt.13851/denovo/2536trio.invcf.sex.ped"
SAMPLE_BAMPATH="/betelgeuse04/analysis/hamanaka/jointgt.13851/denovo/2536trio_bampath.txt"
IDRAREEXAC=/betelgeuse04/analysis/hamanaka/jointgt.13851/id.rare.exac.txt
IDRAREYCU=/betelgeuse04/analysis/hamanaka/jointgt.13851/id.rare.ycu.txt
VCFRAREDNV="/betelgeuse04/analysis/hamanaka/jointgt.13851/vqsr_filter/nrd.vcf"
SANGER="/archive3/hamanaka/resource/denovo.sanger.confirmed20191201.txt"
OUTVARIANTQC=denovoqc.20200206.txt
OUTPNG=denovoqc.20200206.png
OUTSAMPLEQC=sampleqc.denovo.20200206.txt
export PERL5LIB=$PERL5LIB:/usr/local/genome/vcftools-0.1.17/share/perl5/
export LD_LIBRARY_PATH=/usr/local/genome/zlib-1.2.11/lib:$LD_LIBRARY_PATH
export PATH=/betelgeuse01/analysis/hamanaka/function/:$PATH
source VcfFunctions.sh

# trio sample in vcf?
join <(awk '{print $2,$1}' $PED|sort) <(head -n1000 $VCF|awk '$1 ~ /^#CHROM/'|cut -f10-|tr '\t' '\n'|sort)|cut -d" " -f2|sort|uniq -c|awk '$1==3{print $2}'|sort|join -t$'\t' <(sort $PED) - > 2536trio.invcf.ped
join -1 2 -2 1 <(sort -k2,2 2536trio.invcf.ped) <(awk 'NR>1{print $1,$10}' ../sampleqc.20191121.txt|sed -e 's/female/2/' -e's/male/1/'|sort)|awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$4,$6}'|sort|awk '$3==0{print $0"\t"1}$3!=0{print $0"\t"2}' > 2536trio.invcf.sex.ped
#=> seibetu outlier wo manual de naoshita
#=> ryoushin no seibetu ga onaji family:
#Sample_7118; #Sample_5338,7899,7900; #Sample_13765; #Sample_16788
#hottannsha ni kanshite koushiteru: ex. XXY -> M, XO -> M; maaikka

#pheonptype col. wo tuikashita, dnmfilter igai mo ugoku ka kakuninn iru
join <(sort $IDRAREEXAC) <(sort $IDRAREYCU)|sort|join -t$'\t' - <(sort id.dnv_pt_fa_mo.txt) > id.rare.dnv_pt_fa_mo.txt
cut -f1 id.rare.dnv_pt_fa_mo.txt                                                            > id.rare.dnv.txt
subset_vcf $VCF id.rare.dnv.txt                                                             > ../vqsr_filter/nrd.vcf
extract_vqslod2 $VCFRAREDNV                                                                 > id.rare.dnv_vq.txt 
TrioDenovo.sh $VCFRAREDNV $PED2                                                             > id.rare.dnv_td.txt
DNMfilter.sh id.rare.dnv_td.txt $PED2 $SAMPLE_BAMPATH                                       > id.rare.dnv_dn.txt
denovogear.sh $VCFRAREDNV $PED2                                                               id.rare.dnv_dg.txt
denovofilterFormat2.sh id.rare.dnv_pt_fa_mo.txt $VCFRAREDNV ../id_snpeff.txt $PED2 denovofilterFormat #=> denovofilterFormat.snv2,indel2.txt 
denovofilterFormat3.R denovofilterFormat.snv2.txt denovofilterFormat.indel2.txt               id.rare.dnv_df.txt 
merge_denovoqc.R \
    id.rare.dnv_pt_fa_mo.txt \
    id.rare.dnv_vq.txt \
    id.rare.dnv_td.txt \
    id.rare.dnv_dn.txt \
    id.rare.dnv_dg.txt \
    id.rare.dnv_df.txt \
    $SANGER \
    ../id_snpeff2.txt \ #normalized no hougaii??
    ../id_annovar.txt \
    $OUTVARIANTQC 

#snv_vq: orst: -8.178; 2nd: -3.941; my: -28.3  
#snv_td: worst: 5.22;   2nd: 5.72               
#snv_dn: worst: 0.1994; 2nd: 0.46948            
#snv_dg: worst: 0.0026; 2nd: 0.02               
#indel_vq: worst:-3.182;  2nd: -0.7775; my: -5.77 
#indel_td: worst:5.52;    2nd: 6.55               
check_denovoqc_by_plot.R \
    --denovoqc $OUTVARIANTQC \
    --ped $PED2 \
    --sampleqc ../sampleqc.20191202.txt \
    --snv_vq -8.18 \
    --snv_td 5.2 \
    --snv_dn 0.199 \
    --snv_dg 0.0025 \
    --indel_vq -3.19 \
    --indel_td 5.5 \
    --indel_df \
    --filteredtotal_max 10 \
    --out_png $OUTPNG \
    --out_txt $OUTSAMPLEQC 



## select damaging sv
#varid_oi2 = []
#for idx,row in ped.iterrows():
#  pt = row['pt']
#  print(pt)
#  annotsv_path = '/betelgeuse07/analysis/hamanaka/wgs_sv/manta/annotsv/'+pt+'/annotsv.'+pt+'.final.MafFiltered.tsv'
#  try:
#    annotsv = pd.read_table(annotsv_path,index_col=0).query('Annotation_mode=="split"')
#    varids_OneSample = annotsv.apply(lambda se: extract_deleterious_sv(se), axis=1).to_list()
#    varid_oi2 = varid_oi2+varids_OneSample
#  except:
#    pass
#
#varid_oi2 = set(varid_oi2)
#varid_oi2.remove(None)
######################
