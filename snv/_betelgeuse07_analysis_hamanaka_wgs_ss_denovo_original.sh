VCFPATH=/betelgeuse07/analysis/ncgm/FY2020-2021_G011/
REF=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
#REF=/antares01/analysis/hamanaka/wgs/parliament/resources_broad_hg38_v0_Homo_sapiens_assembly38.PrimChrXYM.fasta
BCFTOOLS="/usr/local/genome/bcftools-1.8/bin/bcftools"
VCFTOOLS="/usr/local/genome/vcftools-0.1.16/bin/vcftools"
#PED=/betelgeuse07/analysis/hamanaka/wgs/NanbyoNcgmId20210910.ped
PED=/betelgeuse07/analysis/hamanaka/wgs/sample1211.plus40trio.ped
#PEDTRIO=/betelgeuse07/analysis/hamanaka/wgs/NanbyoNcgmId20210910Trio.ped
#PEDTRIO=/betelgeuse07/analysis/hamanaka/wgs/sample1192trio.ped
PEDTRIO=/betelgeuse07/analysis/hamanaka/wgs/sample1211trio.ped
#HEALTHY=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected751.txt
HEALTHY=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected765.txt
SAMPLE_BAMPATH=/betelgeuse07/analysis/hamanaka/wgs/sample_BamPath.txt
SELECT="java -jar /usr/local/bio/src/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants"
QSUBSPAI=/mira05/analysis/hamanaka/resource/qsub/qsub__spai.sh
export PATH=/usr/local/genome/samtools-1.9/bin/:$PATH
export PATH=/antares01/analysis/hamanaka/function/:$PATH
source Variable.sh
source VcfFunctions.sh
UCHIYAMAANNOVAR=/betelgeuse01/analysis/uchiyama/annovar_20190508/annovarHg38_20210414.sh
DDD=/betelgeuse04/analysis/hamanaka/resource/DDG2P_14_10_2021.tsv
OMIM=/betelgeuse07/analysis/hamanaka/resource/omim.genemap2.edit.txt
PLI=/betelgeuse04/analysis/hamanaka/resource/gnomad.v2.1.1.lof_metrics.by_transcript.txt
SUPANNOVA=/betelgeuse07/analysis/hamanaka/wgs741/supp_AnnovaFile.txt
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.ped
ONEPROBAND=/betelgeuse07/analysis/hamanaka/wgs/one_proband_each_fam2.txt

for CHR in `seq 22 -1 1` X ; do
  TOMMO=/mira05/analysis/hamanaka/resource/tommo/chr$CHR.tommo.hg38.sort.txt
  BRAVO=/mira05/analysis/hamanaka/resource/topmed/chr$CHR.TOPMed_freeze_8.all.sort.txt
  NCBN=/mira05/analysis/hamanaka/resource/ncbn/ncbn.slim.$CHR.sort.txt
  NCGM=/betelgeuse07/analysis/hamanaka/wgs/$CHR/id_af.$CHR.sort.txt
  #NCGMALL=/betelgeuse07/analysis/hamanaka/wgs/$CHR/id_af.$CHR.sort.all.txt
  YCU=/betelgeuse07/analysis/hamanaka/wgs_ycu/$CHR/id_af.$CHR.sort.txt
  WBBC=/mira05/analysis/hamanaka/resource/wbbc/WBBC.chr$CHR.GRCh38.sort.txt
  SPAI=/mira05/analysis/hamanaka/resource/spliceai/SnvInd.$CHR.txt
  mkdir $CHR; cd $CHR 
  #$BCFTOOLS norm -m -any -f $REF ${VCFPATH}G011_jointcall.VQSR.chr$CHR.vcf.gz|bgzip -c > jointcall.VQSR.chr$CHR.n.vcf.gz
  # kokode mouippatu norm sitahougaii kedo hotonndo chigawanaikara sabottemoii!!!!!!!!!
  #awk '$1 ~ /^#/ || ($1 !~ /^#/ && $5!="*"){print $0}' <(less jointcall.VQSR.chr$CHR.n.vcf.gz)|bgzip -c > jointcall.VQSR.chr$CHR.n2.vcf.gz
  #tabix jointcall.VQSR.chr$CHR.n2.vcf.gz
  #$SELECT -R $REF -sf $HEALTHY    -V jointcall.VQSR.chr$CHR.n2.vcf.gz -o jointcall.VQSR.chr$CHR.n2Healthy.vcf.gz > log.GatkSelect.txt 2>&1
  #extract_af jointcall.VQSR.chr$CHR.n2Healthy.vcf.gz|sed 's/AF\=//g'|sort > id_af.$CHR.sort.txt
  ##extract_af jointcall.VQSR.chr$CHR.n2.vcf.gz       |sed 's/AF\=//g'|sort > id_af.$CHR.sort.all.txt
  ## add precomputed annotations & MAFs
  #gzip -dc jointcall.VQSR.chr$CHR.n2.vcf.gz|cut -f1-10 > riker.private.novo.vcf
  #bash $UCHIYAMAANNOVAR > log.annovar.$CHR.txt 2>&1
  #paste exome_test_refGene.hg38_multianno.txt $SUPANNOVA|cut -f1-161 > exome_test_refGene.hg38_multianno.2.txt
  #source /usr/local/genome/python3-venv/env.sh
  #add_info_exome_summary.py -annovar exome_test_refGene.hg38_multianno.2.txt -spliceai $SPAI -topmed $BRAVO -tommo $TOMMO -wbbc $WBBC -ncbn $NCBN -ncgm $NCGM -ycu $YCU -out exome_test_refGene.hg38_multianno.edit.txt    
  ## extract very rare
  #filter_maf.py -annovar exome_test_refGene.hg38_multianno.edit.txt -out id.rare.$CHR.txt
  ## extract maf1%
  #filter_maf.py -annovar exome_test_refGene.hg38_multianno.edit.txt -esp6500siv2_all 0.01 -gnomAD_exome_ALL 0.01 -gnomAD_exome_AFR 0.01 -gnomAD_exome_AMR 0.01 -gnomAD_exome_ASJ 0.01 -gnomAD_exome_EAS 0.01 -gnomAD_exome_FIN 0.01 \
  #  -gnomAD_exome_NFE 0.01 -gnomAD_exome_OTH 0.01 -gnomAD_exome_SAS 0.01 -AF 0.01 -AF_afr 0.01 -AF_ami 0.01 -AF_amr 0.01 -AF_asj 0.01 -AF_eas 0.01 -AF_fin 0.01 -AF_nfe 0.01 -AF_oth 0.01 -AF_sas 0.01 -topmed 0.01 -tommo 0.01 \
  #  -wbbc 0.01 -ncbn_all 0.01 -ncbn_hondo 0.01 -ncbn_ryukyu 0.01 -ncgm 0.01 \
  #  -out id.maf0.01.txt
  ## extract maf1%
  #awk -F"-" '!(length($3)==1 && length($4)==2) && !(length($3)<6 && length($4)==1)' id.maf0.01.txt > id.spai.remain.txt
  ## extract denovo
  #$VCFTOOLS --gzvcf jointcall.VQSR.chr$CHR.n2.vcf.gz --mendel $PEDTRIO --out id.fam.$CHR > log.mendel.$CHR.txt 2>&1
  #awk '($6=="0\/1" || $6=="1\/0" || $6=="1\/1") && $7=="0\/0" && $8=="0\/0"{print $1"-"$2"-"$3"-"$4"\t"$5}' id.fam.$CHR.mendel| awk -F"_" '{print $1"\t"$2"\t"$3}'|awk '$1 !~ /*/'|sort > id.dnv_pt_fa_mo.$CHR.txt
  # extract damaging
  select_damaging_vars.py -annovar exome_test_refGene.hg38_multianno.edit.txt -out id.damaging.txt
  cat <(grep "${CHR}-" ../varid__3utr.ncgm.txt) <(grep "${CHR}-" ../varid__GeneUtrannotator.txt) id.damaging.txt > id.damaging2.txt
  # extract maf1% dmaging
  join <(sort id.maf0.01.txt) <(sort id.damaging2.txt) > id.maf0.01.damaging.txt
  ## add genotype for maf1% damaging
  #$SELECT -R $REF -sf $ONEPROBAND -V jointcall.VQSR.chr$CHR.n2.vcf.gz -o jointcall.VQSR.chr$CHR.n2OneProband.vcf.gz --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES > log.GatkSelect2.txt 2>&1
  #extract_existing_variant_id jointcall.VQSR.chr$CHR.n2OneProband.vcf.gz > id_OneProband.txt
  join <(sort id.maf0.01.damaging.txt) <(sort id_OneProband.txt) > id.maf0.01.damaging_OneProband.txt
  join -1 2 <(sort -k2,2 id.maf0.01.damaging_OneProband.txt) <(cut -f2- $PEDREFORMAT|sort)|awk 'BEGIN{OFS="\t"}{split($2,A,"-"); print $1,A[1],A[2],A[3],A[4]; if(NF>2){for(i=3;i<=NF;i++)print $i,A[1],A[2],A[3],A[4]}}' > sample_chr_pos_ref_alt__ar.txt
  extract_genotype_from_variant_sample_combi2 sample_chr_pos_ref_alt__ar.txt jointcall.VQSR.chr$CHR.n2.vcf.gz > sample_varid_gt__ar.txt
  source /usr/local/genome/gcc-8.5.0/env.sh
  annotate_ar.R id.maf0.01.damaging_OneProband.txt $PED sample_varid_gt__ar.txt exome_test_refGene.hg38_multianno.edit.txt $DDD $OMIM $PLI ar_anno.txt
  ## ids of very rare denovo
  #join -t$'\t' <(sort id.rare.$CHR.txt) <(sort id.dnv_pt_fa_mo.$CHR.txt) > id.rare.dnv_pt_fa_mo.$CHR.txt
  #cut -f1 id.rare.dnv_pt_fa_mo.$CHR.txt|sort -u > id.rare.dnv.$CHR.txt
  ## for annotations computed by ourseleves
  #subset_vcf_by_variant jointcall.VQSR.chr$CHR.n2.vcf.gz id.rare.dnv.$CHR.txt > jointcall.VQSR.chr$CHR.n2RareDnv.vcf 
  #qsub -N spai$CHR $QSUBSPAI jointcall.VQSR.chr$CHR.n2RareDnv.vcf $REF jointcall.VQSR.chr$CHR.n2RareDnvSpai.vcf 300
  #sleep 20s
  #while [ `qstat|grep spai$CHR|wc -l` != 0 ]; do
  #  echo 'waiting'
  #  sleep 30m
  #done
  #extract_spliceai	jointcall.VQSR.chr$CHR.n2RareDnvSpai.vcf      > id.rare.dnv_spai.txt
  ## add our annotations to annovar res
  #source /usr/local/genome/gcc-8.5.0/env.sh
  #merge_denovoqc3.R \
  #  id.rare.dnv_pt_fa_mo.$CHR.txt \
  #  id.rare.dnv_spai.txt \
  #  denovoqc.$CHR.txt 
  #annotate_denovo.R denovoqc.$CHR.txt exome_test_refGene.hg38_multianno.edit.txt denovoqc.anno.$CHR.txt
  ## extract genotype
  #awk 'BEGIN{OFS="\t"}NR>1{split($1,A,"_"); print $2,A[1],A[2],A[3],A[4]; print $3,A[1],A[2],A[3],A[4]; print $4,A[1],A[2],A[3],A[4]}' denovoqc.anno.$CHR.txt|sort -u > sample_chr_pos_ref_alt.txt
  #extract_genotype_from_variant_sample_combi2 sample_chr_pos_ref_alt.txt jointcall.VQSR.chr$CHR.n2.vcf.gz > sample_varid_gt.txt
  ## add genotype + ddd/omim/pli
  #annotate_denovo3.R denovoqc.anno.$CHR.txt sample_varid_gt.txt $DDD $OMIM $PLI denovoqc.anno2.$CHR.txt
  cd ..
done 

# make normalized site vcf
SiteVcf=/betelgeuse07/analysis/ncgm/FY2020-2021_G011/G011_jointcall.siteonly.vcf.gz
$BCFTOOLS norm -m -any -f $REF $SiteVcf|bgzip -c > G011_jointcall.siteonly.n.vcf.gz
awk '$1 ~ /^#/ || ($1 !~ /^#/ && $5!="*"){print $0}' <(less G011_jointcall.siteonly.n.vcf.gz)|bgzip -c > G011_jointcall.siteonly.n2.vcf.gz

# spai 
SpaiThr=0.1
cat <(head -n1 1/denovoqc.anno2.1.txt) <(ls */denovoqc.anno2.*.txt|xargs -I{} tail -n+2 {}) > denovoqc.anno.all.txt
cat <(head -n1 1/ar_anno.txt)          <(ls */ar_anno.txt         |xargs -I{} tail -n+2 {}) > ar_anno.txt
awk -F"\t" -v SpaiThr=${SpaiThr} 'NR==1{print $0}$16=="intronic" && $10!="NoScore" && $10!="NA"{split($10,A,"|"); if(A[3]>SpaiThr || A[4]>SpaiThr || A[5]>SpaiThr || A[6]>SpaiThr)print $0}' denovoqc.anno.all.txt|less
#awk -F"\t" -v SpaiThr=${SpaiThr} 'NR==1{print $0}$13!="intronic" && $7!="NoScore" && $7!="NA"{split($7,A,"|"); if(A[3]>SpaiThr || A[4]>SpaiThr || A[5]>SpaiThr || A[6]>SpaiThr)print $0}' denovoqc.anno.all.txt|less
#awk -F"\t" -v SpaiThr=${SpaiThr} '($15=="exonic" || $15=="exonic;splicing" || $15=="intronic" || $15=="splicing") && $9!="NoScore" && $9!="NA"{split($9,A,"|"); if(A[3]>SpaiThr || A[4]>SpaiThr || A[5]>SpaiThr || A[6]>SpaiThr)print $0}' denovoqc.anno.all.txt|less

# WES low depth
select_low_dep_var.R --var denovoqc.anno.all.txt --pos_dep gnomad.exomes.r2.1.coverage__DepthLt7.tsv --out denovoqc.anno.all.DepthLt7.txt


#
import pandas as pd
import numpy as np

ped = pd.read_table('sample1192.ped',names=('fam','pt','fa','mo','sex','status')).query('fam==pt')
dt = pd.read_table('ar_anno.txt',low_memory=False)
# extract varid OI
varid_spliceai = dt[['varid','Func.refGene','spliceai']].query('`Func.refGene`.str.contains("UTR") or `Func.refGene`.str.contains("intronic")').copy()
varid_spliceai['spliceai'] = varid_spliceai['spliceai'].map(lambda x: extract_spliceai(x))
varid_oi = set(varid_spliceai.query('spliceai>0.15')['varid'].to_list())
##################
dt = dt.fillna({'allelic_requirement':'aaa'}).query('allelic_requirement.str.contains("biallelic")')[['Gene.refGene','varid','sample','membs','gts']] 
dt['gt_pt'] = dt.apply(lambda df: split_gts(df['sample'],df['membs'],df['gts'],ped,'pt'),axis=1)
dt['gt_fa'] = dt.apply(lambda df: split_gts(df['sample'],df['membs'],df['gts'],ped,'fa'),axis=1)
dt['gt_mo'] = dt.apply(lambda df: split_gts(df['sample'],df['membs'],df['gts'],ped,'mo'),axis=1)
dt = pd.concat([dt,annotsv_all])
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

dt_all = pd.concat([dt_SinglGene,dt_MultiGene2,annotsv_all])
sample_gene_ArJudg = dt_all.groupby(['Gene.refGene','sample']).apply(lambda df: judge_ar(   df                     )).reset_index().rename({0:'ar'   },axis=1)
sample_gene_varoi  = dt_all.groupby(['Gene.refGene','sample']).apply(lambda df: judge_varoi(df,varid_oi | varid_oi2)).reset_index().rename({0:'varoi'},axis=1)
pd.merge(pd.merge(dt_all,sample_gene_ArJudge),sample_gene_varoi).query('ar=="ok" & varoi=="ok"')

# prep annotsv
annotsv_all = pd.DataFrame(columns=['Gene.refGene','varid','sample','gt_pt','gt_fa','gt_mo'])
for idx,row in ped.iterrows():
  pt = row['pt']
  fa = row['fa']
  mo = row['mo']
  print(pt)
  annotsv_path = '/betelgeuse07/analysis/hamanaka/wgs_sv/manta/annotsv/'+pt+'/annotsv.'+pt+'.final.MafFiltered.tsv'
  annotsv = pd.read_table(annotsv_path,index_col=0).query('Annotation_mode=="split"').rename({pt:'gt_pt','AnnotSV_ID':'varid','Gene_name':'Gene.refGene'},axis=1)
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
  annotsv_all = pd.concat([annotsv_all,annotsv[['Gene.refGene','varid','sample','gt_pt','gt_fa','gt_mo']]])

annotsv_all['gt_pt'] = annotsv_all['gt_pt'].map(lambda x: correct_gt_annotsv(x))
annotsv_all['gt_fa'] = annotsv_all['gt_fa'].map(lambda x: correct_gt_annotsv(x))
annotsv_all['gt_mo'] = annotsv_all['gt_mo'].map(lambda x: correct_gt_annotsv(x))
annotsv_all['inhe' ] = annotsv_all.apply(lambda df: judge_inhe(df['gt_pt'],df['gt_fa'],df['gt_mo']       ),axis=1)
###############

# select damaging sv
varid_oi2 = []
for idx,row in ped.iterrows():
  pt = row['pt']
  print(pt)
  annotsv_path = '/betelgeuse07/analysis/hamanaka/wgs_sv/manta/annotsv/'+pt+'/annotsv.'+pt+'.final.MafFiltered.tsv'
  annotsv = pd.read_table(annotsv_path,index_col=0).query('Annotation_mode=="split"')
  varids_OneSample = annotsv.apply(lambda se: extract_deleterious_sv(se), axis=1).to_list()
  varid_oi2 = varid_oi2+varids_OneSample

varid_oi2 = set(varid_oi2)
varid_oi2.remove(None)
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
    return se['AnnotSV_ID']

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
  DA_of_who = ped.query('fam==@SAMPLE')[who].iloc[0]
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
