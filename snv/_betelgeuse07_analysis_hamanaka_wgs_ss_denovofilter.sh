# depth, allelic ratio, etc...
import pandas as pd
DENOVO = '/betelgeuse07/analysis/hamanaka/wgs/denovoqc.anno.all.txt'
OUT = '/betelgeuse07/analysis/hamanaka/wgs/denovoqc__depth_filtered.txt'

denovo = pd.read_table(DENOVO, low_memory = False)
denovo['GT']    = denovo['gt'   ].map(lambda x: x.split(':')[0])
denovo['GT_fa'] = denovo['gt.fa'].map(lambda x: x.split(':')[0])
denovo['GT_mo'] = denovo['gt.mo'].map(lambda x: x.split(':')[0])
denovo['AD']    = denovo['gt'   ].map(lambda x: x.split(':')[1])
denovo['AD_fa'] = denovo['gt.fa'].map(lambda x: x.split(':')[1])
denovo['AD_mo'] = denovo['gt.mo'].map(lambda x: x.split(':')[1])
denovo['AD_ref']    = denovo['AD'   ].map(lambda x: int(x.split(',')[0]))
denovo['AD_fa_ref'] = denovo['AD_fa'].map(lambda x: int(x.split(',')[0]))
denovo['AD_mo_ref'] = denovo['AD_mo'].map(lambda x: int(x.split(',')[0]))
denovo['AD_alt']    = denovo['AD'   ].map(lambda x: int(x.split(',')[1]))
denovo['AD_fa_alt'] = denovo['AD_fa'].map(lambda x: int(x.split(',')[1]))
denovo['AD_mo_alt'] = denovo['AD_mo'].map(lambda x: int(x.split(',')[1]))
denovo['DP']    = denovo['AD_ref'   ] + denovo['AD_alt']
denovo['DP_fa'] = denovo['AD_fa_ref'] + denovo['AD_fa_alt']
denovo['DP_mo'] = denovo['AD_mo_ref'] + denovo['AD_mo_alt']
denovo.\
  query('DP > 7 & DP_fa > 7 & DP_mo > 7').\
  query('AD_alt > 2').\
  query('AD_alt / DP > 0.2 & AD_alt / DP < 0.8').\
  query('AD_fa_alt < 3 & AD_mo_alt < 3').\
  query('AD_fa_alt==0 | AD_mo_alt==0')\
  [['sample', 'varid']].to_csv(OUT, sep = '\t', index = False)

# vqslod
for CHR in `seq 1 22` X; do
  gzip -dc $CHR/jointcall.VQSR.chr$CHR.n2.vcf.gz > jointcall.VQSR.chr$CHR.n2.vcf
  awk '{print $1"_"$2"_"$4"_"$5"\t"$7}' jointcall.VQSR.chr$CHR.n2.vcf >> varid_vqslod.txt 
done

# overlap problematic/repeat region
BEDS=(
  /betelgeuse04/analysis/hamanaka/resource/gtex/SimpleRepeatHg38Ucsc.prim.numericsort.bed
  /betelgeuse04/analysis/hamanaka/resource/gtex/RmskHg38Ucsc.prim.numericsort.bed 
  /antares01/analysis/hamanaka/resource/sd__hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/encode_blacklist.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/giab_align_problem.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/giab_allele_imbalance.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/giab_call_conflict.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/giab_genotype_conflict.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/giab_filter_sse.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/giab_quality_issue.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/ncbi_ngs_high_stringency.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/ncbi_ngs_low_stringency.hg38.prim.numericsort.bed
  /mira05/analysis/hamanaka/resource/ncbi_ngs_dead_zone.hg38.prim.numericsort.bed
)
SiteVcfNcgm=G011_jointcall.siteonly.n2.PlusChrX.vcf
SiteVcfYcu=../wgs_ycu/jointcall.VQSR.site.n2.vcf.gz
BEDTOOLS=/usr/local/genome/bedtools2-2.26.0/bin/bedtools

for BED in ${BEDS[@]}; do
  BEDHEAD=`echo $BED|awk -F"/" '{print $NF}'|sed 's/.hg38.prim.numericsort.bed$//'`
  $BEDTOOLS intersect -a $SiteVcfNcgm -b $BED -wa -wb|awk '{print $1"_"$2"_"$4"_"$5}' > varid__$BEDHEAD.txt
done

# select HQ denovo
comm -12 <(sort varid__denovoqc.anno.all.txt) <(sort varid__VqslodPass.txt)             > tmp0.txt
comm -23 <(sort tmp0.txt)                     <(sort varid__RmskHg38Ucsc.txt)           > tmp1.txt
comm -23 <(sort tmp1.txt)                     <(sort varid__SimpleRepeatHg38Ucsc.txt)   > tmp2.txt
comm -23 <(sort tmp2.txt)                     <(sort varid__sd.txt)                     > tmp3.txt
comm -23 <(sort tmp3.txt)                     <(sort varid__giab_align_problem.txt)     > tmp4.txt
comm -23 <(sort tmp4.txt)                     <(sort varid__giab_allele_imbalance.txt)  > tmp5.txt
comm -23 <(sort tmp5.txt)                     <(sort varid__giab_filter_sse.txt)        > tmp6.txt
comm -23 <(sort tmp6.txt)                     <(sort varid__giab_quality_issue.txt)     > tmp7.txt
comm -23 <(sort tmp7.txt)                     <(sort varid__giab_genotype_conflict.txt) > tmp8.txt
comm -23 <(sort tmp8.txt)                     <(sort varid__ncbi_ngs_dead_zone.txt)     > tmp9.txt
join <(sort tmp9.txt) <(awk '{print $2,$1}' denovoqc__depth_filtered.txt|sort)  > denovoqc__final_filtered.txt


comm -23 <(sort varid__denovoqc.anno.all.txt) <(sort varid__sd.txt)                     > tmp23.txt
comm -23 <(sort tmp23.txt)                     <(sort varid__giab_align_problem.txt)     > tmp24.txt
comm -23 <(sort tmp24.txt)                     <(sort varid__giab_allele_imbalance.txt)  > tmp25.txt
comm -23 <(sort tmp25.txt)                     <(sort varid__giab_filter_sse.txt)        > tmp26.txt
comm -23 <(sort tmp26.txt)                     <(sort varid__giab_quality_issue.txt)     > tmp27.txt
comm -23 <(sort tmp27.txt)                     <(sort varid__giab_genotype_conflict.txt) > tmp28.txt
comm -23 <(sort tmp28.txt)                     <(sort varid__ncbi_ngs_dead_zone.txt)     > tmp29.txt
