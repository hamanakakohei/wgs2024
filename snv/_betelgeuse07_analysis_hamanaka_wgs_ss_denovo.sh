#VCFPATH=/betelgeuse07/analysis/ncgm/FY2020-2021_G011/
VCFPATH=/betelgeuse07/analysis/hamanaka/wgs_40/
REF=/mira05/analysis/hamanaka/resource/ref/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
BCFTOOLS="/usr/local/genome/bcftools-1.8/bin/bcftools"
VCFTOOLS="/usr/local/genome/vcftools-0.1.16/bin/vcftools"
PEDTRIO=/betelgeuse07/analysis/hamanaka/wgs/sample1211.plus40trio.ped
HEALTHY=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.txt
SAMPLE_BAMPATH=/betelgeuse07/analysis/hamanaka/wgs/sample_BamPath.txt
SELECT="java -jar /usr/local/bio/src/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants"
QSUBSPAI=/mira05/analysis/hamanaka/resource/qsub/qsub__spai.sh
export PATH=/usr/local/genome/samtools-1.9/bin/:$PATH
export PATH=/antares01/analysis/hamanaka/function/:$PATH
source Variable.sh
source VcfFunctions.sh
UCHIYAMAANNOVAR=/betelgeuse01/analysis/uchiyama/annovar_20190508/annovarHg38_20210414.sh
#DDD=/betelgeuse04/analysis/hamanaka/resource/DDG2P_14_10_2021.tsv
DDD=/betelgeuse04/analysis/hamanaka/resource/DDG2P_20_3_2023.tsv
OMIM=/betelgeuse07/analysis/hamanaka/resource/omim.genemap2.edit.txt
PLI=/betelgeuse04/analysis/hamanaka/resource/gnomad.v2.1.1.lof_metrics.by_transcript.txt
SUPANNOVA=/betelgeuse07/analysis/hamanaka/wgs741/supp_AnnovaFile.txt
ONEPROBAND=/betelgeuse07/analysis/hamanaka/wgs/one_proband_each_fam.txt

#for CHR in `seq 20 -1 13` X; do
for CHR in 22; do
  TOMMO=/mira05/analysis/hamanaka/resource/tommo/chr$CHR.tommo.hg38.sort.txt
  BRAVO=/mira05/analysis/hamanaka/resource/topmed/chr$CHR.TOPMed_freeze_8.all.sort.txt
  NCBN=/mira05/analysis/hamanaka/resource/ncbn/ncbn.slim.$CHR.sort.txt
  NCGM=/betelgeuse07/analysis/hamanaka/wgs/$CHR/id_af.$CHR.sort.txt
  YCU=/betelgeuse07/analysis/hamanaka/wgs_ycu/$CHR/id_af.$CHR.sort.txt
  WBBC=/mira05/analysis/hamanaka/resource/wbbc/WBBC.chr$CHR.GRCh38.sort.txt
  SPAI=/mira05/analysis/hamanaka/resource/spliceai/SnvInd.$CHR.txt
  mkdir $CHR; cd $CHR 
  $BCFTOOLS norm -m -any -f $REF ${VCFPATH}combined120.$CHR.vcf.gz|bgzip -c > jointcall.VQSR.chr$CHR.n.vcf.gz

  # kokode mouippatu norm sitahougaii kedo hotonndo chigawanaikara sabottemoii!!!!!!!!!
  awk '$1 ~ /^#/ || ($1 !~ /^#/ && $5!="*"){print $0}' <(less jointcall.VQSR.chr$CHR.n.vcf.gz)|bgzip -c > jointcall.VQSR.chr$CHR.n2.vcf.gz
  tabix jointcall.VQSR.chr$CHR.n2.vcf.gz
  $SELECT -R $REF -sf $HEALTHY    -V jointcall.VQSR.chr$CHR.n2.vcf.gz -o jointcall.VQSR.chr$CHR.n2Healthy.vcf.gz > log.GatkSelect.txt 2>&1
  extract_af jointcall.VQSR.chr$CHR.n2Healthy.vcf.gz|sed 's/AF\=//g'|sort > id_af.$CHR.sort.txt
 
 # add precomputed annotations & MAFs
  gzip -dc jointcall.VQSR.chr$CHR.n2.vcf.gz|cut -f1-10 > riker.private.novo.vcf
  bash $UCHIYAMAANNOVAR > log.annovar.$CHR.txt 2>&1
  paste exome_test_refGene.hg38_multianno.txt $SUPANNOVA|cut -f1-161 > exome_test_refGene.hg38_multianno.2.txt
  source /usr/local/genome/python3-venv/env.sh
  add_info_exome_summary.py -annovar exome_test_refGene.hg38_multianno.2.txt -spliceai $SPAI -topmed $BRAVO -tommo $TOMMO -wbbc $WBBC -ncbn $NCBN -ncgm $NCGM -ycu $YCU -out exome_test_refGene.hg38_multianno.edit.txt    
  
# extract very rare
  filter_maf.py -annovar exome_test_refGene.hg38_multianno.edit.txt -out id.rare.$CHR.txt
  
# extract denovo
  $VCFTOOLS --gzvcf jointcall.VQSR.chr$CHR.n2.vcf.gz --mendel $PEDTRIO --out id.fam.$CHR > log.mendel.$CHR.txt 2>&1
  awk '($6=="0\/1" || $6=="1\/0" || $6=="1\/1") && $7=="0\/0" && $8=="0\/0"{print $1"-"$2"-"$3"-"$4"\t"$5}' id.fam.$CHR.mendel| awk -F"_" '{print $1"\t"$2"\t"$3}'|awk '$1 !~ /*/'|sort > id.dnv_pt_fa_mo.$CHR.txt
 
 # ids of very rare denovo
  join -t$'\t' <(sort id.rare.$CHR.txt) <(sort id.dnv_pt_fa_mo.$CHR.txt) > id.rare.dnv_pt_fa_mo.$CHR.txt
  cut -f1 id.rare.dnv_pt_fa_mo.$CHR.txt|sort -u > id.rare.dnv.$CHR.txt
  
# for annotations computed by ourseleves
  subset_vcf_by_variant jointcall.VQSR.chr$CHR.n2.vcf.gz id.rare.dnv.$CHR.txt > jointcall.VQSR.chr$CHR.n2RareDnv.vcf 
  qsub -N spai$CHR $QSUBSPAI jointcall.VQSR.chr$CHR.n2RareDnv.vcf $REF jointcall.VQSR.chr$CHR.n2RareDnvSpai.vcf 300
  sleep 20s
  while [ `qstat|grep spai$CHR|wc -l` != 0 ]; do
    echo 'waiting'
    sleep 30m
  done
  extract_spliceai	jointcall.VQSR.chr$CHR.n2RareDnvSpai.vcf      > id.rare.dnv_spai.txt
 
 # add our annotations to annovar res
  source /usr/local/genome/gcc-8.5.0/env.sh
  merge_denovoqc3.R \
    id.rare.dnv_pt_fa_mo.$CHR.txt \
    id.rare.dnv_spai.txt \
    denovoqc.$CHR.txt 
  annotate_denovo.R denovoqc.$CHR.txt exome_test_refGene.hg38_multianno.edit.txt denovoqc.anno.$CHR.txt
  
# extract genotype
  awk 'BEGIN{OFS="\t"}NR>1{split($1,A,"_"); print $2,A[1],A[2],A[3],A[4]; print $3,A[1],A[2],A[3],A[4]; print $4,A[1],A[2],A[3],A[4]}' denovoqc.anno.$CHR.txt|sort -u > sample_chr_pos_ref_alt.txt
  extract_genotype_from_variant_sample_combi2 sample_chr_pos_ref_alt.txt jointcall.VQSR.chr$CHR.n2.vcf.gz > sample_varid_gt.txt
 
 # add genotype + ddd/omim/pli
  annotate_denovo3.R denovoqc.anno.$CHR.txt sample_varid_gt.txt $DDD $OMIM $PLI denovoqc.anno2.$CHR.txt
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

