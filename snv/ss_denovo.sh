REF=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
BCFTOOLS=bcftools
VCFTOOLS=vcftools
PEDTRIO=sample1211.plus40trio.ped
HEALTHY=ncgmid__unaffected844.txt
SAMPLE_BAMPATH=sample_BamPath.txt
SELECT="java -jar GenomeAnalysisTK.jar -T SelectVariants"
QSUBSPAI=qsub__spai.sh
source VcfFunctions.sh
UCHIYAMAANNOVAR=annovar_20190508/annovarHg38_20210414.sh
DDD=DDG2P_20_3_2023.tsv
OMIM=omim.genemap2.edit.txt
PLI=gnomad.v2.1.1.lof_metrics.by_transcript.txt
SUPANNOVA=supp_AnnovaFile.txt
ONEPROBAND=one_proband_each_fam.txt

for CHR in `seq 1 22` X; do
  TOMMO=tommo/chr$CHR.tommo.hg38.sort.txt
  BRAVO=topmed/chr$CHR.TOPMed_freeze_8.all.sort.txt
  NCBN=ncbn/ncbn.slim.$CHR.sort.txt
  NCGM=wgs/$CHR/id_af.$CHR.sort.txt
  YCU=wgs_ycu/$CHR/id_af.$CHR.sort.txt
  WBBC=wbbc/WBBC.chr$CHR.GRCh38.sort.txt
  SPAI=spliceai/SnvInd.$CHR.txt
  mkdir $CHR; cd $CHR 
  $BCFTOOLS norm -m -any -f $REF ${VCFPATH}combined120.$CHR.vcf.gz | bgzip -c > jointcall.VQSR.chr$CHR.n.vcf.gz

  awk '$1 ~ /^#/ || ($1 !~ /^#/ && $5!="*"){print $0}' <(less jointcall.VQSR.chr$CHR.n.vcf.gz) | bgzip -c > jointcall.VQSR.chr$CHR.n2.vcf.gz
  tabix jointcall.VQSR.chr$CHR.n2.vcf.gz
  $SELECT -R $REF -sf $HEALTHY    -V jointcall.VQSR.chr$CHR.n2.vcf.gz -o jointcall.VQSR.chr$CHR.n2Healthy.vcf.gz > log.GatkSelect.txt 2>&1
  extract_af jointcall.VQSR.chr$CHR.n2Healthy.vcf.gz | sed 's/AF\=//g' | sort > id_af.$CHR.sort.txt
 
 # add precomputed annotations & MAFs
  gzip -dc jointcall.VQSR.chr$CHR.n2.vcf.gz | cut -f1-10 > riker.private.novo.vcf
  bash $UCHIYAMAANNOVAR > log.annovar.$CHR.txt 2>&1
  paste exome_test_refGene.hg38_multianno.txt $SUPANNOVA | cut -f1-161 > exome_test_refGene.hg38_multianno.2.txt
  source python3-venv/env.sh
  add_info_exome_summary.py -annovar exome_test_refGene.hg38_multianno.2.txt -spliceai $SPAI -topmed $BRAVO -tommo $TOMMO -wbbc $WBBC -ncbn $NCBN -ncgm $NCGM -ycu $YCU -out exome_test_refGene.hg38_multianno.edit.txt    
  
  # extract very rare
  filter_maf.py -annovar exome_test_refGene.hg38_multianno.edit.txt -out id.rare.$CHR.txt
  
  # extract denovo
  $VCFTOOLS --gzvcf jointcall.VQSR.chr$CHR.n2.vcf.gz --mendel $PEDTRIO --out id.fam.$CHR > log.mendel.$CHR.txt 2>&1
  awk '($6=="0\/1" || $6=="1\/0" || $6=="1\/1") && $7=="0\/0" && $8=="0\/0"{print $1"-"$2"-"$3"-"$4"\t"$5}' id.fam.$CHR.mendel \
    | awk -F"_" '{print $1"\t"$2"\t"$3}' \
    | awk '$1 !~ /*/' \
    | sort \
    > id.dnv_pt_fa_mo.$CHR.txt
 
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
  source gcc-8.5.0/env.sh
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


