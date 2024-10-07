REF=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
PED=sample1211.plus40trio.ped
SELECT="java -jar GenomeAnalysisTK.jar -T SelectVariants"
source VcfFunctions.sh
DDD=DDG2P_14_10_2021.tsv
OMIM=omim.genemap2.edit.txt
PLI=gnomad.v2.1.1.lof_metrics.by_transcript.txt
PEDREFORMAT=sample1211reformat.plus40trio.ped
ONEPROBAND=one_proband_each_fam2.txt

for CHR in `seq 1 22` X ; do
  cd $CHR 
  # extract maf1%
  source python3-venv/env.sh
  filter_maf.py \
    -annovar exome_test_refGene.hg38_multianno.edit.txt \
    -esp6500siv2_all 0.01 \
    -gnomAD_exome_ALL 0.01 \
    -gnomAD_exome_AFR 0.01 \
    -gnomAD_exome_AMR 0.01 \
    -gnomAD_exome_ASJ 0.01 \
    -gnomAD_exome_EAS 0.01 \
    -gnomAD_exome_FIN 0.01 \
    -gnomAD_exome_NFE 0.01 \
    -gnomAD_exome_OTH 0.01 \
    -gnomAD_exome_SAS 0.01 \
    -AF 0.01 \
    -AF_afr 0.01 \
    -AF_ami 0.01 \
    -AF_amr 0.01 \
    -AF_asj 0.01 \
    -AF_eas 0.01 \
    -AF_fin 0.01 \
    -AF_nfe 0.01 \
    -AF_oth 0.01 \
    -AF_sas 0.01 \
    -topmed 0.01 \
    -tommo 0.01 \
    -wbbc 0.01 \
    -ncbn_all 0.01 \
    -ncbn_hondo 0.01 \
    -ncbn_ryukyu 0.01 \
    -ncgm 0.01 \
    -out id.maf0.01.txt
  
  # extract damaging
  select_damaging_vars.py -annovar exome_test_refGene.hg38_multianno.edit.txt -out id.damaging.txt
  cat \
    <(grep "${CHR}-" ../varid__3utr.ncgm.txt) \
    <(grep "${CHR}-" ../varid__GeneUtrannotator.txt) \
    id.damaging.txt \
    > id.damaging2.txt
  
  # extract maf1% dmaging
  join <(sort id.maf0.01.txt) <(sort id.damaging2.txt) > id.maf0.01.damaging.txt
  
  # add genotype for maf1% damaging
  $SELECT -R $REF -sf $ONEPROBAND -V jointcall.VQSR.chr$CHR.n2.vcf.gz -o jointcall.VQSR.chr$CHR.n2OneProband.vcf.gz --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES > log.GatkSelect2.txt 2>&1
  extract_existing_variant_id jointcall.VQSR.chr$CHR.n2OneProband.vcf.gz > id_OneProband.txt
  join <(sort id.maf0.01.damaging.txt) <(sort id_OneProband.txt) > id.maf0.01.damaging_OneProband.txt
  join -1 2 \
    <(sort -k2,2 id.maf0.01.damaging_OneProband.txt) \
    <(cut -f2- $PEDREFORMAT | sort) \
    | awk 'BEGIN{OFS="\t"}{split($2,A,"-"); print $1,A[1],A[2],A[3],A[4]; if(NF>2){for(i=3;i<=NF;i++)print $i,A[1],A[2],A[3],A[4]}}' \
    > sample_chr_pos_ref_alt__ar.txt
  extract_genotype_from_variant_sample_combi2 sample_chr_pos_ref_alt__ar.txt jointcall.VQSR.chr$CHR.n2.vcf.gz > sample_varid_gt__ar.txt
  
  source gcc-8.5.0/env.sh
  annotate_ar.R id.maf0.01.damaging_OneProband.txt $PED sample_varid_gt__ar.txt exome_test_refGene.hg38_multianno.edit.txt $DDD $OMIM $PLI ar_anno.txt
  
  cd ..
done 

