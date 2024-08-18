EHDN=/usr/local/genome/ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo 
REF__gtexhg38=/antares01/analysis/hamanaka/resource/gtex_ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta 
#REF__1kghg38=/betelgeuse04/analysis/hamanaka/resource/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa 
#REF__YcuHg38=/antares01/analysis/hamanaka/wgs/ref__Mira04Seyama/Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix.fa
source /betelgeuse01/analysis/hamanaka/function/Variable.sh

# ehdn
do_ehdn_profile(){
    EHDN=/usr/local/genome/ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo
    BAMPATH=$1; REF=$2; WESWGS=$3
    #SAMPLE=`basename $BAMPATH|cut -d"." -f1`
    SAMPLE=`basename $BAMPATH`
    $EHDN profile --log-reads --reads $BAMPATH --reference $REF --output-prefix $SAMPLE.$WESWGS.ehdn --min-anchor-mapq 50 --max-irr-mapq 40 > log.$SAMPLE.$WESWGS.ehdn.txt 2>&1
    mkdir $SAMPLE.$WESWGS.ehdn; mv $SAMPLE.$WESWGS.ehdn.str_profile.json $SAMPLE.$WESWGS.ehdn.locus.tsv $SAMPLE.$WESWGS.ehdn.reads.tsv $SAMPLE.$WESWGS.ehdn.motif.tsv log.$SAMPLE.$WESWGS.ehdn.txt $SAMPLE.$WESWGS.ehdn/
}

less mira03_cram.txt|tail -n+4|while read BAMPATH;  do do_ehdn_profile $BAMPATH  $REF__gtexhg38 wgs; done 

# qc (locus file no gyousuu
5203 DA0000000479: epi fa
5220 DA0000000397: DD? pt
5224 DA0000000447: epi mo
5631 DA0000000604: taijisuishu pt
5633 DA0000002148: coffin pt
5638 DA0000000517: epi pt
5794 DA0000002267: epi mo
5838 DA0000000515: epi fa
5920 DA0000000646: tanshinshitu pt
6152 DA0000000628: shounouzenkesson pt
9279 DA0000002674: aicardi mo


# ehdn merge -> maf filtering -> candidate listing
source /usr/local/genome/python3-venv/env.sh
OUTLIER=/betelgeuse04/analysis/hamanaka/resource/ExpansionHunterDenovo/scripts/outlier.py
ANNOTATE=/betelgeuse04/analysis/hamanaka/resource/ExpansionHunterDenovo/scripts/annotate_ehdn.sh
HELPER=/betelgeuse04/analysis/hamanaka/resource/Fazal2020Scripts/EHDn-v0.6.2_HelperScripts
#REF=/mira04/analysis/seyama/genotyping/gatk4-hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix.fa
REF=/betelgeuse04/analysis/hamanaka/resource/gtex/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
#MANIFEST=/betelgeuse07/analysis/hamanaka/ehdn_ycu/manifest.txt
MANIFEST=/betelgeuse07/analysis/hamanaka/ehdn/manifest1333.noDA0000002674.txt
#CONTROL=/mira05/analysis/hamanaka/svtool/svtools/156healthy.txt
CONTROL=/betelgeuse07/analysis/hamanaka/wgs/ncgmid__unaffected844.noDA0000002674.txt
#PEDREFORMAT=/mira05/analysis/hamanaka/svtool/svtools/WgsPedigree.txt
PEDREFORMAT=/betelgeuse07/analysis/hamanaka/wgs/sample1211reformat.plus40trio.noDA0000002674.ped 
PREFIX=sample1333
EHDN=/usr/local/genome/ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo
HumanDb=/betelgeuse01/system/db/annovar/201612/hg38
ANNOVAR=/usr/local/bio/src/annovar2016oct/annotate_variation.pl
$EHDN merge --reference $REF --manifest $MANIFEST --output-prefix $PREFIX
$OUTLIER locus --manifest $MANIFEST --multisample-profile $PREFIX.multisample_profile.json --output $PREFIX.outlier_locus.tsv
$OUTLIER motif --manifest $MANIFEST --multisample-profile $PREFIX.multisample_profile.json --output $PREFIX.outlier_motif.tsv
$ANNOTATE --ehdn-results $PREFIX.outlier_locus.tsv --ehdn-annotated-results $PREFIX.outlier_locus.annotated.tsv --annovar-annotate-variation $ANNOVAR --annovar-humandb $HumanDb --annovar-buildver hg38
less $PEDREFORMAT|while read LINE; do
  Proband=`echo $LINE|awk '{print $2}'`
  echo $LINE|awk '{for(i=3;i<=NF;i++)print $i}' > tmp.others.txt
  python3 /antares01/analysis/hamanaka/function/ehdn_summary.py -proband $Proband -others tmp.others.txt -controls $CONTROL -locus_summary $PREFIX.outlier_locus.annotated.tsv -motif_summary $PREFIX.outlier_motif.tsv -out $Proband.txt    
done


# annotate ddg2p & genomad
less $PEDREFORMAT|while read LINE; do
  Proband=`echo $LINE|awk '{print $2}'`
  /antares01/analysis/hamanaka/function/annotate_ehdn.R $Proband.txt $Proband.anno.txt
  awk -F"\t" 'NR==1 || ($10<2 && $11>5)' $Proband.anno.txt > $Proband.anno.filt.txt
  sleep 3s
done

