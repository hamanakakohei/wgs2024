HGMD=/betelgeuse01/analysis/resource/HGMD/2022.1/HGMD_Pro_2022.1_hg38.vcf.gz 
DDG2P=/betelgeuse04/analysis/hamanaka/resource/DDG2P_14_10_2021.tsv

less $HGMD|awk '$8~/CLASS=DM;/ && $8~/%3Ac.-/ && $8!~/SVLEN=/ && $8!~/%3Ac.-[0-9]*\+/ && $8!~/%3Ac.-[0-9]*\-/' > 5utr_897.vcf # dattakadouka kaitenakute atode souzoudekaita
grep -f <(cut -f1 $DDG2P|tail -n+2|sort -u|awk '{print "GENE="$0";"}') \
  <(grep -v -e uAUG -e uSTOP -e uFrameShift /betelgeuse04/analysis/UTR/UTRannotator/5utr_897_test_centauri03.vcf) \
  > 5utr_897__ddg2p.vcf

PLI=/betelgeuse04/analysis/hamanaka/resource/gnomad.v2.1.1.lof_metrics.by_transcript.txt
GENCODE=/mira05/analysis/hamanaka/gtex/gencode.v26.annotation.gtf 
TommoHg38=/mira05/analysis/hamanaka/resource/tommo/chrAll.tommo.hg38.vcf 
TommoPass=/mira05/analysis/hamanaka/resource/tommo/chrAll.tommo.hg38.Pass.AfMt0.001.vcf
TommoPassUtr=/mira05/analysis/hamanaka/resource/tommo/chrAll.tommo.hg38.Pass.AfMt0.001.utr.vcf
BED=utr__PliMt0.99.bed
BEDTOOLS=/usr/local/genome/bedtools2-2.26.0/bin/bedtools

awk -F"\t" '$22>0.99 && $3=="true"{print $2}' $PLI > enst__PliMt0.99.txt

awk 'NR==FNR{ENST[$1]}NR!=FNR && $1!~/^#/ && $3=="UTR"{
  gsub("\"","",$12)
  sub("\\..*","",$12)
  if($12 in ENST)print $1"\t"$4"\t"$5
}' enst__PliMt0.99.txt $GENCODE|sort -k1,1 -k2,2n > $BED

awk 'BEGIN{OFS="\t"}$1~/^#/{print $0}$1!~/^#/ && $7=="PASS"{sub("AF=","",$8);if(($8+0)>0.001)print $0}' $TommoHg38 > $TommoPass
$BEDTOOLS intersect -a $TommoPass -b $BED -wa -header > $TommoPassUtr 

# extract uAUG, uSTOP, uFrameShift
#VCF=/betelgeuse07/analysis/UTR/data/G011_jointcall.siteonly.n2.UTR.denovoqc.vcf
#VCF=/betelgeuse07/analysis/UTR/UTRannotator/denovoqc.anno.all.n2.UTR.vcf
VCF=/betelgeuse07/analysis/UTR/UTRannotator/YCU_jointcall.denovoqc.siteonly.n2.UTR.vcf
#DENOVO=/betelgeuse07/analysis/hamanaka/wgs/denovoqc.anno.all.txt
DENOVO=/betelgeuse07/analysis/hamanaka/wgs_ycu/denovoqc.anno.all.txt
grep -e uAUG -e uSTOP -e uFrameShift $VCF|awk 'BEGIN{print "varid\t5utr"}NR>1{print $1"_"$2"_"$4"_"$5"\t"$8}' > varid_5utrAno.txt
join -t$'\t' -1 1 -2 2 <(sort varid_5utrAno.txt) <(sort -k2,2 $DENOVO)|tac > denovoqc.anno.all.5utr_ycu.txt




 awk '$1!~/^#/ && $NF!="-"' /betelgeuse04/analysis/UTR/UTRannotator/G011_jointcall.noStar.chrX.vcf.gz_MANE > tmp.x.txt
 awk '$1!~/^#/ && $NF!="-"' /betelgeuse04/analysis/UTR/UTRannotator/G011_jointcall.siteonly.n2._MANE_centauri05.vcf > tmp.auto.txt
 paste <(awk -F"/" '{print $1}' tmp.auto.txt) <(cut -f 3,4,22 tmp.auto.txt) > tmp.auto2.txt
 paste <(awk -F"/" '{print $1}' tmp.x.txt)    <(cut -f 3,4,22 tmp.x.txt)    > tmp.x2.txt
 cat tmp.auto2.txt tmp.x2.txt|awk '{gsub("_","\t",$1); print $0}'  > tmp.auto.x.txt
 awk 'BEGIN{OFS="\t"
}$4=="-"{
  print $1,$2-1,"A"$3,"A",$1"_"$2"_"$3"_"$4,$(NF-1),$NF
  print $1,$2-1,"T"$3,"T",$1"_"$2"_"$3"_"$4,$(NF-1),$NF
  print $1,$2-1,"G"$3,"G",$1"_"$2"_"$3"_"$4,$(NF-1),$NF
  print $1,$2-1,"C"$3,"C",$1"_"$2"_"$3"_"$4,$(NF-1),$NF
}$3=="-"{
  print $1,$2-1,"A","A"$4,$1"_"$2"_"$3"_"$4,$(NF-1),$NF
  print $1,$2-1,"T","T"$4,$1"_"$2"_"$3"_"$4,$(NF-1),$NF
  print $1,$2-1,"G","G"$4,$1"_"$2"_"$3"_"$4,$(NF-1),$NF
  print $1,$2-1,"C","C"$4,$1"_"$2"_"$3"_"$4,$(NF-1),$NF
}$3!="-" && $4!="-"{
  print $1,$2,$3,$4,$1"_"$2"_"$3"_"$4,$(NF-1),$NF
}' tmp.auto.x.txt > tmp.auto.x.2.txt

 sort -u tmp.auto.x.2.txt > tmp.auto.x.3.txt
 join <(sort tmp.varid.txt) <(awk '{print $1"_"$2"_"$3"_"$4"\t"$0}' tmp.auto.x.3.txt|sort) > tmp.auto.x.4.txt
 join -1 2 -2 1 <(awk '{print $1,$7,$8}' tmp.auto.x.4.txt|sort -k2,2) <(awk -F"[.\t]" '{print $1,$3}' /antares01/analysis/hamanaka/str/ensg_gene.txt|sort)|awk '{print $2,$4,$3}' > tmp.auto.x.5.txt
 awk '{A[$1" "$2]=A[$1" "$2]";"$3}END{for(i in A)print i,A[i]}' tmp.auto.x.5.txt|tr ' ' '\t' > varid_gene_utrannotator.txt



## select PASS & MAF0.1% from gnomad
#source /usr/local/genome/hail-0.2.34/env.sh
#python
#
#import os
#os.environ['PYSPARK_SUBMIT_ARGS'] = '--jars \
# /usr/local/miniconda3/lib/python3.7/site-packages/hail/hail-all-spark.jar \
# --conf spark.executor.cores=30 \
# --conf spark.executor.instances=4 \
# --conf spark.driver.memory=20G \
# --conf spark.executor.memory=300G \
#pyspark-shell '
#
#import hail as hl
#hl.init()
#
#gnomad211 = "/mira05/analysis/hamanaka/gnomad/gnomad.genomes.r2.1.1.sites.vcf.bgz"
#mt = hl.import_vcf(gnomad211)
#ht = mt.select_rows(mt.filters,mt.info.AF).rows()
#ht = ht.explode(ht.AF)
#ht = ht.filter(ht.filters==hl.empty_set(hl.tstr))
#ht = ht.filter(ht.AF>0.001)
#ht = ht.key_by()
#ht.select(ht.locus,ht.alleles).export("gnomad.genomes.r2.1.1.sites.pass.AfMt0.001.tsv",delimiter="\t")
#
## lift over tp hg38
#CHAIN=/betelgeuse04/analysis/hamanaka/resource/hg19ToHg38.over.chain.gz
#REF=/betelgeuse04/analysis/hamanaka/resource/hg38.fa
#GatkLift="/usr/local/genome/gatk-4.2.2.0/gatk LiftoverVcf"
#VcfIn=gnomad.genomes.r2.1.1.sites.pass.AfMt0.001.vcf # kono vcf no header ha xxx kara tottekita
#VcfOut=gnomad.genomes.r2.1.1.sites.pass.AfMt0.001.Hg38.vcf
#
#$GatkLift I=$VcfIn O=$VcfOut REJECT=tmp.reject.txt CHAIN=$CHAIN R=$REF > log.txt 2>&1
