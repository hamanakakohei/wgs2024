POLYA=/mira05/analysis/hamanaka/resource/gencode.v41.metadata.PolyA_feature.gz 
SiteVcf=G011_jointcall.siteonly.n2.vcf.gz
SiteVcfX=jointcall.VQSR.chrX.n2.site.vcf
SiteVcfNcgm=G011_jointcall.siteonly.n2.PlusChrX.vcf
SiteVcfYcu=../wgs_ycu/jointcall.VQSR.site.n2.vcf.gz
BEDTOOLS=/usr/local/genome/bedtools2-2.26.0/bin/bedtools

less $POLYA|cut -f4,5,6,8|awk '$1~"^chr"'|sed 's/chr//' > tmp.3utr.txt
grep -v -e X -e Y tmp.3utr.txt|sort -k1,1n -k2,2n|awk '{print "chr"$0}'  > tmp.polya.bed
grep    -e X -e Y tmp.3utr.txt|sort -k1,1  -k2,2n|awk '{print "chr"$0}' >> tmp.polya.bed
awk '!a[$0]++' tmp.polya.bed > tmp.polya2.bed
awk '{print $1"\t"$2-1"\t"$3"\t"$4}' tmp.polya2.bed > tmp.polya3.bed

$BEDTOOLS intersect -a $SiteVcfNcgm -b tmp.polya3.bed  -wa -wb|awk '{print $1"_"$2"_"$4"_"$5}' > varid__3utr.ncgm.txt
$BEDTOOLS intersect -a $SiteVcfYcu  -b tmp.polya3.bed  -wa -wb|awk '{print $1"_"$2"_"$4"_"$5}' > varid__3utr.ycu.txt

#grep -f <(awk '{print $1"_"$2"_"$4"_"$5}' tmp.txt)     denovoqc.anno.all.txt
#grep -f <(awk '{print $1"_"$2"_"$4"_"$5}' tmp_ycu.txt) ../wgs_ycu/denovoqc.anno.all.YcuAf0.008.txt

