
a
a
a




# manta
ls ar_anno.plus_sv.DA*.omim.txt | while read FILE; do
  awk -F"\t" '
    ($4 ~ /biallelic/ || $6 ~ /recessive/) && \
    $(NF-188) == "ok" && \
    ($(NF-181) ~ /INV/ || ($(NF-181) ~ /DEL/ && ($(NF-179) + 0) < 100)){
      split($(NF-180), A, "-")
      if(A[1] == A[2] && A[1] ~ /intron/){print FILENAME, $0}
  }' $FILE
done

ls ar_anno.plus_sv.DA*.omim.txt | while read FILE; do
  awk -F"\t" '
    ($4 ~ /biallelic/ || $6 ~ /recessive/) && \
    $(NF-189) == "ok" && \
    ($(NF-182) ~ /INV/ || ($(NF-182) ~ /DEL/ && ($(NF-180) + 0) < 100)){
      split($(NF-181), A, "-")
      if(A[1] == A[2] && A[1] ~ /intron/){print FILENAME, $0}
  }' $FILE
done


$7 ~ /INV/ && $(NF-2) < 3{
    split($(NF-74), A, "-")
    if(A[1] == A[2] && A[1] ~ /intron/){print $2"\t"$8}
  }' $FILE
done | awk -F"[\t,]" '{print $1"\t"$NF}' > tmp.txt

cat <(echo -e 'proband\tSvId') tmp.txt > /betelgeuse05/analysis/WGS/manta_all_rare_and_SUPP2_intronic_inv/proband_svid.txt


## splicing intronic DEL near splice site
# manta
ls /betelgeuse07/analysis/hamanaka/wgs_sv/manta/annotsv/DA*/annotsv.DA*.final.Maf0Filtered.tsv|while read FILE; do
  awk -F"\t" '$7 ~ /DEL/ && $(NF-72) < 200 && $(NF-2) < 3{
    split($(NF-74), A, "-")
    if(A[1] == A[2] && A[1] ~ /intron/){print $2"\t"$8}
  }' $FILE
done | awk -F"[\t,]" '{print $1"\t"$NF}' > tmp.txt

cat <(echo -e 'proband\tSvId') tmp.txt > /betelgeuse05/analysis/WGS/manta_all_rare_and_SUPP2_ss_near_del/proband_svid.txt

