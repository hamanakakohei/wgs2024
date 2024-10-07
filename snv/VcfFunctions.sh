
extract_af(){ #1st arg: vcf; output: status"\t"vqslod score
    awk -F"\t" '$1 !~ /^#/{
        n=split($8,A,";")
        for(i=1;i<=n;i++){
            if(A[i] ~ /^AF=/){AF=A[i]; break}
            if(i == n){AF="NoScore"}
        }
        print $1"-"$2"-"$4"-"$5"\t"AF
    }' <(zcat $1)
}

extract_spliceai(){ #1st arg: vcf; output: status"\t"vqslod score
    awk -F"\t" '$1 !~ /^#/{
        printf $1"\t"$2"\t"$4"\t"$5"\t"
        n=split($8,A,";")
        for(i=1;i<=n;i++){
            if(A[i] ~ /^SpliceAI=/){print A[i]; break}
            if(i == n){print "NoScore"}
        }
    }' <(zcat $1) | sed -e 's/SpliceAI=//g'
}

subset_vcf_by_variant(){
    #1st arg: vcf
    #2nd arg: variant id(chr-pos-ref-alt
    if [[ $1 =~ gz$ ]]; then
        awk 'NR == FNR{VARID[$1]}NR != FNR && $1 ~ /^#/{print $0}NR != FNR && $1 !~ /^#/ && $1"-"$2"-"$4"-"$5 in VARID{print $0}' $2 <(gzip -dc $1) 
    else
        awk 'NR == FNR{VARID[$1]}NR != FNR && $1 ~ /^#/{print $0}NR != FNR && $1 !~ /^#/ && $1"-"$2"-"$4"-"$5 in VARID{print $0}' $2 $1 
    fi
}

extract_genotype_from_variant_sample_combi2(){
    #1st:sample,chr,pos,ref,alt(whose genotype you want to extract from vcf, no header); 2nd: vcf
    #output: the sample, the variant, the genotype
    awk -F"\t" 'NR==FNR{VARIANT[$2":"$3":"$4":"$5];SAMPLE_VARIANT[$1":"$2":"$3":"$4":"$5]}NR!=FNR && $1=="#CHROM"{
        for(i=10;i<=NF;i++){SAMPLE[i]=$i}
    }NR!=FNR && $1 !~ /^#/{
        if($1":"$2":"$4":"$5 in VARIANT){
            for(i=10;i<=NF;i++){
                if(SAMPLE[i]":"$1":"$2":"$4":"$5 in SAMPLE_VARIANT){
                    print SAMPLE[i],$1,$2,$4,$5,$i
                }
            }
        }
    }' $1 <(gzip -dc $2) > tmp${FUNCNAME[0]}.sample_chr_pos_ref_alt_gt.txt

    join -a 1 -o 1.1 2.2 -e "notexist" <(awk '{print $1":"$2":"$3":"$4":"$5}' $1|sort) <(awk '{print $1":"$2":"$3":"$4":"$5,$6}' tmp${FUNCNAME[0]}.sample_chr_pos_ref_alt_gt.txt|sort)|awk -F"[: ]" '{print $1,$2"_"$3"_"$4"_"$5,$0}'|cut -d" " -f1,2,4|tr ' ' '\t'
}

extract_existing_variant_id(){
  local VCF=$1
  awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^#CHROM/{
    for(i=10;i<=NF;i++){SAMPLE[i]=$i}
  }$1 !~ /^#/{
    for(i=10;i<=NF;i++){
      if($i~/^0\/1/ || $i~/^1\/0/ || $i~/^1\/1/)print $1"-"$2"-"$4"-"$5"\t"SAMPLE[i]
    }
  }' <(zcat $VCF)
}

extract_genotype_from_variant_sample_combi2(){
    #1st:sample,chr,pos,ref,alt(whose genotype you want to extract from vcf, no header); 2nd: vcf
    #output: the sample, the variant, the genotype
    awk -F"\t" 'NR==FNR{VARIANT[$2":"$3":"$4":"$5];SAMPLE_VARIANT[$1":"$2":"$3":"$4":"$5]}NR!=FNR && $1=="#CHROM"{
        for(i=10;i<=NF;i++){SAMPLE[i]=$i}
    }NR!=FNR && $1 !~ /^#/{
        if($1":"$2":"$4":"$5 in VARIANT){
            for(i=10;i<=NF;i++){
                if(SAMPLE[i]":"$1":"$2":"$4":"$5 in SAMPLE_VARIANT){
                    print SAMPLE[i],$1,$2,$4,$5,$i
                }
            }
        }
    }' $1 <(gzip -dc $2) \
    > tmp${FUNCNAME[0]}.sample_chr_pos_ref_alt_gt.txt

    join -a 1 -o 1.1 2.2 -e "notexist" \
      <(awk '{print $1":"$2":"$3":"$4":"$5}' $1 | sort) \
      <(awk '{print $1":"$2":"$3":"$4":"$5,$6}' tmp${FUNCNAME[0]}.sample_chr_pos_ref_alt_gt.txt | sort) \
      | awk -F"[: ]" '{print $1,$2"_"$3"_"$4"_"$5,$0}' \
      | cut -d" " -f1,2,4 \
      | tr ' ' '\t'
}

