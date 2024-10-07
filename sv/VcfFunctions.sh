select_samples(){ #1st:vcf; 2nd:samplelist
    awk 'BEGIN{OFS="\t"}NR==FNR{
            SAMPLEOI[$1]
    }NR!=FNR && $1 ~ /^#/ && $1!="#CHROM"{
        print $0
    }NR!=FNR && $1=="#CHROM"{
        printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
        for(i=10;i<=NF;i++)SAMPLE[i]=$i
        for(i=10;i<=NF;i++){
            if(SAMPLE[i] in SAMPLEOI){
                printf "\t"$i
            }
        }
        printf "\n"
    }NR!=FNR && $1 !~ /^#/{
        printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
        for(i=10;i<=NF;i++){
            if(SAMPLE[i] in SAMPLEOI){
                printf "\t"$i
            }
        }
        printf "\n"
    }' $2 <(zcat $1)
}

exclude_no_carrier_variant(){
  local VCF=$1
  awk '$1~/^#/{print $0}$1!~/^#/{for(i=10;i<=NF;i++){if($i~/0\/1/ || $i~/1\/0/ || $i~/1\/1/){print $0; next} ; if(i==NF){next}}}' <(less $VCF)
}

filter_by_supp(){
  local VCF=$1; local SUPP_THR=$2
  awk -v SUPP_THR=${SUPP_THR} '$1~/^#/{print $0}$1!~/^#/{
    N=split($8,A,";")
    for(i=1;i<=N;i++){
      if(A[i]~/SUPP=/){
        sub("SUPP=","",A[i])
        if((A[i]+0)<(SUPP_THR+0)){print $0}
      }
    }
  }' <(zcat $VCF)
}

split_manta_family_vcf(){
  local VCF=$1
  NCOL=`less $VCF|head -n5000|awk -F"\t" '$1=="#CHROM"{print NF}'`
  for i in `seq 10 $NCOL`; do
    SAMPLE=`less $VCF|head -n5000|awk -F"\t" 'BEGIN{OFS="\t"}$1=="#CHROM"'|cut -f $i`
    cat $VCF \
      | cut -f1-9,$i \
      | awk -v SAMPLE="${SAMPLE}" 'BEGIN{OFS="\t"
          }$1 ~ /^#/{
	    print $0
          }$1 !~ /^#/{
	    $3 = SAMPLE"_"NR
	    print $0
          }' > $SAMPLE.split.vcf
    exclude_no_carrier_variant $SAMPLE.split.vcf > $SAMPLE.split.2.vcf
  done
}

rename_id_and_sort(){
  local VCF=$1
  less $VCF \
    | awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^##/{print $0}$1 !~ /^#/{$3=SAMPLE"_"NR;print $0}$1 ~ /^#CHROM/{SAMPLE=$10;print $0}' \
    | vcf-sort
}

reformat_cnvnator_vcf(){
  local VCF=$1
  awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^##/{
    print $0
  }$1 == "#CHROM"{
    N = split($10, A, "\/")
    $10 = A[(N-2)]
    print $0
  }$1 !~ /^#/ && $1 !~ /^HLA/ && $1 !~ /_/ && $1 !~ /EBV/{
    for(i=10;i<=NF;i++){
      if($i ~ /^\.\/1/){
        sub("^\./1", "0/1", $i)
      }
    }
    print $0
  }' <(zcat $VCF)
}

reformat_canvas_vcf(){
  #1: CNVLEN -> SVLEN
  #2: SVTYPE no CNV -> DEL or DUP
  #3: ALT col ga "," wo fukundehainai
  #4: ID col ga "COMPLEXCNV" wo fukundehainai
  local VCF=$1
  awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^#/{
    print $0
  }$1 !~ /^#/ && $3 ~ /LOSS/ && $5 !~ /,/{
    sub("CNVLEN", "SVLEN", $8)
    sub("SVTYPE=CNV", "SVTYPE=DEL", $8)
    for(i=10;i<=NF;i++){
      if($i ~ /^\.\/1/){
        sub("^\./1", "0/1", $i)
      }
    }
    print $0
  }$1 !~ /^#/ && $3 ~ /GAIN/ && $5 !~ /,/{
    sub("CNVLEN", "SVLEN", $8)
    sub("SVTYPE=CNV", "SVTYPE=DUP", $8)
    for(i=10;i<=NF;i++){
      if($i ~ /^\.\/1/){
        sub("^\./1", "0/1", $i)
      }
    }
    print $0
  }' <(less $VCF)
}

split_canvas_family_vcf(){
  local VCF=$1
  NCOL=`less $VCF|head -n5000|awk -F"\t" '$1=="#CHROM"{print NF}'`
  for i in `seq 10 $NCOL`; do
    SAMPLE=`less $VCF|head -n5000|awk -F"\t" 'BEGIN{OFS="\t"}$1=="#CHROM"'|cut -f $i`
    echo $SAMPLE
    less $VCF|cut -f1-9,$i | awk -v SAMPLE="${SAMPLE}" 'BEGIN{OFS="\t"}$1 ~ /^#/{print $0}$1 !~ /^#/{$3 = SAMPLE"_"NR; print $0}' > $SAMPLE.vcf
    exclude_no_carrier_variant3 $SAMPLE.vcf > $SAMPLE.2.vcf
  done
}
