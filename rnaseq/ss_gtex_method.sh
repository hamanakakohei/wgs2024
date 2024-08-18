OH=150
star_index=star_index_oh$OH
path_to_data=/betelgeuse05/analysis/RNAseq_20230216
genes_gtf=gencode.v26.GRCh38.genes.gtf
genome_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
DOCK=mira:5000/broadinstitute/gtex_rnaseq:latest

# build the STAR index:
mkdir $star_index
docker run --rm -v $path_to_data:/data -t $DOCK \
    /bin/bash -c "STAR \
        --runMode genomeGenerate \
        --genomeDir /data/$star_index \
        --genomeFastaFiles /data/$genome_fasta \
        --sjdbGTFfile /data/gencode.v26.GRCh38.annotation.gtf \
        --sjdbOverhang $OH \
        --runThreadN 10"

# build the RSEM index:
mkdir rsem_reference
docker run --rm -v $path_to_data:/data -t $DOCK \
    /bin/bash -c "rsem-prepare-reference \
        /data/$genome_fasta \
        /data/rsem_reference/rsem_reference \
        --gtf /data/gencode.v26.GRCh38.annotation.gtf \
        --num-threads 10"

# STAR alignment
grep -v -f <(ls -lrt */*.tar.gz|awk -F"[ /]" '{print $(NF-1)}') sample_ReadLength.txt | while read LINE; do
#less sample_ReadLength.txt | tail -n+2 | while read LINE; do

  sample_id=`echo $LINE|awk '{print $1}'`
  OH=$(expr $(echo $LINE|awk '{print $2}') - 1)
  star_index=star_index_oh$OH
  #mkdir $sample_id

  #docker run --rm -v $path_to_data:/data -t $DOCK \
  #    /bin/bash -c "/src/run_STAR.py \
  #        /data/$star_index \
  #        /data/fastq/${sample_id}_1.fastq.gz \
  #        /data/fastq/${sample_id}_2.fastq.gz \
  #        ${sample_id} \
  #        --threads 10 \
  #        --output_dir /data/$sample_id" \
  #        > $sample_id/log.star.txt 2>&1
  #
  ### sync BAMs (optional; copy QC flags and read group IDs)
  ##docker run --rm -v $path_to_data:/data -t $DOCK \
  ##    /bin/bash -c "/src/run_bamsync.sh \
  ##        /data/$input_bam \
  ##        /data/star_out/${sample_id}.Aligned.sortedByCoord.out.bam \
  ##        /data/star_out/${sample_id}"
  #
  ## mark duplicates (Picard)
  #docker run --rm -v $path_to_data:/data -t $DOCK \
  #    /bin/bash -c "/src/run_MarkDuplicates.py \
  #        /data/$sample_id/$sample_id.Aligned.sortedByCoord.out.bam \
  #        $sample_id.Aligned.sortedByCoord.out.md \
  #        --output_dir /data/$sample_id/ --memory 20 " \
  #        > $sample_id/log.markdup.txt 2>&1
   
  # RNA-SeQC
  docker run --rm -v $path_to_data:/data -t $DOCK \
      /bin/bash -c "/src/run_rnaseqc.py \
      /data/$sample_id/$sample_id.Aligned.sortedByCoord.out.md.bam \
      /data/$genes_gtf \
      /data/$genome_fasta \
      $sample_id \
      --java_path /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java \
      --memory 32 \
      --output_dir /data/$sample_id"
  
  ## RSEM transcript quantification
  #docker run --rm -v $path_to_data:/data -t $DOCK \
  #    /bin/bash -c "/src/run_RSEM.py \
  #        /data/rsem_reference \
  #        /data/$sample_id/$sample_id.Aligned.toTranscriptome.out.bam \
  #        /data/$sample_id/$sample_id \
  #        --threads 6" \
  #        > $sample_id/log.rsem.txt 2>&1

done 


#Aggregating outputs
#Sample-level outputs in GCT format can be concatenated using combine_GCTs.py:
docker run --rm -v $path_to_data:/data -t $DOCK \
  /bin/bash -c "python3 /src/combine_GCTs.py \
    /data/gene_rpkm.list \
    /data/sample27.rnaseqc_rpkm"
    #/data/sample30.rnaseqc_reads"


# TMM norm
RPKM=sample27.rnaseqc_rpkm.gct.gz
READ=sample27.rnaseqc_reads.gct.gz
GTF=gencode.v26.GRCh38.genes.gtf
SAMPLELIST=sample_participant27.txt
VCFCHR=vcf_chr_list.txt
DOCK=mira:5000/broadinstitute/gtex_eqtl:latest

docker run --rm -v $path_to_data:/data -t $DOCK \
  /bin/bash -c "cp /data/eqtl_prepare_expression_mod.py /src/; /src/eqtl_prepare_expression_mod.py \
  /data/$RPKM /data/$READ /data/$GTF \
  /data/$SAMPLELIST /data/$VCFCHR /data/sample27.tmm-norm_count \
  --convert_tpm \
  --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 --normalization_method tmm"


ENSG00000173575.20
ENSG00000272888.5

