#!/usr/bin/env bash
###############################################################################
# mrna_lncrna_pipeline.sh
# Reproducible mRNA & lncRNA RNA-seq pipeline (fastp → STAR → featureCounts)
# Usage: ./mrna_lncrna_pipeline.sh SAMPLE1 [SAMPLE2 SAMPLE3 ...]
###############################################################################
set -Eeuo pipefail

############################### USER SETTINGS #################################
THREADS=8
GENOME_DIR="ref/index_GRCh38_STAR_mRNA"         # will be created if absent
GENOME_FA="ref/GRCh38.primary_assembly.genome.fa"
GTF_FILE="ref/gencode.v47.annotation.gtf"
RAW_DIR="raw"                                   # where SAMPLE_1/2.fq.gz live
WORK_DIR="work"                                 # results will be put here
###############################################################################

## ---------- 1. build STAR genome index once ---------------------------------
if [[ ! -d "$GENOME_DIR" ]]; then
  echo "• Building STAR index → $GENOME_DIR"
  mkdir -p "$GENOME_DIR"
  STAR --runThreadN "$THREADS" \
       --runMode genomeGenerate \
       --genomeDir "$GENOME_DIR" \
       --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$GTF_FILE" \
       --sjdbOverhang 149
  echo "  STAR index done ✔"
fi

## ---------- 2. iterate over all samples sent on the CLI ---------------------
for SAMPLE in "$@"; do
  echo -e "\n===================  Processing $SAMPLE  ==================="

  # create per-sample workspace
  OUTDIR="$WORK_DIR/$SAMPLE"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  R1="$PWD/../../$RAW_DIR/${SAMPLE}_1.fq.gz"
  R2="$PWD/../../$RAW_DIR/${SAMPLE}_2.fq.gz"

  ######################## fastp QC & trimming ################################
  echo "• fastp QC"
  fastp -w "$THREADS" \
        -i "$R1" -I "$R2" \
        -o "${SAMPLE}_1.clean.fq.gz" \
        -O "${SAMPLE}_2.clean.fq.gz" \
        --html  ${SAMPLE}.fastp.html \
        --json  ${SAMPLE}.fastp.json

  ######################## STAR alignment ####################################
  echo "• STAR alignment"
  STAR --runThreadN "$THREADS" \
       --genomeDir "$PWD/../../$GENOME_DIR" \
       --readFilesCommand zcat \
       --readFilesIn "${SAMPLE}_1.clean.fq.gz" "${SAMPLE}_2.clean.fq.gz" \
       --outFileNamePrefix "${SAMPLE}" \
       --outSAMtype BAM SortedByCoordinate \
       --outBAMsortingThreadN "$THREADS" \
       --quantMode TranscriptomeSAM GeneCounts

  ######################## featureCounts quantification #######################
  echo "• featureCounts quantification"
  featureCounts -T "$THREADS" \
                -p --countReadPairs \
                -t gene -f -g gene_id \
                -a "$PWD/../../$GTF_FILE" \
                -o "${SAMPLE}counts.txt" \
                   "${SAMPLE}Aligned.sortedByCoord.out.bam"

  echo "  ${SAMPLE} done ✔"
  cd - >/dev/null
done

echo -e "\n⚡ Pipeline completed for all samples!"
