#!/usr/bin/env bash
###############################################################################
# mrna_lncrna_pipeline.sh
# Reproducible mRNA & lncRNA RNA-seq pipeline
# fastp  →  STAR  →  featureCounts  (single combined run)
#
# Usage examples
#   ./mrna_lncrna_pipeline.sh  SAMPLE1 SAMPLE2 …
#   ./mrna_lncrna_pipeline.sh  all        # process every pair in raw/
#
# The script is idempotent:
#   • STAR genome index is built once and then reused.
#   • Per-sample steps create a dedicated sub-folder under work/.
#   • All BAM files are collected and counted together to produce
#     a wide gene × sample matrix (CSV) in work/.
###############################################################################
set -Eeuo pipefail        # E: fail in subshells, e: stop on error,
                          # u: undefined var is error, o pipefail: catch pipe errors

########################### USER-TUNABLE SETTINGS #############################
THREADS=8                                 # number of CPU threads for every tool
GENOME_DIR="ref/index_GRCh38_STAR_mRNA"   # STAR index will be created here
GENOME_FA="ref/GRCh38.primary_assembly.genome.fa"
GTF_FILE="ref/gencode.v48.annotation.gtf" # annotation (GTF)
RAW_DIR="raw"                             # FASTQ.gz files live here
WORK_DIR="work"                           # all output written here
###############################################################################

# -----------------------------------------------------------------------------#
# 1. Build STAR genome index (only if missing or incomplete)                   #
# -----------------------------------------------------------------------------#
# Check for the core file “genomeParameters.txt”; absence means index is
# missing or incomplete → regenerate.
if [[ ! -s "$GENOME_DIR/genomeParameters.txt" ]]; then
  echo "• Building STAR index → $GENOME_DIR"
  mkdir -p "$GENOME_DIR"
  STAR --runThreadN "$THREADS" \
       --runMode genomeGenerate \
       --genomeDir "$GENOME_DIR" \
       --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$GTF_FILE" \
       --sjdbOverhang 149          # read length minus 1 (150 bp reads)
  echo "  STAR index done ✔"
fi

###############################################################################
# 2. Handle the special keyword “all”                                          #
#    If user passes exactly one argument: all                                  #
#    → scan $RAW_DIR for every *_1.fq.gz file and build the sample list.       #
#    Otherwise use the arguments given on the command line.                    #
###############################################################################
if [[ "$#" -eq 1 && "$1" == "all" ]]; then
  echo "• 'all' flag detected – scanning $RAW_DIR for sample pairs"
  mapfile -t ARG_SAMPLES < <(
      find "$RAW_DIR" -maxdepth 1 -type f -name '*_1.f*q.gz' |
      sed -E 's#.*/##; s/_1\.f(ast)?q\.gz$//' | sort -u
  )
  echo "  Found ${#ARG_SAMPLES[@]} samples ➜ ${ARG_SAMPLES[*]}"
  [[ ${#ARG_SAMPLES[@]} -eq 0 ]] && {
      echo "❌  No *_1.fq.gz files in '$RAW_DIR' – aborting"; exit 1; }
else
  ARG_SAMPLES=("$@")
fi

###############################################################################
# 3. Per-sample loop: fastp trim  →  STAR alignment                            #
###############################################################################
declare -a BAM_LIST=()    # will collect path to each sorted BAM for counting

for SAMPLE in "${ARG_SAMPLES[@]}"; do
  echo -e "\n===================  Processing $SAMPLE  ==================="

  # Create per-sample working directory
  OUTDIR="$WORK_DIR/$SAMPLE"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  # Input FASTQ paths (paired-end)
  R1="$PWD/../../$RAW_DIR/${SAMPLE}_1.fq.gz"
  R2="$PWD/../../$RAW_DIR/${SAMPLE}_2.fq.gz"

  # ---------- fastp: QC + trimming ------------------------------------------
  echo "• fastp QC"
  fastp -w "$THREADS" \
        -i "$R1" -I "$R2" \
        -o "${SAMPLE}_1.clean.fq.gz" \
        -O "${SAMPLE}_2.clean.fq.gz" \
        --html  ${SAMPLE}.fastp.html \
        --json  ${SAMPLE}.fastp.json

  # ---------- STAR alignment ------------------------------------------------
  echo "• STAR alignment"
  # ensure any previous tmp directory is removed
  rm -rf "/tmp/${SAMPLE}_STARtmp"

  STAR --runThreadN "$THREADS" \
       --genomeDir "$PWD/../../$GENOME_DIR" \
       --readFilesCommand zcat \
       --readFilesIn "${SAMPLE}_1.clean.fq.gz" "${SAMPLE}_2.clean.fq.gz" \
       --outFileNamePrefix "${SAMPLE}" \
       --outSAMtype BAM SortedByCoordinate \
       --outBAMsortingThreadN "$THREADS" \
       --quantMode TranscriptomeSAM GeneCounts \
       --outTmpDir "/tmp/${SAMPLE}_STARtmp"

  # Collect BAM path for later featureCounts
  BAM_LIST+=("$OUTDIR/${SAMPLE}Aligned.sortedByCoord.out.bam")

  echo "  ${SAMPLE} done ✔"
  cd - >/dev/null
done

###############################################################################
# 4. Run featureCounts once on all BAMs to get a wide count matrix            #
###############################################################################
echo -e "\n• Running featureCounts on ${#BAM_LIST[@]} BAM files"

MERGED_TSV="$WORK_DIR/mRNA.lncRNA.counts.txt"
MERGED_CSV="$WORK_DIR/mRNA.lncRNA.counts.csv"

featureCounts -T "$THREADS" \
              -p --countReadPairs \
              -t gene -f -g gene_id \
              -a "$GTF_FILE" \
              -o "$MERGED_TSV" \
              "${BAM_LIST[@]}"

# Convert tab-separated output (skip 2 comment lines) to CSV
tail -n +2 "$MERGED_TSV" | tr '\t' ',' > "$MERGED_CSV"
echo "  ➜ Combined matrix saved to $MERGED_CSV"

echo -e "\n⚡ Pipeline completed for all samples!"
