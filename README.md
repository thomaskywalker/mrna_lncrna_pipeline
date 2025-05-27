# mRNA & lncRNA STAR2-featureCounts Pipeline

## Overview
A reproducible Bash workflow that trims raw paired-end RNA‑seq reads with **fastp**, aligns them to the human GRCh38 reference using **STAR2**, and quantifies gene‑level counts with **featureCounts**. The script is designed to be idempotent—builds the STAR2 genome index only once and processes an arbitrary number of samples passed on the command line.

---

## Directory Structure
```text
project/
├── mrna_lncrna_pipeline.sh        # the main Bash script
├── raw/                           # raw FASTQ.gz files (R1/R2 for each sample)
│   └── SAMPLE_1.fq.gz
│   └── SAMPLE_2.fq.gz
├── ref/                           # reference genome and annotation
│   ├── GRCh38.primary_assembly.genome.fa
│   └── gencode.v47.annotation.gtf
└── work/                          # created automatically; per-sample outputs
```
*Feel free to rename folders—just update the corresponding variables at the top of the script.*

---

## Hardware Requirements
| Resource | Minimum | Recommended |
|----------|---------|-------------|
| Operating system | 64‑bit Linux or macOS | Same |
| CPU cores        | 4 threads     | ≥ 8 threads for faster runs |
| RAM              | 16 GB         | 32 GB for STAR2 index step |
| Disk space       | 50 GB free    | 100 GB + for large cohorts |

> **Tip 💡** STAR2 index creation is the most memory‑intensive step. Once the index is built, per‑sample alignment typically stays under ~8 GB RAM.

---

## Software Dependencies
All tools can be installed via **conda** (recommended) or any other package manager of your choice.

| Tool | Version tested | Conda install command |
|------|---------------|-----------------------|
| Bash | ≥ 4.0 | *(pre‑installed on most systems)* |
| fastp | ≥ 0.23.2 | `conda install -c bioconda fastp` |
| STAR2 | ≥ 2.7.10a | `conda install -c bioconda star2` |
| Subread / featureCounts | ≥ 2.0.3 | `conda install -c bioconda subread` |
| pigz *(optional)* | ≥ 2.4 | `conda install -c conda-forge pigz` |

### Create a dedicated environment
```bash
conda create -n rnaseq -c bioconda -c conda-forge fastp star2 subread pigz
conda activate rnaseq
```

---

## Reference Files
1. **Genome FASTA** – GRCh38 primary assembly (e.g. `GRCh38.primary_assembly.genome.fa`).  
2. **Annotation GTF** – GENCODE v47 comprehensive annotation (`gencode.v47.annotation.gtf`).

Place both files in the `ref/` directory (or adjust `GENOME_FA` / `GTF_FILE` in the script).

> STAR2 will generate its own index folder inside `ref/` the first time the script is run (`index_GRCh38_STAR_mRNA` by default).

---

## Quick Start
```bash
# 1. clone or copy this repository
cd project

# 2. make the pipeline executable (once)
chmod +x mrna_lncrna_pipeline.sh

# 3. drop your FASTQ pairs into raw/ using SAMPLE_1.fq.gz / SAMPLE_2.fq.gz naming

# 4. run the pipeline on one or multiple samples
./mrna_lncrna_pipeline.sh SAMPLE1 SAMPLE2 SAMPLE3
```
Per-sample progress is echoed to the terminal:
```text
• fastp QC
• STAR2 alignment
• featureCounts quantification
… done ✔
```

---

## Output Files (per sample)
| File | Description |
|------|-------------|
| `*_1/2.clean.fq.gz` | Adapter‑trimmed, quality‑filtered reads (fastp) |
| `*Aligned.sortedByCoord.out.bam` | Genome‑sorted BAM with alignments (STAR2) |
| `*Aligned.toTranscriptome.out.bam` | BAM in transcriptome coordinates (optional downstream quantification) |
| `*ReadsPerGene.out.tab` | STAR2‑generated raw counts table |
| `*counts.txt` & `*counts.txt.summary` | featureCounts gene‑level matrix & summary |
| `*.fastp.html` / `*.fastp.json` | QC reports |

All outputs live in `work/SAMPLE/` so you can delete or archive entire folders without affecting other samples.

---

## Customisation
* **Change read length** – adjust `sjdbOverhang` in the index‑building block (`readLength − 1`).
* **Turn on extra fastp filters** – add flags such as `--cut_front` or `--poly_g_min_len`.
* **Use UMI workflows** – insert UMI‑tools steps between fastp and STAR2.
* **Different organism** – substitute FASTA/GTF files and rebuild the STAR2 index.

---

## Troubleshooting
| Symptom | Possible fix |
|---------|--------------|
| `STAR2: command not found` | Verify STAR2 is on `$PATH` or inside the conda env. |
| STAR2 exits with out‑of‑memory error | Build index on a machine with ≥ 32 GB RAM, then copy the index folder back. |
| `featureCounts` zero counts | Check that the GTF matches the genome build and that `-t gene -g gene_id` suit your annotation. |
| Slow decompression | Install **pigz** and set `readFilesCommand pigz -dc` in the script. |

---

## Citations
If you use mrna_lncrna_pipeline.sh in your research, talks, or any published work, please cite it as:

**Huang, P.H. (2025). mRNA & lncRNA STAR2‑featureCounts pipeline. GitHub. https://github.com/thomaskywalker/mrna_lncrna_pipeline**

Also, please cite the original tools when publishing results: **fastp**, **STAR2**, **featureCounts**

---

## License
Distributed under the MIT License. See `LICENSE` file for details.
