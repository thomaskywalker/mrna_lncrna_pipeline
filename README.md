
# mRNA & lncRNA STAR-featureCounts Pipeline

## Overview
A fully bash-based workflow that  
1. **trims** and QC-checks paired-end reads with **fastp**;  
2. **aligns** them to the human GRCh38 reference genome with **STAR**;  
3. **counts** fragments per gene with **featureCounts** (single combined run);  
4. Exports a ready‑to‑use **CSV count matrix** for all samples.

The pipeline is **idempotent**: the STAR genome index is built only once, and reruns automatically skip that step if the index already exists.

---

## Directory Layout
```text
project/
├── mrna_lncrna_pipeline.sh        # this script (make it executable)
├── raw/                           # raw FASTQ(.gz) files
│   ├── SAMPLE_A_1.fq.gz           # pair naming: _1 / _2
│   └── SAMPLE_A_2.fq.gz
├── ref/                           # reference genome & annotation
│   ├── GRCh38.primary_assembly.genome.fa
│   └── gencode.v48.annotation.gtf
└── work/                          # created automatically; all outputs
```
*Keep the above structure or adjust the variables at the top of the script.*

---

## Hardware
| Resource | Minimum | Recommended |
|----------|---------|-------------|
| OS       | 64‑bit Linux / WSL2 / macOS | Same |
| CPU      | 4 threads | ≥ 8 threads |
| RAM      | 16 GB | 32 GB for smoother STAR indexing |
| Disk     | 60 GB | 150 GB+ for multi‑sample projects |

---

## Software (conda / mamba)
| Tool | Version | Install command |
|------|---------|-----------------|
| fastp | ≥ 0.24 | `mamba install -c bioconda fastp` |
| STAR  | ≥ 2.7  | `mamba install -c bioconda star` |
| Subread / featureCounts | ≥ 2.0 | `mamba install -c bioconda subread` |
| pigz *(optional)* | ≥ 2.4 | `mamba install -c conda-forge pigz` |

### One‑shot environment creation
```bash
# install mamba once (if not present)
conda install -n base -c conda-forge mamba

# create RNA‑seq environment
mamba create -n rnaseq -c bioconda -c conda-forge       fastp star subread pigz samtools
conda activate rnaseq
```
*(The script can auto‑activate this env if you forget—see header.)*

---

## Reference Files
Download from GENCODE (v48) and place in `ref/`:

* `GRCh38.primary_assembly.genome.fa.gz` → **gunzip** → `.fa`  
* `gencode.v48.annotation.gtf.gz` → **gunzip** → `.gtf`

The first run will create `ref/index_GRCh38_STAR_mRNA/` (~25 GB).

---

## Running the Pipeline

### 1. Make executable 
```bash
# one-time (if edited on Windows): convert CRLF to LF
dos2unix mrna_lncrna_pipeline.sh
# VScode users: click the "CRLF" indicator at the bottom-right corner → choose "LF" → save

chmod +x mrna_lncrna_pipeline.sh   # make the script executable
```
### 2. Examples
```bash
# run specific samples
./mrna_lncrna_pipeline.sh SAMPLE_A SAMPLE_B

# OR run all samples (defined by *_1.fq.gz) under raw/
./mrna_lncrna_pipeline.sh all
```
Progress log looks like:
```text
• fastp QC
• STAR alignment
• Running featureCounts on 12 BAM files
➜ Combined matrix saved to work/mRNA.lncRNA.counts.csv
```

---

## Outputs
| Path | Description |
|------|-------------|
| `work/SAMPLE/…clean.fq.gz` | trimmed reads (fastp) |
| `work/SAMPLE/*Aligned.sortedByCoord.out.bam` | coordinate‑sorted genome BAM |
| `work/SAMPLE/*ReadsPerGene.out.tab` | STAR per‑sample raw counts |
| **`work/mRNA.lncRNA.counts.csv`** | gene × sample count matrix |

---

## Customisation Tips
* **Other read length** → change `sjdbOverhang` (readLen − 1) in the index block.  
* **Different filename pattern** (`_R1.fastq.gz`) → adjust both the `find` pattern and `R1/R2` variable lines.  
* **Disk‑tight?** → drop `TranscriptomeSAM`, delete `*.clean.fq.gz` at the end of each loop, or remove the intermediate TSV after CSV is generated.  

---

## How to Cite
If you use **mrna_lncrna_pipeline.sh** in your research, talks or any published works, please cite:

> Huang, P. H. (2024). mRNA-lncRNA Pipeline (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.15542604

Please additionally cite the original tools:

* **fastp** 
* **STAR**
* **featureCounts**  

---

## License
MIT — see `LICENSE`.
