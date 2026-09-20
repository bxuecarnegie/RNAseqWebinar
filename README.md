# RNAseqWebinar — non-SLURM + SE/PE rewrite

This folder is a non-SLURM rewrite of the workflow in `MattStata/RNAseqWebinar`.
It keeps the same overall sequence of operations:

1. FastQC + MultiQC
2. Trimmomatic
3. HISAT2 + SAMtools
4. Salmon
5. BAM/Salmon matrix generation and visualization

The main changes are:

- no `sbatch`, `#SBATCH`, environment modules, or SLURM variables;
- tools are supplied by Conda/Mamba;
- CPU use is controlled in `config.sh`;
- one TSV manifest describes both paired-end (PE) and single-end (SE) samples;
- Trimmomatic, HISAT2, and Salmon choose PE/SE commands per sample;
- mixed projects containing both PE and SE samples are allowed;
- a native `linux-aarch64` Conda environment and 20-core resource preset are supplied for NVIDIA DGX Spark;
- HISAT2 and Salmon indices are reused if already present.

## 1. Install Conda/Mamba

Miniforge/Mambaforge is recommended because the workflow uses `conda-forge` and `bioconda`.
If you already have Conda/Mamba, skip this step.

## 2. Create the normal environment

```bash
cd RNAseqWebinar_nonSLURM
mamba env create -f environment.yml
conda activate rnaseq-webinar
bash validate_environment.sh
```

If you only have `conda`, replace `mamba` with `conda`.

The environment pins Python 3.12 so an unrelated Python 3.13/3.14 pin in another environment does not force incompatible bioinformatics packages.

## 3. NVIDIA DGX Spark setup

DGX Spark is ARM64 Linux. Check:

```bash
uname -m
# expected: aarch64
```

If Conda is not installed, install the ARM64 Miniforge build:

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
bash Miniforge3-Linux-aarch64.sh
source ~/miniforge3/etc/profile.d/conda.sh
```

Then create the Spark environment:

```bash
cd RNAseqWebinar_nonSLURM
mamba env create -f environment_dgx_spark.yml
conda activate rnaseq-webinar-dgx-spark
bash validate_environment.sh
```

For the supplied 20-core preset, run any stage with `CONFIG_FILE=config_dgx_spark.sh`, for example:

```bash
CONFIG_FILE=config_dgx_spark.sh bash run_all.sh
```

`config_dgx_spark.sh` first loads your paths and references from `config.sh`, then overrides only concurrency/thread settings. No CUDA packages are installed because this workflow is CPU-based. The DGX Spark GPU is not used by FastQC, Trimmomatic, HISAT2, SAMtools, or Salmon.

## 4. Edit `config.sh`

At minimum set:

```bash
FASTQ_DIR="FASTQ"
HISAT2_REFERENCE_FASTA="reference.fasta"
HISAT2_REFERENCE_MODE="transcriptome"   # or genome
```

For Salmon use either:

### Option A: existing transcript FASTA + transcript-to-gene map

```bash
SALMON_TRANSCRIPTS_FASTA="transcripts.fa"
SALMON_TX2GENE="tx2gene.tsv"
```

where `tx2gene.tsv` contains two tab-delimited columns without a header:

```text
transcript_1    gene_1
transcript_2    gene_1
```

### Option B: genome FASTA + GFF3

```bash
SALMON_TRANSCRIPTS_FASTA=""
SALMON_TX2GENE=""
GENOME_FASTA="genome.fasta"
GFF3_FILE="annotation.gff3"
```

`prepare_salmon_reference.py` uses `gffread` to extract correctly spliced/strand-aware transcript sequences and then creates `tx2gene.tsv` from the annotation. It accepts a gzip-compressed genome/annotation by temporarily decompressing it in the Salmon reference directory; for a large genome, supplying uncompressed files avoids that extra I/O.

## 5. FASTQ naming and the sample manifest

Generate the manifest:

```bash
python 00_make_manifest.py \
  --fastq-dir FASTQ \
  --output samples.tsv
```

It recognizes common paired conventions:

```text
sample_R1.fastq.gz / sample_R2.fastq.gz
sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
sample_1.fastq.gz / sample_2.fastq.gz
sample.R1.fastq.gz / sample.R2.fastq.gz
sample.1.fastq.gz / sample.2.fastq.gz
```

A file without a recognized mate marker, for example:

```text
SRR123456.fastq.gz
```

is treated as single-end.

The output format is:

```text
sample  layout  read1  read2
SRR1    SE      /path/SRR1.fastq.gz
S1      PE      /path/S1_R1.fastq.gz  /path/S1_R2.fastq.gz
```

If a true SE file is unfortunately named like `sample_R1.fastq.gz`, the program assumes a mate may be missing and stops. Override that intentionally with:

```bash
python 00_make_manifest.py \
  --fastq-dir FASTQ \
  --output samples.tsv \
  --orphan-policy single
```

If your naming is unusual, you can simply create/edit `samples.tsv` manually.

## 6. Run FastQC + MultiQC

```bash
bash 01_run_fastqc_multiqc.sh
```

This QC step examines every FASTQ file regardless of SE/PE layout.

Main output:

```text
fastqc_multiqc_output/
  fastqc_reports/
  multiqc_report/multiqc_FASTQ_report.html
```

## 7. Run Trimmomatic

```bash
bash 02_run_trimmomatic.sh
```

PE samples use Trimmomatic `PE`; SE samples use Trimmomatic `SE`.
The script creates:

```text
trimmomatic_output/
  paired/                 # surviving PE mates
  single/                 # trimmed SE samples
  unpaired/               # orphan reads produced while trimming PE data
  logs/
  trimmed_samples.tsv     # downstream manifest
```

The PE orphan files are retained, but intentionally excluded from `trimmed_samples.tsv`. This keeps the downstream HISAT2/Salmon library definition internally consistent.

If you want to skip trimming, set these in `config.sh`:

```bash
HISAT2_INPUT_MANIFEST="samples.tsv"
SALMON_INPUT_MANIFEST="samples.tsv"
```

## 8. Run HISAT2

```bash
bash 03_run_hisat2.sh
```

For PE samples the script runs HISAT2 with `-1/-2`; for SE samples it uses `-U`.
Sorted/indexed BAM files are written to:

```text
hisat2_mapping/bam/
```

### Transcript/CDS reference

If `HISAT2_REFERENCE_FASTA` contains transcript or CDS sequences:

```bash
HISAT2_REFERENCE_MODE="transcriptome"
```

The script adds `--no-spliced-alignment`, matching the intent of the original webinar workflow.

### Genomic reference

If the reference is a genome assembly:

```bash
HISAT2_REFERENCE_MODE="genome"
```

This leaves spliced alignment enabled.

**Important:** `QuantifyAllBAM.py` below counts alignments by FASTA record ID. That is suitable when your HISAT2 reference records are genes/transcripts/CDSs. If you align against chromosomes/contigs in a genome, use an annotation-aware counter such as featureCounts/HTSeq instead of treating chromosome IDs as genes.

## 9. Run Salmon

```bash
bash 04_run_salmon.sh
```

PE samples use:

```text
-1 read1 -2 read2
```

SE samples use:

```text
-r reads --fldMean <mean> --fldSD <sd>
```

Unlike PE data, Salmon cannot infer the fragment-length distribution from paired mappings for SE libraries. The defaults in `config.sh` are `SALMON_SE_FLD_MEAN=250` and `SALMON_SE_FLD_SD=25`; if your library preparation provides better fragment-size estimates, change these values because they affect effective lengths and TPM estimates.

Output:

```text
salmon_quant/quants/<sample>/quant.sf
salmon_quant/quants/<sample>/quant.genes.sf
```

## 10. Generate matrices from HISAT2 BAMs

For a transcript/CDS-reference HISAT2 run:

```bash
python QuantifyAllBAM.py \
  --fasta reference.fasta \
  --bam_dir hisat2_mapping/bam \
  --manifest trimmomatic_output/trimmed_samples.tsv \
  --output hisat2_expression \
  --threads 8 \
  --progress
```

Outputs:

```text
hisat2_expression.counts.csv
hisat2_expression.tpm.csv
```

The manifest lets the program count PE fragments once while counting SE alignments individually. For a uniform SE-only BAM directory, `--single` remains available as a compatibility shortcut.

## 11. Generate matrices from Salmon

```bash
python QuantifySalmon.py \
  --quant-root salmon_quant/quants \
  --metric Both \
  --which both \
  --outdir salmon_matrices
```

Typical outputs include:

```text
*.transcripts.TPM.csv
*.transcripts.NumReads.csv
*.genes.TPM.csv
*.genes.NumReads.csv
```

## 12. Visualization

Example:

```bash
python VisualizeResults.py \
  --data salmon_matrices/quants.genes.TPM.csv \
  --output plots/salmon_genes \
  --top_var_pct 10 \
  --tsne
```

The rewritten visualization script produces PCA, Pearson correlation, Spearman correlation, and optional t-SNE images.

## 13. Run the main preprocessing/quantification workflow in one command

After editing `config.sh`:

```bash
bash run_all.sh
```

This runs manifest generation, QC, trimming, HISAT2, and Salmon. Comment out the HISAT2 or Salmon line in `run_all.sh` if you only need one quantification strategy.

## 14. Resource tuning on a workstation

The original SLURM scripts requested large fixed allocations. On a workstation, concurrency must instead be chosen based on actual CPU and RAM.

`config.sh` separates **jobs** from **threads per job**. For example:

```bash
TOTAL_THREADS=32
HISAT2_JOBS=2
HISAT2_THREADS=12
SAMTOOLS_THREADS=2
```

Because the HISAT2 → `samtools view` → `samtools sort` pipeline can overlap, budget roughly `HISAT2_THREADS + 2 × SAMTOOLS_THREADS` per simultaneous sample. With `HISAT2_JOBS=2`, `HISAT2_THREADS=12`, and `SAMTOOLS_THREADS=2`, the conservative budget is `2 × (12 + 2×2) = 32` threads.

For a smaller laptop, use fewer simultaneous jobs, for example:

```bash
FASTQC_JOBS=2
TRIMMOMATIC_JOBS=1
TRIMMOMATIC_THREADS=4
HISAT2_JOBS=1
HISAT2_THREADS=6
SAMTOOLS_THREADS=2
SALMON_JOBS=1
SALMON_THREADS=6
```

For DGX Spark, you can raise concurrency, but memory consumption still depends on the reference/index and read depth. Increase job counts gradually rather than starting with all cores occupied.

## 15. Notes about the original Python scripts

The original repository's `QuantifySalmon.py` and `VisualizeResults.py` are not intrinsically tied to SLURM or FASTQ layout; they operate on generated quantification matrices/results. The original `QuantifyAllBAM.py` also already has a `--single` option. This package includes rewritten versions so the folder is self-contained and so mixed SE/PE BAM directories can use a manifest explicitly.

## 16. Lane/chunk behavior

The automatic manifest builder deliberately does **not** merge multiple lane/chunk files for the same biological sample. For example, `sample_R1_001.fastq.gz`/`sample_R2_001.fastq.gz` and `_002` are treated as separate entries. This prevents silent merging. If those files are sequencing lanes for one sample, either concatenate corresponding R1 and R2 lanes before manifest generation or create a manifest appropriate for your desired lane handling.

## 17. Using an alternate config

Every shell stage accepts the `CONFIG_FILE` environment variable. This is useful for machine-specific resource presets without duplicating the workflow scripts:

```bash
CONFIG_FILE=config_dgx_spark.sh bash 03_run_hisat2.sh
```

If `CONFIG_FILE` is omitted, the scripts use `config.sh`.
