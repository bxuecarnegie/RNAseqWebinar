#!/usr/bin/env bash
# Central configuration for the non-SLURM RNAseqWebinar rewrite.
# Paths may be relative to the directory where you launch the scripts.

FASTQ_DIR="FASTQ"
RAW_MANIFEST="samples.tsv"
TRIMMED_MANIFEST="trimmomatic_output/trimmed_samples.tsv"

# General CPU budget. Override at launch, e.g. TOTAL_THREADS=32 bash 03_run_hisat2.sh
TOTAL_THREADS="${TOTAL_THREADS:-$(python - <<'PY'
import os
print(os.cpu_count() or 1)
PY
)}"

# 01 FastQC / MultiQC
FASTQC_OUTPUT_DIR="fastqc_multiqc_output"
FASTQC_JOBS="${FASTQC_JOBS:-4}"
FASTQC_THREADS="${FASTQC_THREADS:-1}"

# 02 Trimmomatic
TRIMMOMATIC_OUTPUT_DIR="trimmomatic_output"
TRIMMOMATIC_JOBS="${TRIMMOMATIC_JOBS:-2}"
TRIMMOMATIC_THREADS="${TRIMMOMATIC_THREADS:-4}"
TRIMMOMATIC_JAVA_XMX="${TRIMMOMATIC_JAVA_XMX:-4g}"
PE_ADAPTER_FILE_NAME="TruSeq3-PE.fa"
SE_ADAPTER_FILE_NAME="TruSeq3-SE.fa"
PE_ADAPTER_FASTA=""   # Optional explicit path. Blank = search inside the Conda env.
SE_ADAPTER_FASTA=""   # Optional explicit path. Blank = search inside the Conda env.
TRIM_STEPS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"
# ILLUMINACLIP settings
SEED_MISMATCHES=2
PALINDROME_CLIP_THRESHOLD=30
SIMPLE_CLIP_THRESHOLD=10
MIN_ADAPTER_LENGTH=2
KEEP_BOTH_READS=true

# 03 HISAT2
HISAT2_INPUT_MANIFEST="$TRIMMED_MANIFEST"  # Change to samples.tsv to skip trimming.
HISAT2_REFERENCE_FASTA="reference.fasta"
HISAT2_OUTPUT_DIR="hisat2_mapping"
HISAT2_JOBS="${HISAT2_JOBS:-2}"
HISAT2_THREADS="${HISAT2_THREADS:-8}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-2}"
HISAT2_INDEX_THREADS="${HISAT2_INDEX_THREADS:-8}"
# transcriptome/CDS: add --no-spliced-alignment. genome: leave normal spliced alignment enabled.
HISAT2_REFERENCE_MODE="transcriptome"  # transcriptome or genome
HISAT2_EXTRA_ARGS=""

# 04 Salmon
SALMON_INPUT_MANIFEST="$TRIMMED_MANIFEST" # Change to samples.tsv to skip trimming.
SALMON_OUTPUT_DIR="salmon_quant"
SALMON_JOBS="${SALMON_JOBS:-2}"
SALMON_THREADS="${SALMON_THREADS:-8}"
SALMON_INDEX_THREADS="${SALMON_INDEX_THREADS:-8}"
# Single-end Salmon cannot estimate fragment-length distribution from mate pairs.
SALMON_SE_FLD_MEAN="${SALMON_SE_FLD_MEAN:-250}"
SALMON_SE_FLD_SD="${SALMON_SE_FLD_SD:-25}"

# Option A: provide transcript FASTA + tx2gene mapping directly.
SALMON_TRANSCRIPTS_FASTA=""
SALMON_TX2GENE=""
# Option B: if Option A is blank, provide genome FASTA + GFF3 and the helper will build them.
GENOME_FASTA="genome.fasta"
GFF3_FILE="annotation.gff3"
