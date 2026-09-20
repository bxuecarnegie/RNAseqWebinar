#!/usr/bin/env bash
# DGX Spark preset for the non-SLURM workflow.
# Source the base config first so file/reference paths only need to be edited once.
SCRIPT_DIR_DGX="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR_DGX/config.sh"

# NVIDIA DGX Spark has a 20-core Arm CPU and 128 GB unified memory.
# These tools are CPU programs, so the preset targets CPU throughput rather than CUDA.
TOTAL_THREADS="${TOTAL_THREADS:-20}"
FASTQC_JOBS="${FASTQC_JOBS:-10}"
FASTQC_THREADS="${FASTQC_THREADS:-1}"
TRIMMOMATIC_JOBS="${TRIMMOMATIC_JOBS:-2}"
TRIMMOMATIC_THREADS="${TRIMMOMATIC_THREADS:-4}"
HISAT2_JOBS="${HISAT2_JOBS:-1}"
HISAT2_THREADS="${HISAT2_THREADS:-16}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-2}"
HISAT2_INDEX_THREADS="${HISAT2_INDEX_THREADS:-16}"
SALMON_JOBS="${SALMON_JOBS:-2}"
SALMON_THREADS="${SALMON_THREADS:-8}"
SALMON_INDEX_THREADS="${SALMON_INDEX_THREADS:-16}"
