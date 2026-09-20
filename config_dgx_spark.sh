#!/usr/bin/env bash
# NVIDIA DGX Spark resource preset (20 ARM CPU cores).
# Source the main config first so paths/reference settings remain centralized.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Total CPU budget available to the workflow.
TOTAL_THREADS=20

# QC / trimming
FASTQC_JOBS=4
TRIMMOMATIC_JOBS=2
TRIMMOMATIC_THREADS=8

# HISAT2 -> samtools view -> samtools sort
# Conservative budget per concurrent sample:
#   HISAT2_THREADS + 2 * SAMTOOLS_THREADS
# Here: 6 + 2*2 = 10 threads/sample; 2 samples = 20 threads total.
HISAT2_JOBS=2
HISAT2_THREADS=6
SAMTOOLS_THREADS=2

# Salmon: two samples concurrently, 10 threads each.
SALMON_JOBS=2
SALMON_THREADS=10

export TOTAL_THREADS
export FASTQC_JOBS
export TRIMMOMATIC_JOBS TRIMMOMATIC_THREADS
export HISAT2_JOBS HISAT2_THREADS SAMTOOLS_THREADS
export SALMON_JOBS SALMON_THREADS
