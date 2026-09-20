#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.sh}"
source "$CONFIG_FILE"
python "$SCRIPT_DIR/00_make_manifest.py" --fastq-dir "$FASTQ_DIR" --output "$RAW_MANIFEST"
bash "$SCRIPT_DIR/01_run_fastqc_multiqc.sh"
bash "$SCRIPT_DIR/02_run_trimmomatic.sh"
# Both quantification paths are run by default. Comment one out if you only want the other.
bash "$SCRIPT_DIR/03_run_hisat2.sh"
bash "$SCRIPT_DIR/04_run_salmon.sh"
