#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.sh}"
source "$CONFIG_FILE"

for x in fastqc multiqc parallel; do command -v "$x" >/dev/null || { echo "ERROR: $x not found in PATH" >&2; exit 1; }; done
[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: FASTQ_DIR not found: $FASTQ_DIR" >&2; exit 1; }

OUT="$FASTQC_OUTPUT_DIR"
FQ_OUT="$OUT/fastqc_reports"
MQ_OUT="$OUT/multiqc_report"
LIST="$OUT/fastq_files.txt"
mkdir -p "$FQ_OUT" "$MQ_OUT"
find "$FASTQ_DIR" -maxdepth 1 -type f \( -iname '*.fastq' -o -iname '*.fq' -o -iname '*.fastq.gz' -o -iname '*.fq.gz' \) | sort > "$LIST"
[[ -s "$LIST" ]] || { echo "ERROR: no FASTQ files found in $FASTQ_DIR" >&2; exit 1; }

if (( FASTQC_JOBS * FASTQC_THREADS > TOTAL_THREADS )); then
  echo "ERROR: FASTQC_JOBS * FASTQC_THREADS exceeds TOTAL_THREADS ($TOTAL_THREADS)" >&2; exit 1
fi

parallel --jobs "$FASTQC_JOBS" --halt soon,fail=1 --joblog "$OUT/fastqc_parallel_joblog.tsv" \
  fastqc --quiet --threads "$FASTQC_THREADS" --outdir "$FQ_OUT" :::: "$LIST"
multiqc "$FQ_OUT" --outdir "$MQ_OUT" --filename multiqc_FASTQ_report.html --force

echo "FastQC/MultiQC complete: $MQ_OUT/multiqc_FASTQ_report.html"
