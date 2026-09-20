#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.sh}"
source "$CONFIG_FILE"
for x in hisat2 hisat2-build samtools parallel; do command -v "$x" >/dev/null || { echo "ERROR: $x not found" >&2; exit 1; }; done
[[ -s "$HISAT2_INPUT_MANIFEST" ]] || { echo "ERROR: manifest not found: $HISAT2_INPUT_MANIFEST" >&2; exit 1; }
[[ -s "$HISAT2_REFERENCE_FASTA" ]] || { echo "ERROR: reference FASTA not found: $HISAT2_REFERENCE_FASTA" >&2; exit 1; }

if (( HISAT2_JOBS * (HISAT2_THREADS + 2 * SAMTOOLS_THREADS) > TOTAL_THREADS )); then
  echo "ERROR: HISAT2_JOBS * (HISAT2_THREADS + 2*SAMTOOLS_THREADS) exceeds TOTAL_THREADS ($TOTAL_THREADS)" >&2
  echo "       The two SAMtools stages are budgeted separately because pipeline stages may overlap." >&2
  exit 1
fi
OUT="$HISAT2_OUTPUT_DIR"; IDX="$OUT/index"; BAM="$OUT/bam"; LOG="$OUT/logs"; SCRIPTS="$OUT/sample_scripts"
PREFIX="$IDX/reference"
mkdir -p "$IDX" "$BAM" "$LOG" "$SCRIPTS"

if [[ ! -s "${PREFIX}.1.ht2" && ! -s "${PREFIX}.1.ht2l" ]]; then
  echo "Building HISAT2 index..."
  hisat2-build -p "$HISAT2_INDEX_THREADS" "$HISAT2_REFERENCE_FASTA" "$PREFIX"
else
  echo "Reusing existing HISAT2 index: $PREFIX"
fi

MODE_ARGS=""
if [[ "$HISAT2_REFERENCE_MODE" == "transcriptome" ]]; then MODE_ARGS="--no-spliced-alignment"; fi
rm -f "$SCRIPTS"/*.hisat2.sh

tail -n +2 "$HISAT2_INPUT_MANIFEST" | while IFS=$'\t' read -r sample layout r1 r2; do
  [[ -n "$sample" ]] || continue
  outbam="$BAM/${sample}.bam"; log="$LOG/${sample}.hisat2.log"; s="$SCRIPTS/${sample}.hisat2.sh"
  if [[ "$layout" == "PE" ]]; then READ_ARGS="-1 \"$r1\" -2 \"$r2\""; else READ_ARGS="-U \"$r1\""; fi
  cat > "$s" <<EOS
#!/usr/bin/env bash
set -euo pipefail
hisat2 -x "$PREFIX" $READ_ARGS -p "$HISAT2_THREADS" $MODE_ARGS $HISAT2_EXTRA_ARGS --no-unal 2> "$log" | \
  samtools view -@ "$SAMTOOLS_THREADS" -b -F 4 - | \
  samtools sort -@ "$SAMTOOLS_THREADS" -o "$outbam" -
samtools index -@ "$SAMTOOLS_THREADS" "$outbam"
EOS
  chmod +x "$s"
done

find "$SCRIPTS" -maxdepth 1 -type f -name '*.hisat2.sh' | sort | \
  parallel --jobs "$HISAT2_JOBS" --halt soon,fail=1 --joblog "$OUT/parallel_joblog.tsv" bash {}
echo "HISAT2 complete: $BAM"
