#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.sh}"
source "$CONFIG_FILE"
for x in salmon parallel python; do command -v "$x" >/dev/null || { echo "ERROR: $x not found" >&2; exit 1; }; done
[[ -s "$SALMON_INPUT_MANIFEST" ]] || { echo "ERROR: manifest not found: $SALMON_INPUT_MANIFEST" >&2; exit 1; }
if (( SALMON_JOBS * SALMON_THREADS > TOTAL_THREADS )); then echo "ERROR: SALMON_JOBS * SALMON_THREADS exceeds TOTAL_THREADS ($TOTAL_THREADS)" >&2; exit 1; fi

OUT="$SALMON_OUTPUT_DIR"; REF="$OUT/reference"; IDX="$OUT/index"; Q="$OUT/quants"; LOG="$OUT/logs"; SCRIPTS="$OUT/sample_scripts"
mkdir -p "$REF" "$IDX" "$Q" "$LOG" "$SCRIPTS"
TX="$SALMON_TRANSCRIPTS_FASTA"; MAP="$SALMON_TX2GENE"
if [[ -z "$TX" || -z "$MAP" ]]; then
  for x in gffread samtools; do command -v "$x" >/dev/null || { echo "ERROR: $x not found (needed for genome+annotation Salmon reference preparation)" >&2; exit 1; }; done
  [[ -s "$GENOME_FASTA" ]] || { echo "ERROR: GENOME_FASTA missing: $GENOME_FASTA" >&2; exit 1; }
  [[ -s "$GFF3_FILE" ]] || { echo "ERROR: GFF3_FILE missing: $GFF3_FILE" >&2; exit 1; }
  TX="$REF/transcripts.fa"; MAP="$REF/tx2gene.tsv"
  if [[ ! -s "$TX" || ! -s "$MAP" ]]; then
    python "$SCRIPT_DIR/prepare_salmon_reference.py" --genome "$GENOME_FASTA" --gff3 "$GFF3_FILE" --transcripts "$TX" --tx2gene "$MAP"
  fi
else
  [[ -s "$TX" && -s "$MAP" ]] || { echo "ERROR: supplied Salmon transcript FASTA or tx2gene file is missing" >&2; exit 1; }
fi

if [[ ! -s "$IDX/versionInfo.json" ]]; then
  salmon index -t "$TX" -i "$IDX" -p "$SALMON_INDEX_THREADS"
else
  echo "Reusing Salmon index: $IDX"
fi

rm -f "$SCRIPTS"/*.salmon.sh
tail -n +2 "$SALMON_INPUT_MANIFEST" | while IFS=$'\t' read -r sample layout r1 r2; do
  [[ -n "$sample" ]] || continue
  s="$SCRIPTS/${sample}.salmon.sh"; q="$Q/$sample"; log="$LOG/${sample}.salmon.log"
  if [[ "$layout" == "PE" ]]; then
    READ_ARGS="-1 \"$r1\" -2 \"$r2\""
    LAYOUT_ARGS=""
  elif [[ "$layout" == "SE" ]]; then
    READ_ARGS="-r \"$r1\""
    LAYOUT_ARGS="--fldMean $SALMON_SE_FLD_MEAN --fldSD $SALMON_SE_FLD_SD"
  else
    echo "ERROR: unknown layout '$layout' for $sample" >&2
    exit 1
  fi
  cat > "$s" <<EOS
#!/usr/bin/env bash
set -euo pipefail
salmon quant -i "$IDX" -l A $READ_ARGS $LAYOUT_ARGS -p "$SALMON_THREADS" \
  --validateMappings --seqBias --gcBias --geneMap "$MAP" -o "$q" > "$log" 2>&1
EOS
  chmod +x "$s"
done
find "$SCRIPTS" -maxdepth 1 -type f -name '*.salmon.sh' | sort | \
  parallel --jobs "$SALMON_JOBS" --halt soon,fail=1 --joblog "$OUT/parallel_joblog.tsv" bash {}

missing=0
while IFS=$'\t' read -r sample layout r1 r2; do
  [[ "$sample" == "sample" || -z "$sample" ]] && continue
  [[ -s "$Q/$sample/quant.sf" ]] || { echo "ERROR: missing $Q/$sample/quant.sf" >&2; missing=1; }
done < "$SALMON_INPUT_MANIFEST"
(( missing == 0 )) || exit 1
echo "Salmon complete: $Q"
