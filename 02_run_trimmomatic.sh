#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.sh}"
source "$CONFIG_FILE"
for x in trimmomatic parallel python; do command -v "$x" >/dev/null || { echo "ERROR: $x not found" >&2; exit 1; }; done
[[ -s "$RAW_MANIFEST" ]] || { echo "ERROR: manifest not found: $RAW_MANIFEST" >&2; exit 1; }

OUT="$TRIMMOMATIC_OUTPUT_DIR"
PE_DIR="$OUT/paired"; SE_DIR="$OUT/single"; U_DIR="$OUT/unpaired"; LOG_DIR="$OUT/logs"; SCRIPT_OUT="$OUT/sample_scripts"
TRIM_MANIFEST="$OUT/trimmed_samples.tsv"
mkdir -p "$PE_DIR" "$SE_DIR" "$U_DIR" "$LOG_DIR" "$SCRIPT_OUT"

find_adapter() {
  local explicit="$1" name="$2"
  if [[ -n "$explicit" ]]; then [[ -s "$explicit" ]] && { printf '%s\n' "$explicit"; return; }; echo "ERROR: adapter not found: $explicit" >&2; return 1; fi
  local hit=""
  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    hit="$(find "$CONDA_PREFIX" -type f -path '*/adapters/*' -name "$name" -print -quit 2>/dev/null || true)"
  fi
  [[ -n "$hit" ]] || { echo "ERROR: cannot locate $name; set its explicit path in config.sh" >&2; return 1; }
  printf '%s\n' "$hit"
}
PE_ADAPTER="$(find_adapter "$PE_ADAPTER_FASTA" "$PE_ADAPTER_FILE_NAME")"
SE_ADAPTER="$(find_adapter "$SE_ADAPTER_FASTA" "$SE_ADAPTER_FILE_NAME")"

if (( TRIMMOMATIC_JOBS * TRIMMOMATIC_THREADS > TOTAL_THREADS )); then
  echo "ERROR: TRIMMOMATIC_JOBS * TRIMMOMATIC_THREADS exceeds TOTAL_THREADS ($TOTAL_THREADS)" >&2; exit 1
fi

rm -f "$SCRIPT_OUT"/*.sh
printf 'sample\tlayout\tread1\tread2\n' > "$TRIM_MANIFEST"
# Skip header. Preserve spaces in paths via tab-separated read.
tail -n +2 "$RAW_MANIFEST" | while IFS=$'\t' read -r sample layout r1 r2; do
  [[ -n "$sample" ]] || continue
  s="$SCRIPT_OUT/${sample}.trim.sh"
  log="$LOG_DIR/${sample}.trimmomatic.log"
  if [[ "$layout" == "PE" ]]; then
    o1="$PE_DIR/${sample}_R1.fastq.gz"; o2="$PE_DIR/${sample}_R2.fastq.gz"
    u1="$U_DIR/${sample}_R1.unpaired.fastq.gz"; u2="$U_DIR/${sample}_R2.unpaired.fastq.gz"
    printf '%s\tPE\t%s\t%s\n' "$sample" "$o1" "$o2" >> "$TRIM_MANIFEST"
    cat > "$s" <<EOS
#!/usr/bin/env bash
set -euo pipefail
JAVA_TOOL_OPTIONS="-Xmx$TRIMMOMATIC_JAVA_XMX" trimmomatic PE -threads $TRIMMOMATIC_THREADS -phred33 \
  "$r1" "$r2" "$o1" "$u1" "$o2" "$u2" \
  ILLUMINACLIP:"$PE_ADAPTER":$SEED_MISMATCHES:$PALINDROME_CLIP_THRESHOLD:$SIMPLE_CLIP_THRESHOLD:$MIN_ADAPTER_LENGTH:$KEEP_BOTH_READS \
  $TRIM_STEPS > "$log" 2>&1
EOS
  elif [[ "$layout" == "SE" ]]; then
    o1="$SE_DIR/${sample}.fastq.gz"
    printf '%s\tSE\t%s\t\n' "$sample" "$o1" >> "$TRIM_MANIFEST"
    cat > "$s" <<EOS
#!/usr/bin/env bash
set -euo pipefail
JAVA_TOOL_OPTIONS="-Xmx$TRIMMOMATIC_JAVA_XMX" trimmomatic SE -threads $TRIMMOMATIC_THREADS -phred33 \
  "$r1" "$o1" \
  ILLUMINACLIP:"$SE_ADAPTER":$SEED_MISMATCHES:$PALINDROME_CLIP_THRESHOLD:$SIMPLE_CLIP_THRESHOLD \
  $TRIM_STEPS > "$log" 2>&1
EOS
  else
    echo "ERROR: unknown layout '$layout' for $sample" >&2; exit 1
  fi
  chmod +x "$s"
done

find "$SCRIPT_OUT" -maxdepth 1 -type f -name '*.trim.sh' | sort | \
  parallel --jobs "$TRIMMOMATIC_JOBS" --halt soon,fail=1 --joblog "$OUT/parallel_joblog.tsv" bash {}

echo "Trimming complete. Downstream manifest: $TRIM_MANIFEST"
echo "PE orphan/unpaired reads are retained in $U_DIR but are not included in the downstream manifest."
