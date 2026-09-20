#!/usr/bin/env bash
set -euo pipefail
for x in python fastqc multiqc trimmomatic hisat2 hisat2-build samtools gffread salmon parallel java; do
  printf '%-14s ' "$x"
  command -v "$x" >/dev/null && echo "OK ($(command -v "$x"))" || { echo "MISSING"; exit 1; }
done
python - <<'PY'
import numpy,pandas,scipy,sklearn,matplotlib,seaborn,rich
print('Python scientific stack: OK')
PY
hisat2 --version | head -n 1
samtools --version | head -n 1
salmon --version
fastqc --version
multiqc --version
