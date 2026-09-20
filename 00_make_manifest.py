#!/usr/bin/env python3
"""Build a sample manifest that supports both paired-end and single-end FASTQ files."""
from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path

FASTQ_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.I)
PAIR_PATTERNS = [
    # sample_R1.fastq.gz / sample_R2.fastq.gz and sample_R1_001.fastq.gz
    re.compile(r"^(?P<sample>.+)_R(?P<mate>[12])(?P<chunk>_[0-9]+)?(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
    # sample_1.fastq.gz / sample_2.fastq.gz
    re.compile(r"^(?P<sample>.+)_(?P<mate>[12])(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
    # sample.R1.fastq.gz / sample.R2.fastq.gz
    re.compile(r"^(?P<sample>.+)\.R(?P<mate>[12])(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
    # sample.1.fastq.gz / sample.2.fastq.gz
    re.compile(r"^(?P<sample>.+)\.(?P<mate>[12])(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
]


def pair_key(path: Path):
    name = path.name
    for i, rx in enumerate(PAIR_PATTERNS):
        m = rx.match(name)
        if m:
            chunk = m.groupdict().get("chunk") or ""
            # chunk is part of the key so lane/chunk files are not silently merged.
            key = (i, m.group("sample"), chunk, m.group("ext").lower())
            return key, int(m.group("mate")), m.group("sample") + chunk
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--output", default="samples.tsv")
    ap.add_argument("--recursive", action="store_true")
    ap.add_argument("--orphan-policy", choices=["error", "single"], default="error",
                    help="How to treat a file named like R1/1 that has no mate (default: error).")
    args = ap.parse_args()

    root = Path(args.fastq_dir)
    if not root.is_dir():
        raise SystemExit(f"ERROR: FASTQ directory not found: {root}")
    it = root.rglob("*") if args.recursive else root.glob("*")
    files = sorted(p.resolve() for p in it if p.is_file() and FASTQ_RE.search(p.name))
    if not files:
        raise SystemExit(f"ERROR: no FASTQ files found in {root}")

    paired = {}
    unclassified = []
    for p in files:
        parsed = pair_key(p)
        if parsed is None:
            unclassified.append(p)
            continue
        key, mate, sample = parsed
        rec = paired.setdefault(key, {"sample": sample})
        if mate in rec:
            raise SystemExit(f"ERROR: duplicate mate {mate} for inferred sample {sample}: {p} and {rec[mate]}")
        rec[mate] = p

    rows = []
    consumed = set()
    for rec in paired.values():
        if 1 in rec and 2 in rec:
            rows.append((rec["sample"], "PE", str(rec[1]), str(rec[2])))
            consumed.update([rec[1], rec[2]])
        else:
            p = rec.get(1) or rec.get(2)
            if args.orphan_policy == "error":
                raise SystemExit(
                    f"ERROR: {p.name} looks like one mate of a pair but its mate is missing. "
                    "If this is truly single-end data, rerun with --orphan-policy single."
                )
            rows.append((rec["sample"], "SE", str(p), ""))
            consumed.add(p)

    for p in unclassified:
        sample = FASTQ_RE.sub("", p.name)
        rows.append((sample, "SE", str(p), ""))
        consumed.add(p)

    rows.sort(key=lambda x: x[0])
    names = [r[0] for r in rows]
    dup = sorted({x for x in names if names.count(x) > 1})
    if dup:
        raise SystemExit(
            "ERROR: duplicate inferred sample names: " + ", ".join(dup) +
            ". Rename files or edit the manifest manually so every sample name is unique."
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "layout", "read1", "read2"])
        w.writerows(rows)

    n_pe = sum(r[1] == "PE" for r in rows)
    n_se = sum(r[1] == "SE" for r in rows)
    print(f"Wrote {out}: {len(rows)} samples ({n_pe} PE, {n_se} SE)")


if __name__ == "__main__":
    main()
