#!/usr/bin/env python3
"""Prepare a Salmon transcript FASTA and transcript-to-gene map from genome + GFF3/GTF.

Transcript sequence extraction is delegated to gffread, which correctly handles exon
order and strand. This helper then derives a two-column transcript->gene map from the
annotation and validates it against the transcript IDs written by gffread.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def parse_attrs(text: str) -> dict[str, str]:
    """Parse common GFF3 key=value and GTF key "value" attributes."""
    out: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key.strip()] = value.strip().strip('"')
        else:
            fields = item.split(None, 1)
            if len(fields) == 2:
                out[fields[0]] = fields[1].strip().strip('"')
    return out


def transcript_gene_map(annotation: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open_text(annotation) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2]
            attrs = parse_attrs(fields[8])
            if feature in {"mRNA", "transcript"}:
                tid = attrs.get("ID") or attrs.get("transcript_id")
                gid = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("gene")
                if tid:
                    mapping[tid] = (gid.split(",")[0] if gid else tid)
            # Some GTFs have no explicit transcript feature. Exon/CDS lines still carry IDs.
            elif feature in {"exon", "CDS"}:
                tid = attrs.get("transcript_id")
                gid = attrs.get("gene_id") or attrs.get("gene")
                if tid and tid not in mapping:
                    mapping[tid] = gid or tid
    return mapping


def fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def materialize_if_gz(src: Path, workdir: Path) -> tuple[Path, bool]:
    if src.suffix != ".gz":
        return src, False
    # gffread/samtools work most predictably with seekable, uncompressed FASTA/GFF files.
    name = src.name[:-3]
    dst = workdir / name
    with gzip.open(src, "rb") as inp, dst.open("wb") as out:
        shutil.copyfileobj(inp, out, length=8 * 1024 * 1024)
    return dst, True


def run(cmd: list[str]) -> None:
    print("+", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff3", required=True, help="GFF3 or GTF annotation; .gz is accepted")
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--tx2gene", required=True)
    args = ap.parse_args()

    genome = Path(args.genome).resolve()
    annotation = Path(args.gff3).resolve()
    transcripts = Path(args.transcripts).resolve()
    tx2gene = Path(args.tx2gene).resolve()
    for p, label in [(genome, "genome FASTA"), (annotation, "annotation")]:
        if not p.is_file() or p.stat().st_size == 0:
            raise SystemExit(f"ERROR: {label} not found or empty: {p}")

    transcripts.parent.mkdir(parents=True, exist_ok=True)
    tx2gene.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="salmon_ref_", dir=str(transcripts.parent)) as td:
        workdir = Path(td)
        genome_work, genome_was_gz = materialize_if_gz(genome, workdir)
        annotation_work, ann_was_gz = materialize_if_gz(annotation, workdir)
        if genome_was_gz:
            print(f"Temporarily decompressed genome: {genome_work}")
        if ann_was_gz:
            print(f"Temporarily decompressed annotation: {annotation_work}")

        # gffread can use a FASTA index and is much more memory-efficient than loading a
        # large plant genome into a Python dictionary.
        fai = Path(str(genome_work) + ".fai")
        if not fai.exists():
            run(["samtools", "faidx", str(genome_work)])

        run([
            "gffread", "-E", str(annotation_work),
            "-g", str(genome_work),
            "-w", str(transcripts),
        ])

    tids = fasta_ids(transcripts)
    if not tids:
        raise SystemExit("ERROR: gffread produced no transcript sequences")

    mapping = transcript_gene_map(annotation)
    missing = 0
    with tx2gene.open("w") as out:
        for tid in tids:
            gid = mapping.get(tid)
            if gid is None:
                # Safer fallback than guessing a gene by stripping version suffixes.
                gid = tid
                missing += 1
            out.write(f"{tid}\t{gid}\n")

    print(f"Wrote {len(tids)} transcript sequences: {transcripts}")
    print(f"Wrote transcript-to-gene map: {tx2gene}")
    if missing:
        print(
            f"WARNING: {missing} transcript IDs had no explicit gene mapping; "
            "they were mapped to themselves. Inspect the annotation if gene-level "
            "aggregation is important."
        )


if __name__ == "__main__":
    main()
