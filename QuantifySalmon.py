#!/usr/bin/env python3
"""Convert per-sample Salmon quant.sf / quant.genes.sf files into wide CSV matrices."""
from __future__ import annotations
import argparse, csv
from pathlib import Path

def read_sf(p, metric):
    with p.open(newline='') as fh:
        r=csv.DictReader(fh,delimiter='\t'); return {x['Name']:x[metric] for x in r}

def build(root,samples,filename,id_header,metric,out,missing):
    data={}; ids=set()
    for s in samples:
        p=root/s/filename
        if not p.exists(): raise SystemExit(f'ERROR: missing {p}')
        data[s]=read_sf(p,metric); ids.update(data[s])
    with out.open('w',newline='') as fh:
        w=csv.writer(fh); w.writerow([id_header]+samples)
        for i in sorted(ids): w.writerow([i]+[data[s].get(i,missing) for s in samples])

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--quant-root',required=True); ap.add_argument('--samples',nargs='+')
    ap.add_argument('--metric',choices=['Both','TPM','NumReads'],default='Both'); ap.add_argument('--outdir')
    ap.add_argument('--prefix'); ap.add_argument('--which',choices=['both','tx','gene'],default='both'); ap.add_argument('--missing-value',default='0')
    a=ap.parse_args(); root=Path(a.quant_root); samples=a.samples or sorted(p.name for p in root.iterdir() if p.is_dir() and (p/'quant.sf').exists())
    out=Path(a.outdir) if a.outdir else root.parent/'salmon_matrices'; out.mkdir(parents=True,exist_ok=True)
    prefix=a.prefix or root.name.removesuffix('_salmon_quant').removesuffix('_quant')
    metrics=['TPM','NumReads'] if a.metric=='Both' else [a.metric]
    for m in metrics:
        if a.which in ('both','tx'): build(root,samples,'quant.sf','transcript_id',m,out/f'{prefix}.transcripts.{m}.csv',a.missing_value)
        if a.which in ('both','gene'):
            # Salmon --geneMap normally creates quant.genes.sf.
            build(root,samples,'quant.genes.sf','gene_id',m,out/f'{prefix}.genes.{m}.csv',a.missing_value)
    print(f'Wrote matrices to {out}')
if __name__=='__main__': main()
