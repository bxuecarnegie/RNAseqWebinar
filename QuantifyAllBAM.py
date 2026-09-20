#!/usr/bin/env python3
"""Build fragment-count and TPM matrices from HISAT2 BAMs, with per-sample SE/PE handling."""
from __future__ import annotations
import argparse, csv, subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def fasta_lengths(path):
    out, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                cur=line[1:].split()[0]; out[cur]=0
            elif cur is not None: out[cur]+=len(line)
    return out


def layouts_from_manifest(path):
    if not path: return {}
    with open(path, newline='') as fh:
        return {r['sample']: r['layout'].upper() for r in csv.DictReader(fh, delimiter='\t')}


def count_bam(args):
    bam, layout, valid = args; sample=bam.stem; counts=defaultdict(int)
    p=subprocess.Popen(['samtools','view',str(bam)], stdout=subprocess.PIPE, text=True)
    assert p.stdout is not None
    for line in p.stdout:
        f=line.rstrip('\n').split('\t')
        if len(f)<7: continue
        flag=int(f[1]); rname=f[2]; rnext=f[6]
        if rname=='*' or rname not in valid or flag & 0x100 or flag & 0x800: continue
        if layout=='PE':
            if (flag & 0x2) and (flag & 0x40) and rnext=='=': counts[rname]+=1
        else:
            if not (flag & 0x4): counts[rname]+=1
    rc=p.wait()
    if rc: raise RuntimeError(f'samtools view failed for {bam}')
    return sample, dict(counts)


def write_matrix(path, ids, samples, data, value_fn):
    with open(path,'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['GeneID']+samples)
        for gid in ids: w.writerow([gid]+[value_fn(s,gid) for s in samples])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--fasta',required=True); ap.add_argument('--bam_dir',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--manifest',help='TSV with sample/layout/read1/read2; recommended for mixed SE/PE projects')
    ap.add_argument('--layout',choices=['PE','SE'],default='PE',help='Fallback layout for samples absent from manifest')
    ap.add_argument('--single',action='store_true',help='Backward-compatible alias for --layout SE')
    ap.add_argument('--threads',type=int,default=4); ap.add_argument('--decimals',type=int,default=4); ap.add_argument('--progress',action='store_true')
    a=ap.parse_args(); fallback='SE' if a.single else a.layout
    lengths=fasta_lengths(a.fasta); ids=list(lengths); valid=set(ids); lm=layouts_from_manifest(a.manifest)
    bams=sorted(Path(a.bam_dir).glob('*.bam'))
    if not bams: raise SystemExit(f'ERROR: no BAM files in {a.bam_dir}')
    jobs=[(b, lm.get(b.stem,fallback), valid) for b in bams]
    results={}
    with ProcessPoolExecutor(max_workers=a.threads) as ex:
        futs={ex.submit(count_bam,j):j[0].stem for j in jobs}
        for i,f in enumerate(as_completed(futs),1):
            s,c=f.result(); results[s]=c
            if a.progress: print(f'[{i}/{len(futs)}] {s}')
    samples=sorted(results)
    write_matrix(a.output+'.counts.csv',ids,samples,results,lambda s,g:results[s].get(g,0))
    # Same normalization convention as the original script: count/length then scale to 1e6.
    rpk={s:{g:(results[s].get(g,0)/lengths[g] if lengths[g] else 0.0) for g in ids} for s in samples}
    denom={s:sum(rpk[s].values())/1e6 for s in samples}
    write_matrix(a.output+'.tpm.csv',ids,samples,results,
                 lambda s,g:round(rpk[s][g]/denom[s],a.decimals) if denom[s] else 0.0)
    print(a.output+'.counts.csv'); print(a.output+'.tpm.csv')
if __name__=='__main__': main()
