#!/usr/bin/env python3
"""PCA, optional t-SNE, and Pearson/Spearman correlation plots for an expression matrix."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

def load(path):
    d=pd.read_csv(path); idc=next((c for c in ['GeneID','gene_id','transcript_id'] if c in d.columns),None)
    if idc is None: raise SystemExit('ERROR: expected GeneID, gene_id, or transcript_id column')
    x=d.drop(columns=[idc]).apply(pd.to_numeric,errors='coerce').fillna(0.0).T
    return x

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--categories'); ap.add_argument('--cat_names'); ap.add_argument('--min',type=float); ap.add_argument('--top_var_pct',type=float)
    ap.add_argument('--tsne',action='store_true'); ap.add_argument('--pcs',default='1,2'); a=ap.parse_args()
    x=load(a.data)
    if a.min is not None: x=x.loc[:,x.mean(0)>=a.min]
    lx=np.log1p(x)
    if a.top_var_pct is not None:
        n=max(1,int(np.ceil(lx.shape[1]*a.top_var_pct/100))); keep=lx.var(0).nlargest(n).index; lx=lx[keep]
    if lx.shape[1]<2 or lx.shape[0]<2: raise SystemExit('ERROR: not enough data after filtering')
    z=StandardScaler().fit_transform(lx)
    pc=PCA(n_components=min(lx.shape[0],lx.shape[1])).fit_transform(z)
    p1,p2=[int(i)-1 for i in a.pcs.split(',')]
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    cats=[x.strip() for x in a.categories.split(',')] if a.categories else []
    labels=[]
    for sample in lx.index:
        hits=[c for c in cats if c in sample]
        labels.append(max(hits,key=len) if hits else 'Other')
    uniq=list(dict.fromkeys(labels))
    cmap=plt.get_cmap('tab20')
    cmap_map={u:cmap(i % 20) for i,u in enumerate(uniq)}
    plt.figure(figsize=(7,6))
    for u in uniq:
        idx=[i for i,v in enumerate(labels) if v==u]
        plt.scatter(pc[idx,p1],pc[idx,p2],label=u if cats else None)
    for i,s in enumerate(lx.index): plt.annotate(s,(pc[i,p1],pc[i,p2]),fontsize=7)
    if cats: plt.legend(fontsize=8)
    plt.xlabel(f'PC{p1+1}'); plt.ylabel(f'PC{p2+1}'); plt.tight_layout(); plt.savefig(a.output+'_PCA.jpg',dpi=300); plt.close()
    pear=lx.T.corr(method='pearson'); spea=lx.T.corr(method='spearman')
    for name,mat in [('Pearson',pear),('Spearman',spea)]:
        g=sns.clustermap(mat,cmap='vlag',vmin=-1,vmax=1,center=0,figsize=(8,8)); g.savefig(a.output+f'_{name}_correlation.jpg',dpi=300); plt.close(g.fig)
    if a.tsne:
        perplex=max(1,min(30,lx.shape[0]-1)); emb=TSNE(n_components=2,perplexity=perplex,init='pca',learning_rate='auto',random_state=1).fit_transform(z)
        plt.figure(figsize=(7,6)); plt.scatter(emb[:,0],emb[:,1]);
        for i,s in enumerate(lx.index): plt.annotate(s,(emb[i,0],emb[i,1]),fontsize=7)
        plt.tight_layout(); plt.savefig(a.output+'_TSNE.jpg',dpi=300); plt.close()
if __name__=='__main__': main()
