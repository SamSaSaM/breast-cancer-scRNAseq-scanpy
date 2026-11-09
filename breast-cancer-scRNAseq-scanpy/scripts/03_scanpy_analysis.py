#!/usr/bin/env python3
"""Pilot Scanpy QC, PCA, UMAP, Leiden and marker plots for Wu et al. 2021."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from scipy import sparse

DEFAULT_DATA = Path.home()/'Downloads'/'Wu_etal_2021_BRCA_scRNASeq'

def save_figure(path):
    plt.savefig(path,dpi=180,bbox_inches='tight')
    plt.close('all')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--data-dir',type=Path,default=DEFAULT_DATA)
    p.add_argument('--output-dir',type=Path,default=None)
    p.add_argument('--min-genes',type=int,default=200)
    p.add_argument('--max-mito',type=float,default=20.0)
    p.add_argument('--hvg',type=int,default=2000)
    args=p.parse_args()
    out=args.output_dir.expanduser() if args.output_dir else args.data_dir.expanduser()/'analysis'
    fig=out/'figures'
    fig.mkdir(parents=True,exist_ok=True)
    adata=sc.read_h5ad(out/'wu2021_pilot_raw.h5ad')
    print('Loaded:',adata.shape,flush=True)
    if not sparse.issparse(adata.X):
        raise ValueError('Expected sparse input matrix')
    # Count-based QC; use author-provided percent.mito because gene identifiers may vary.
    sc.pp.calculate_qc_metrics(adata,percent_top=None,inplace=True)
    qc_cols=['total_counts','n_genes_by_counts','percent.mito']
    adata.obs[qc_cols].describe().to_csv(out/'qc_summary_before_filter.csv')
    sc.pl.violin(adata,qc_cols,multi_panel=True,jitter=0.1,show=False)
    save_figure(fig/'01_qc_violin.png')
    # These are pilot thresholds, not a reproduction of the original paper's QC.
    keep=(adata.obs['n_genes_by_counts']>=args.min_genes)&(adata.obs['percent.mito']<args.max_mito)
    adata=adata[keep].copy()
    if adata.n_obs<100 or adata.n_vars<100:
        raise ValueError('Too few cells or genes after filtering')
    sc.pp.filter_genes(adata,min_cells=3)
    # Keep a sparse copy of integer counts before log-normalization.
    adata.layers['counts']=adata.X.copy()
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=min(args.hvg,adata.n_vars-1),flavor='seurat')
    # Store log-normalized all-gene values for marker plots; remains sparse.
    adata.raw=adata.copy()
    adata=adata[:,adata.var['highly_variable']].copy()
    n_pcs=min(40,adata.n_obs-1,adata.n_vars-1)
    if n_pcs<3:
        raise ValueError('Not enough cells or genes for PCA')
    sc.pp.pca(adata,n_comps=n_pcs,zero_center=False,random_state=42)
    sc.pp.neighbors(adata,n_neighbors=min(15,adata.n_obs-1),n_pcs=min(30,n_pcs),random_state=42)
    sc.tl.umap(adata,random_state=42)
    sc.tl.leiden(adata,resolution=0.5,key_added='leiden',flavor='igraph',directed=False)
    for col,name in [('celltype_major','02_umap_author_celltype.png'),
                     ('subtype','03_umap_subtype.png'),('leiden','04_umap_leiden.png')]:
        sc.pl.umap(adata,color=col,show=False,frameon=False)
        save_figure(fig/name)
    markers=['CD3D','MS4A1','LYZ','COL1A1','PECAM1','EPCAM']
    markers=[m for m in markers if m in adata.raw.var_names]
    if markers:
        sc.pl.dotplot(adata,var_names=markers,groupby='celltype_major',use_raw=True,show=False)
        save_figure(fig/'05_marker_dotplot.png')
    # Pilot counts per selected sample; not valid for subtype-level inference.
    composition=(adata.obs.groupby(['orig.ident','subtype','celltype_major'],observed=True)
                 .size().reset_index(name='n_cells'))
    composition.to_csv(out/'pilot_cell_composition.csv',index=False)
    pd.crosstab(adata.obs['leiden'],adata.obs['celltype_major']).to_csv(out/'leiden_vs_author_celltypes.csv')
    adata.write_h5ad(out/'wu2021_pilot_processed.h5ad',compression='gzip')
    print('Completed:',adata.shape,'Figures:',fig,flush=True)

if __name__=='__main__':
    main()
