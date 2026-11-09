#!/usr/bin/env python3
"""Stream Wu et al. Matrix Market counts to construct a small cells x genes sparse AnnData."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.io import mminfo

DEFAULT_DATA = Path.home()/'Downloads'/'Wu_etal_2021_BRCA_scRNASeq'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    p.add_argument('--output-dir', type=Path, default=None)
    p.add_argument('--cells-per-sample', type=int, default=2500)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--chunk-size', type=int, default=500000)
    args = p.parse_args()
    data = args.data_dir.expanduser()
    out = args.output_dir.expanduser() if args.output_dir else data/'analysis'
    out.mkdir(parents=True, exist_ok=True)
    if args.cells_per_sample < 1 or args.chunk_size < 1:
        p.error('cells-per-sample and chunk-size must be positive')
    meta = pd.read_csv(data/'metadata.csv', index_col=0, low_memory=False)
    meta.index = meta.index.astype(str)
    barcodes = pd.read_csv(data/'count_matrix_barcodes.tsv', sep='\t', header=None, dtype=str)[0].to_numpy()
    genes = pd.read_csv(data/'count_matrix_genes.tsv', sep='\t', header=None, dtype=str)[0].to_numpy()
    n_genes,n_cells,nnz,fmt,field,sym = mminfo(data/'count_matrix_sparse.mtx')
    if (n_genes,n_cells)!=(len(genes),len(barcodes)) or fmt!='coordinate' or sym!='general':
        raise ValueError('Unexpected Matrix Market dimensions or format')
    if len(set(barcodes)) != len(barcodes) or meta.index.has_duplicates:
        raise ValueError('Duplicate barcode IDs')
    meta = meta.reindex(barcodes)
    if meta[['orig.ident','subtype','celltype_major']].isna().any().any():
        raise ValueError('Missing required metadata for one or more barcodes')
    rng = np.random.default_rng(args.seed)
    selected = []
    for subtype in ['ER+','HER2+','TNBC']:
        candidates = meta.loc[meta['subtype'].eq(subtype)]
        if candidates.empty:
            raise ValueError(f'No cells for subtype {subtype}')
        sample_id = candidates['orig.ident'].value_counts().index[0]
        cells = candidates.index[candidates['orig.ident'].eq(sample_id)].to_numpy()
        chosen = rng.choice(cells, size=min(args.cells_per_sample,len(cells)), replace=False)
        selected.extend(chosen)
        print(f'{subtype}: {sample_id}; selected {len(chosen)} cells', flush=True)
    selected_mask = np.isin(barcodes, selected)
    indices = np.flatnonzero(selected_mask)
    names = barcodes[indices]
    lookup = np.full(n_cells,-1,dtype=np.int32)
    lookup[indices] = np.arange(len(indices),dtype=np.int32)
    # A Matrix Market coordinate file has a header, optional % comments, then dimensions.
    # Count the header lines robustly instead of assuming exactly three lines.
    with open(data/'count_matrix_sparse.mtx','rt') as f:
        header = f.readline().strip()
        if not header.lower().startswith('%%matrixmarket matrix coordinate'):
            raise ValueError('Not a coordinate Matrix Market file')
        skiprows = 1
        while True:
            line = f.readline()
            if not line:
                raise ValueError('Missing Matrix Market dimensions')
            skiprows += 1
            if not line.startswith('%'):
                dims = list(map(int,line.split()))
                if dims != [n_genes,n_cells,nnz]:
                    raise ValueError('Matrix Market header dimensions disagree')
                break
    row_parts, col_parts, val_parts = [], [], []
    reader = pd.read_csv(data/'count_matrix_sparse.mtx', sep=r'\s+', skiprows=skiprows,
                         header=None, names=['gene','cell','count'],
                         dtype={'gene':np.int32,'cell':np.int32,'count':np.int32},
                         chunksize=args.chunk_size)
    for i,chunk in enumerate(reader,1):
        source_cols = chunk['cell'].to_numpy()-1
        if source_cols.min()<0 or source_cols.max()>=n_cells:
            raise ValueError('Matrix cell index out of bounds')
        mapped = lookup[source_cols]
        keep = mapped>=0
        if keep.any():
            row_parts.append(chunk['gene'].to_numpy()[keep]-1)
            col_parts.append(mapped[keep])
            val_parts.append(chunk['count'].to_numpy()[keep])
        if i%50==0:
            print(f'Scanned approximately {i*args.chunk_size:,} entries',flush=True)
    if not row_parts:
        raise RuntimeError('No matrix entries matched selected cells')
    rows = np.concatenate(row_parts)
    cols = np.concatenate(col_parts)
    vals = np.concatenate(val_parts)
    if rows.min()<0 or rows.max()>=n_genes or (vals<0).any():
        raise ValueError('Invalid gene index or negative count')
    X = sparse.coo_matrix((vals,(cols,rows)),shape=(len(indices),n_genes),dtype=np.int32).tocsr()
    X.sum_duplicates()
    del row_parts,col_parts,val_parts,rows,cols,vals
    adata = ad.AnnData(X=X, obs=meta.loc[names].copy(), var=pd.DataFrame(index=pd.Index(genes,dtype=str)))
    adata.obs_names = pd.Index(names)
    adata.var_names_make_unique()
    adata.uns['source'] = 'Wu et al. 2021 GSE176078; sampled pilot, not full cohort'
    adata.uns['sampling_seed'] = args.seed
    output = out/'wu2021_pilot_raw.h5ad'
    adata.write_h5ad(output,compression='gzip')
    print(f'Saved: {output}; shape={adata.shape}; nnz={adata.X.nnz:,}',flush=True)

if __name__=='__main__':
    main()
