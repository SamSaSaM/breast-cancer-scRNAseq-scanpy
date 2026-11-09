#!/usr/bin/env python3
"""Inspect the Wu et al. 2021 matrix and metadata without loading expression counts."""
import argparse
from pathlib import Path
import pandas as pd
from scipy.io import mminfo

FILES = ['count_matrix_sparse.mtx','count_matrix_barcodes.tsv','count_matrix_genes.tsv','metadata.csv']

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, default=Path.home()/'Downloads'/'Wu_etal_2021_BRCA_scRNASeq')
    args = p.parse_args()
    d = args.data_dir.expanduser()
    for f in FILES:
        if not (d/f).is_file():
            raise FileNotFoundError(f'Missing: {d/f}')
    rows, cols, nnz, fmt, field, sym = mminfo(d/FILES[0])
    barcodes = pd.read_csv(d/FILES[1], sep='\t', header=None, dtype=str)[0]
    genes = pd.read_csv(d/FILES[2], sep='\t', header=None, dtype=str)[0]
    meta = pd.read_csv(d/FILES[3], index_col=0, low_memory=False)
    print(f'Matrix: {rows:,} x {cols:,}; nnz={nnz:,}; format={fmt}; field={field}; symmetry={sym}')
    print(f'Barcodes: {len(barcodes):,}; genes: {len(genes):,}; metadata: {meta.shape}')
    print('Matrix orientation:', 'genes x cells' if (rows,cols)==(len(genes),len(barcodes)) else 'UNKNOWN — STOP')
    if (rows,cols)!=(len(genes),len(barcodes)):
        raise ValueError('Matrix dimensions do not match genes x barcodes')
    if barcodes.duplicated().any() or meta.index.duplicated().any():
        raise ValueError('Duplicate barcode identifiers')
    print('Metadata columns:', meta.columns.tolist())
    print('Barcode overlap:', len(set(barcodes).intersection(meta.index)), '/', len(barcodes))
    print('Sample counts:')
    print(meta.groupby(['subtype','orig.ident']).size().to_string())
    print('Cell types:')
    print(meta['celltype_major'].value_counts().to_string())

if __name__ == '__main__':
    main()
