#!/usr/bin/env python3
"""Compare Leiden clusters with author annotations and sample identities."""
from pathlib import Path
import argparse
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_DATA = Path.home() / 'Downloads' / 'Wu_etal_2021_BRCA_scRNASeq'

def heatmap(table, output, title, xlabel, figsize=(12, 8)):
    plt.figure(figsize=figsize)
    sns.heatmap(table, cmap='Blues', vmin=0, vmax=1, annot=True, fmt='.2f')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Leiden cluster')
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    analysis = args.data_dir / 'analysis'
    figs = analysis / 'figures'
    figs.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(analysis / 'wu2021_pilot_processed.h5ad')

    counts = pd.crosstab(adata.obs['leiden'], adata.obs['celltype_major'])
    proportions = counts.div(counts.sum(axis=1), axis=0)
    counts.to_csv(analysis / 'leiden_celltype_counts.csv')
    proportions.to_csv(analysis / 'leiden_celltype_proportions.csv')
    heatmap(proportions, figs / '06_leiden_annotation_heatmap.png',
            'Leiden Clusters vs Author Cell Types', 'Author cell type')

    sc.pl.umap(adata, color='orig.ident', show=False)
    plt.savefig(figs / '07_umap_sample.png', dpi=300, bbox_inches='tight')
    plt.close('all')

    sample_counts = pd.crosstab(adata.obs['leiden'], adata.obs['orig.ident'])
    sample_props = sample_counts.div(sample_counts.sum(axis=1), axis=0)
    sample_counts.to_csv(analysis / 'leiden_sample_counts.csv')
    sample_props.to_csv(analysis / 'leiden_sample_proportions.csv')
    heatmap(sample_props, figs / '08_leiden_sample_heatmap.png',
            'Sample Composition of Leiden Clusters', 'Sample ID', figsize=(10, 7))
    print(f'Validation figures saved to {figs}')

if __name__ == '__main__':
    main()
