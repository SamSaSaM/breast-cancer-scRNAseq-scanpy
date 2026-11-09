# Figure legends

**Figure 1. Quality-control metrics of the pilot scRNA-seq dataset.** Violin plots of total UMI counts (`total_counts`), detected genes (`n_genes_by_counts`), and mitochondrial transcript percentage (`percent.mito`). Individual dots represent cells.

**Figure 2. UMAP by published major cell-type annotation.** The embedding is colored by the nine broad cell populations provided in the Wu et al. dataset: B cells, CAFs, cancer epithelial cells, endothelial cells, myeloid cells, normal epithelial cells, perivascular-like cells (PVL), plasmablasts, and T cells.

**Figure 3. UMAP by breast cancer subtype label.** Cells are colored according to the subtype label of their source specimen (ER+, HER2+, or TNBC). Each subtype is represented by only one specimen in this pilot, so subtype and specimen effects cannot be separated.

**Figure 4. Unsupervised Leiden clustering.** The UMAP is colored by the 17 Leiden clusters computed in the Scanpy pilot workflow. Multiple clusters can map to the same published broad cell type.

**Figure 5. Canonical marker gene expression.** Dot plot of CD3D, MS4A1, LYZ, COL1A1, PECAM1, and EPCAM across published major cell types. Dot size denotes the fraction of expressing cells; color reflects the expression summary produced by Scanpy.

**Figure 6. Leiden cluster concordance with published cell-type annotations.** Heatmap of the within-cluster fraction assigned to each published major cell type. Each row sums to approximately one, subject to rounding. This measures agreement with provided labels, not independent annotation accuracy.

**Figure 7. UMAP by specimen identity.** The same UMAP colored by `orig.ident` (CID3586, CID4471, CID44971), showing shared and specimen-enriched regions.

**Figure 8. Sample composition within Leiden clusters.** Heatmap of the within-cluster fraction of cells from each specimen. Each row sums to approximately one, subject to rounding. Sample enrichment cannot by itself distinguish biological heterogeneity from technical batch effects.
