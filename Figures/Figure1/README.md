# Figure 1

This folder contains the notebooks/script to **generate Figure 1** for the manuscript, focused on  consensus cell types across all anatomical sites.

## Contents

### `Figure1_a.ipynb` —  **Broad cell types across all anatomical sites of skin**.

**Input :** Integrated data after batch removal.(atals_Granular.h5ad)

**Output :** UMAP plot of all broad cell type and dotplot of marker genes.

### `Figure1_b_plot_dendrogram.r` — generates a dendrogram of cell types based on top 10 marker gene expression

**Input :** The gene expression matrix of top 10 marker genes and meta data of all cells.

**Output :** The dendrogram of all sub types.

### `Figure1_c.ipynb` — The 8 sub types of Keratinocyte.

**Input :** The anndata of Keratinocyte（kc.h5ad）

**Output :** UMAP plot of 8 sub type for Keratinocyte and dotplot of marker genes.
