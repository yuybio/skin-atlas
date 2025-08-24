# Figure 2 — Site heterogeneity metrics & gene modules

This folder contains notebooks and scripts to **generate and assemble Figure 2** panels focused on **anatomical-site heterogeneity** (e.g., site-biased genes, clustering quality, and PCA views).

---

## Contents

### `1. cal_ratio.py`— **Anatomically differential composition analysis for Figure2_a.**

For broad cell types, the abundance was defined as the percentage relative to the whole cell population.

For subtypes of a broad cell type, the abundance is defined as the percentage relative to the entire broad cell type population.

**Input :** The metadata of all cells (adata.obs)

**Output :** The csv report and the heatmap of the relative abundance of broad cell types and sub cell types.

### `2.cal_Silhouette_Coefficient_site.ipynb` — computes **Silhouette Coefficient by site**.

**Input :** The pseudo-bulk gene expression profiles for each sample.

**Output :** The silhouette coefficient for each sample.

### `3.plot_pca_3d.R` — renders a **3D PCA**  for site-level separation.

**Input :** The pseudo-bulk gene expression profiles for each sample(same as 2.)

**Output :** PCA plot.

### `4.run_site_gene_module.py` — computes **site gene modules** (e.g., aggregates ASBGs per site, calculates module scores per cell/sample) and exports tables.

**Input :** The pseudo-bulk difference gene expression profiles for each sample.

**Output :** Hotspot_gene_modules.csv

---

## For figure2_f-h, the scripts used to identify super-modules are located in `processing/Super-Module`
