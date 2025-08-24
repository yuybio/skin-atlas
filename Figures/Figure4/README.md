# Figure 4 — Integrating skin GWAS data with the cross-anatomical skin atlas

This folder contains notebooks and scripts used to generate the panels for **Figure 4**, focusing on **skin-related genetic traits and disorders**.

---

## Contents

- `magma/` — The gene-level p-values/z-scores for a given trait (pval_file/zscore_file) using **MAGMA** from GWAS summary statistics.
- `score_subsample/` — scDRS score file.
- `gwas_analysis.ipynb` — main notebook.

  - **Input** : scDRS full score file(score_subsample/) and cell annotations stored in `adata.obs` of `atlas_Granular.h5ad`
  - **Output** : test statistics/p-values based on the scDRS MC tests.
