# Figure 5 — BCC & EMPD single‑cell multi‑omics

This folder contains the notebooks used to generate **Figure 5** for the manuscript, covering **Basal Cell Carcinoma (BCC)** and **Extramammary Paget’s Disease (EMPD)** analyses across scRNA‑seq, scRNA+ATAC‑seq, gene‑regulatory network inference (**scenicplus**), and lineage/trajectory analysis (**scFates**).

---

## Contents

-`1.BCC_scRNA_1.ipynb` — BCC scRNA‑seq analysis: initial QC → normalization → clustering → annotation; focuse on tumor. For Figure5a,b.

-`2.BCC_scRNA+ATAC_macs2.ipynb` — BCC scRNA+ATAC-seq analysis. For Figure5g.

-`3.BCC_scenicplus.ipynb` — BCC gene‑regulatory network inference by **scenicplus**. For Figure5e.

-`4.BCC_vs_Normal_HHKC_pseudobulk_analysis.R` — For Figure 5d.f.

-`5.EMPD_scRNA.ipynb` — EMPD scRNA‑seq: QC → processing → clustering → annotation; exports AnnData and markers. For Figure5i.

-`6.EMPD_scFates.ipynb` — EMPD lineage/trajectory analysis with **scFates**; For Figure5h.

-`7.EMPD_scRNA+ATAC_macs2.ipynb` — EMPD scRNA+ATAC-seq analysis. For Figure5n.

-`8.EMPD_scenicplus.ipynb` — EMPD gene‑regulatory network inference by **scenicplus**. For Figure5l.

-`9.EMPD_vs_Normale_Mesenchymal_analysis.r` — For Figure 5k.m.
