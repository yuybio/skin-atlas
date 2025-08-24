# Cross-anatomical single-cell Skin Atlas

This repository contains preprocessing and integration scripts for single-cell RNA-seq and single-cell multiome analysis. The results support the figures presented in the manuscript.

---

## 📁 Project Structure

```
├── Figures
│   ├── Figure1
│   ├── Figure2
│   ├── Figure3
│   ├── Figure4
│   ├── Figure5
│   └── Figure6
├── LICENSE
├── preprocessing
│   ├── 0.run_cellranger_script.sh
│   ├── 1.Data_import_and_Processing.ipynb
│   ├── 2.ambientRNA_removede_DecontX_singleSample.R
│   ├── 3.Basic_Process_for_the_remove_ambientRNA_count.ipynb
│   ├── 4.scVI_integrated.ipynb
│   ├── 5.Pseudobulk_DE_limma.R
│   ├── 6.run_cpdb_method2.py
│   ├── 7.celltype_ligand.R
│   ├── 8.Super-Module
│   └── README.md
├── README.md
└── requirements.txt
```

---

## 🖥️ System Requirements

- Operating System: CentOS 7.9
- Python ≥ 3.9
- R ≥ 4.2
- GPU is optional for scVI but highly recommended for performance

---

## 📦 Requirements

All required Python and R packages are listed in `requirements.txt` and can be installed using pip, CRAN, or Bioconductor as appropriate.

---

## 📜 License

This project is released under the **MIT License**. See [LICENSE](./LICENSE) for details.
