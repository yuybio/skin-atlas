This **folder** implements an end‑to‑end **scRNA‑seq** workflow from raw reads to downstream multi‑analysis. Follow the numbered files in order to:

1. Run Cell Ranger to align FASTQs and generate count matrices.
2. Import raw matrices and perform basic QC/processing2–3) Remove ambient RNA contamination with **DecontX** and write back adjusted counts
3. Cross‑batch integration with **scVI.**
4. **Pseudobulk** differential expression gene with limma‑voom.
5. **CellPhoneDB** ligand–receptor analysis at the cell‑type and sub-type level.
6. Multicellular **super-module** analysis.

---

# `0.run_cellranger_script.sh`

**Purpose**: Run `cellranger count` to produce gene×cell matrices.

**Output**: For each sample, `outs/filtered_feature_bc_matrix`.

---

# `1.Data_import and Processing.ipynb`

**Purpose**: Import Cell Ranger matrices, perform basic QC (nFeature/nCount/mito%), initial t‑SNE/UMAP, and save the object.

---

# `2.ambientRNA removede DecontX singleSample.R`

**Purpose**: Run **DecontX** per sample to estimate and annotate ambient RNA contamination.

**Output**: Object with `decontX_contamination` h5ad

---

# `3.Basic Process for the remove_ambientRNA_count.ipynb`

**Purpose**: Use DecontX contamination estimates to adjust the count matrix (e.g., per‑cell proportion subtraction) and re‑QC.

---

# `4.scVI integrated.ipynb`

**Purpose**: Batch integration with **scVI** to learn a shared latent space and corrected expression.

**Steps**: AnnData → `scvi.model.SCVI.setup_anndata` → train → `get_latent_representation` → UMAP → save.

---

# `5.Pseudobulk_DE_limma.R`

**Purpose**: Aggregate counts to **pseudobulk** per cell type × sample and run `limma‑voom` DE analysis.

---

# `6.run_cpdb_method2.py`

**Purpose**: Run **CellPhoneDB** (`method=statistical_analysis`) at the **cell‑type** and **sub-type** level.

# `7.celltype_ligand.R`

**Purpose**: Perform Ligand–target links from **NicheNet**.

# 8. Super-Module/--- **Codes and scripts to perform gene-set-features coordination (super-module) analysis.**

describes how to reproduce the **gene‑set features coordination (super‑module)** analysis, from building module–expression matrices and correlations to cell–cell communication enrichment and network visualizations. All paths are **relative to this folder** unless noted.

## 0) Build the C++ tools

Compile C++ codes using GCC compiler:

```bash
cd src
make -f Makefile
```

Resulting executables used below (source files in `src/`):

- `generate_module_sc_matrix`  — build module expression matrix from per‑cell‑type CSVs (`generate_module_sc_matrix.cpp`)
- `filter_matrix`              — filter correlation matrices by thresholds and node lists (`filter_matrix.cpp`)
- `ana_cellphone_db_results`   — parse CellPhoneDB outputs to significant L–R tables (`ana_cellphone_db_results.cpp`)
- `integrate_gene_program`     — integrate modules/markers with L–R and ligand–target links (`integrate_gene_program.cpp`)
- `combine_link_metrics`       — combine multiple link metrics to shared‑target tables (`combine_link_metrics.cpp`)
- `generate_dot_cotargets`     — build DOT graphs for co‑targets (`generate_dot_cotargets.cpp`)

---

## Part I. Basic analysis

### Step 1 — Calculate ADEG‑module expression matrix

**Input:** `modulefile` (a text file listing `*.module_mean_exp.csv`), one file per **cell type** (Female/Male). Example: `Female_KC_module_mean_exp.csv`.
**Executable:** `generate_module_sc_matrix ` (C++ compiled; source code: src/ generate_module_sc_matrix.cpp)
**Output:** `module_exp_matrix.txt`

```bash
generate_module_sc_matrix modulefile module_exp_matrix.txt
```

### Step 2 — Prepare cell‑type/subtype proportion matrices

Provide the following text files:

- `celltype_sample_proportion.txt`
- `granular_celltype_sample_proportion.txt`

###  Step 3 — Calculate pair-wise correlation between all possible pairs of values from cell-type proportion and ADEG-module expression.

**Input:** `module_exp_matrix.txt`, `celltype_sample_proportion.txt`, `granular_celltype_sample_proportion.txt`

**Script (Python):** `cal_correlation.py`

**Output:**

- `moduleexp_composition_correlation.txt` (correlation coefficients)
- `moduleexp_composition_correlation_p.txt` (P values)

```bash
python cal_correlation.py
```

### Step 4 — Draw heatmap of correlation matrix with hierarchical clustering using R script.

**R script:**

```r
library(pheatmap)
tb <- read.table("moduleexp_composition_correlation.txt", header=TRUE, sep="	", row.names=1)
pheatmap(as.matrix(tb))
```

### Step 5 — Define super‑modules

From the hierarchical clustering results of correlation matrix, two highly correlated clusters/modules were selected as super-modules, named: **G1 (Super-module #1)** and **G2 (Super-module #2)**.

### Step 6 — Correlation networks for G1 and G2 (main Fig. 3a and e)

#### 6a) Filter edges by thresholds

get correlation coefficient from correlation matrix as follows : **Keep edges** within G1/G2 with **correlation > 0.5** and **P < 0.001**.
**Inputs:** `moduleexp_composition_correlation.txt`, `moduleexp_composition_correlation_p.txt`, `integrated_G1.txt` / `integrated_G2.txt`
**Executable:** `filter_matrix` (C++ compiled; source code: src/filter_matrix.cpp)
**Outputs:** `G1.link.c0.5.p0.001.txt`, `G1.ver.c0.5.p0.001.txt`, `G2.link.c0.5.p0.001.txt`, `G2.ver.c0.5.p0.001.txt`

```bash
# G1
filter_matrix moduleexp_composition_correlation.txt   moduleexp_composition_correlation_p.txt   integrated_G1.txt   0.5 0.001 G1.link.c0.5.p0.001.txt G1.ver.c0.5.p0.001.txt

# G2
filter_matrix moduleexp_composition_correlation.txt   moduleexp_composition_correlation_p.txt   integrated_G2.txt   0.5 0.001 G2.link.c0.5.p0.001.txt G2.ver.c0.5.p0.001.txt
```

> **Notes:**
> • For the vertex (`*.ver.*`) files, add two columns: **cell type** (e.g., KC, FB) and **feature type** (C = composition proportion; M = ADEG‑module).
> • For visualization, you may remove the connections **within the same broad same types in link files**.

#### 6b) Draw correlation network plots. (R)

```r
library(igraph)
library(RColorBrewer)
library(tidyverse)

links <- read.table("G1.link.c0.5.p0.001_removesameCT.txt", sep="	", header=TRUE) |> as.data.frame()
nodes <- read.table("G1.ver.c0.5.p0.001.txt",            sep="	", header=TRUE) |> as.data.frame()

network <- graph_from_data_frame(d=links, vertices=nodes, directed=FALSE)
palette <- colorRampPalette(brewer.pal(n=7, name="YlOrRd"))(100)

layout <- layout_in_circle(network, order=nodes$name)

# customize colors by cell type
G1_cols <- c("#1F77B4","#FF7F0E","#279E68","#D62728",
             "#E377C2","#FFBB78","#FF9896","#C49C94")
node_cols <- G1_cols[as.numeric(factor(V(network)$celltype,
                         levels=c("KC","Channel","SGC","HFC",
                                  "FB","MEL","Imm","LC")))]

shapes <- c("circle","rectangle")
node_shapes <- shapes[as.numeric(as.factor(V(network)$feature))]

plot(network,
     layout=layout,
     vertex.size=V(network)$aggsize,
     vertex.color=node_cols,
     vertex.label=NA,
     vertex.shape=node_shapes,
     edge.width=pmax(E(network)$cor*10-5, 0.5),
     edge.color=palette[pmin(pmax(round(E(network)$cor*100),1),100)],
     edge.curved=0.3)
```

### Step 7 — Use boxplot to show relative expression of elements in super-modules at different anatomical sites (main Fig. 3a and e)

#### 7a) Build expression tables per group

**Inputs:** `module_exp_matrix.txt`, `celltype_sample_proportion.txt`, `granular_celltype_sample_proportion.txt`, `G1.ver.c0.5.p0.001.txt` / `G2.ver.c0.5.p0.001.txt`
**Script (Python):** `filter_exp_matrix_for_group.py`
**Outputs:** `G1.ver.c0.5.p0.001.exptable.txt`, `G2.ver.c0.5.p0.001.exptable.txt`

```bash
python filter_exp_matrix_for_group.py G1.ver.c0.5.p0.001.txt G1.ver.c0.5.p0.001.exptable.txt
python filter_exp_matrix_for_group.py G2.ver.c0.5.p0.001.txt G2.ver.c0.5.p0.001.exptable.txt
```

#### 7b) Z‑score transform and ggplot‑ready tables

```bash
python z_trans.py G1.ver.c0.5.p0.001.exptable.txt G1.ver.c0.5.p0.001.exptable.z.txt
python cal_zs_ave_gen_table.py G1.ver.c0.5.p0.001.exptable.z.txt G1.ver.c0.5.p0.001.aveexp_zs.table.txt
# Repeat for G2 by replacing the prefix.
```

#### 7c) Draw boxplots (R)

```r
library(ggplot2)
tb <- read.table("G1.ver.c0.5.p0.001.aveexp_zs.table.txt", sep="	", header=TRUE, row.names=1) |> as.data.frame()
tb$gender   <- metadata$gender[idx]   # ensure 'metadata' and 'idx' exist in your environment
tb$newgroup <- paste(tb$gender, tb$site, sep="_")

ggplot(tb, aes(x=newgroup, y=z)) +
  geom_boxplot(fill="#5780AD") +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_x_discrete(limits=c("F_NF","F_NP","F_NE","F_NT","M_PS","M_NE","M_NT")) +
  theme_classic()
```

> Repeat for **G2** using its corresponding tables.

---

## Part II. Cell–cell communication analysis

### Step 8 — Prepare L–R tables from CellPhoneDB results(extract and format tables): Significant ligand-receptor pairs with cell type/subtype information were loaded.

**Pre‑requisite:** Run CellPhoneDB on your scRNA‑seq data.
**Inputs:** `statistical_analysis_significant_means_all.txt`, `statistical_analysis_interaction_scores_all.txt`, `statistical_analysis_deconvoluted_all.txt`
**Executable:** `ana_cellphone_db_results` (C++ compiled; source code: src/ ana_cellphone_db_results.cpp)
**Output:** `all_sig_lrp.txt`

```bash
ana_cellphone_db_results   statistical_analysis_significant_means_all.txt   statistical_analysis_interaction_scores_all.txt   statistical_analysis_deconvoluted_all.txt   all_sig_lrp.txt
```

### Step 9 — Calculate enrichment of L-R pairs within any two pairs of gene-set-features within super-module. Perform a Monta-Carlo simulation to estimate significance.

**Inputs:** `all_conserved_markers.txt`, `all_module.txt`, `G1.ver.c0.5.p0.001.txt`, `G1.link.c0.5.p0.001.txt`, `all_conserved_markers.unique.txt`, `all_module_genes.unique.txt`, `all_sig_lrp.txt`
**Outputs:** `G1.L2R.P_matrix.txt`, `G1.L2R.list.txt`, `G1.link_L2R.txt`
**Executable:** `integrate_gene_program` (C++ compiled; source code: src/ integrate_gene_program.cpp)

```bash
integrate_gene_program -C all_conserved_markers.txt   -M all_module.txt   -N G1.ver.c0.5.p0.001.txt   -f siglrp   -O G1   -L G1.link.c0.5.p0.001.txt   -n all_conserved_markers.unique.txt   -m all_module_genes.unique.txt   -l all_sig_lrp.txt
```

#### Step 9a— Draw heatmap to show significance (P-value) of enrichment of L-R pairs within pairs in super-module.

```r
library(pheatmap); library(RColorBrewer)
tb  <- read.table("G1.L2R.P_matrix.txt", sep="	", header=TRUE, row.names=1)
col <- colorRampPalette(brewer.pal(n=7, "OrRd"))(100)
bk  <- seq(0, 3, by=0.03)

pheatmap(as.matrix(tb), color=col, breaks=bk, scale="none",
         cluster_rows=FALSE, cluster_cols=FALSE)

# Focus on selected rows/columns (edit lists to match your labels)
frow <- c("Basal.Ker.ADGRL3","Spinous.Ker.ADGRL3",
          "Female_KC_M1","Female_KC_M4","HFC","HF-MFCs","HF-MSCs","HF-PCs",
          "Mesenchymal.Fb.COCH","Mesenchymal.Fb.INHBA","Female_FB_M5","Mel_1","LC","LC1","MigLC")
fcol <- c("Basal.Ker.ADGRL3","Spinous.Ker.ADGRL3","Female_KC_M1",
          "Female_KC_Channel_M1","HFC","HF.MFCs",
          "Mesenchymal.Fb.COCH","Mesenchymal.Fb.INHBA","Mel_1","LC","LC1","MigLC")

tb_f <- tb[frow, fcol]
pheatmap(as.matrix(tb_f), color=col, breaks=bk, scale="none",
         cluster_rows=FALSE, cluster_cols=FALSE)
```

> Repeat for **G2** (Super‑module #2).

### Step 10 — Ligand–target links from **NicheNet**

**10a).** For each gene‑set feature, we collect predicted **active ligand–target** links and generate a table file named `XXX_active_ligand_target_links.csv` (XXX represents gene-set feature, e.g., `Female_FB_M5_active_ligand_target_links.csv`). For G1/G2, generate a text file recording all the `*.csv` files, named `ltfile`.

**10b).** Get ligand-targets for each pair in super-module (e.g., G1)
**Inputs:** `all_conserved_markers.txt`, `all_module.txt`, `G1.ver.c0.5.p0.001.txt`, `G1.link.c0.5.p0.001.txt`, `all_conserved_markers.unique.txt`, `all_module_genes.unique.txt`, `ltfile`, `all_sig_lrp.txt`
**Outputs:** `G1.ligand_target_ratio.txt`, `G1.sorted_ligand.txt`, `G1.target_ratio.txt`, `G1.sorted_target.stat.txt`, `G1.ligand_count.txt`, `G1.G*.dot` (network files in DOT format)
**Executable:** `integrate_gene_program` (C++ compiled; source code: src/ integrate_gene_program.cpp)

```bash
integrate_gene_program -C all_conserved_markers.txt   -M all_module.txt   -N G1.ver.c0.5.p0.001.txt   -L G1.link.c0.5.p0.001.txt   -f ligandtarget   -O G1   -n all_conserved_markers.unique.txt   -m all_module_genes.unique.txt   -t ltfile   -l all_sig_lrp.txt
```

> DOT format files can be used to draw network graph using dot from Graphviz program (See [https://graphviz.org/download/](https://graphviz.org/download/) to download and install Graphviz)
> `dot -Tsvg G1.G1.dot > G1.G1.svg`

### Step 11 — Calculate shared genes within pairs in super-module.

**Inputs:** `all_conserved_markers.txt`, `all_module.txt`, `G1.ver.c0.5.p0.001.txt`, `G1.link.c0.5.p0.001.txt`, `all_conserved_markers.unique.txt`, `all_module_genes.unique.txt`
**Outputs:** `G1.link_sig_overlappedgene`, `G1.link_overlapped_gene.txt`
**Executable:** `integrate_gene_program` (C++ compiled; source code: src/ integrate_gene_program.cpp)

```bash
integrate_gene_program -C all_conserved_markers.txt   -M all_module.txt   -N G1.ver.c0.5.p0.001.txt   -f linkgeneoverlap   -O G1   -L G1.link.c0.5.p0.001.txt   -n all_conserved_markers.unique.txt   -m all_module_genes.unique.txt
```

> Repeat for **G2**.

### Step 12 — Get shared targets from ligand within super-modules.

**12a).** Get shared targets from ligand from results of Steps 9-11.
**Inputs:** `G1.link.c0.5.p0.001.txt`, `G1.link_overlapped_gene.txt`, `G1.ligand_target_ratio.txt`
**Output:** `G1.ovlg.sharedtarget.txt`
**Executable:** `combine_link_metrics` (C++ compiled, source code: src/ combine_link_metrics.cpp)

```bash
combine_link_metrics G1.link.c0.5.p0.001.txt G1.link_overlapped_gene.txt   G1.ligand_target_ratio.txt G1.ovlg.sharedtarget.txt tmp
```

**12b.** Generate network for shared targets(DOT graphs):
**Inputs:** `all_conserved_markers.txt`, `all_module.txt`, `G1.ver.c0.5.p0.001.txt`, `G1.link.c0.5.p0.001.txt`, `G1.ovlg.sharedtarget.txt`, `G1.sorted_ligand.txt`
**Output prefix:** `G1-cotarget.G*.dot`
**Executable:** `generate_dot_cotargets` (C++ compiled; source code: src/generate_dot_cotargets.cpp)

```bash
generate_dot_cotargets all_conserved_markers.txt all_module.txt   G1.ver.c0.5.p0.001.txt G1.link.c0.5.p0.001.txt   G1.ovlg.sharedtarget.txt G1.sorted_ligand.txt G1-cotarget

# Render with Graphviz
dot -Tsvg G1-cotarget.G0.dot > G1-cotarget.G0.svg
```
