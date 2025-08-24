# =============================
# Pseudobulk limma-treat pipeline
# Author: YuYang
# =============================
suppressPackageStartupMessages({
  library(DaMiRseq)
  library(tidyverse)
  library(limma)
  library(ggplot2)
  library(ggsci)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(extrafont)
})

set.seed(1)                
pdf.options(family = "Arial")

folder_path <- "/home/yangyu/Project/skinAtlas/70samples/Result/1.celltype_pseudobulkDE/all_celltype"
meta_path   <- "/data/yangyu/Project/skinAtlas/70samples/sample_info.txt"

part_color_female <- c(NE = "#79a339", NF = "#a7cee2", NP = "#24467c", NT = "#e1c7a7")
part_color_male   <- c(NE = "#79a339", NT = "#e1c7a7", PS = "#9999CC")

EXCLUDE_SAMPLES <- list(
  F = character(0),
  M = c("NP16")   # only one sample
)

# 
read_expr <- function(expr_file) {
  x <- read.delim(expr_file, check.names = FALSE, row.names = 1)
  x <- x[!duplicated(rownames(x)), , drop = FALSE]
  x
}

ensure_dir <- function(...) {
  d <- file.path(...)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

safe_write_table <- function(x, file, ...) {
  if (is.null(x)) return(invisible(NULL))
  if (is.matrix(x)) x <- as.data.frame(x, check.names = FALSE)
  if (is.data.frame(x) && !nrow(x)) return(invisible(NULL))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE, ...)
}

safe_bitr <- function(genes) {
  if (!length(genes)) return(character(0))
  ids <- suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
  unique(na.omit(as.character(ids)))
}

# dt: decideTests (genes x contrasts)
# contrast_vec: ，names = colnames(dt)， "A-B"
union_up_genes_by_group <- function(dt, genes, contrast_vec) {
  stopifnot(all(colnames(dt) %in% names(contrast_vec)))
  ab <- tibble(
    colname = colnames(dt),
    pair    = unname(contrast_vec[colname]),
    A = sub("-.*$", "", pair),
    B = sub("^.*-", "", pair)
  )
  groups <- sort(unique(c(ab$A, ab$B)))
  out <- setNames(vector("list", length(groups)), groups)
  for (g in groups) {
    pos_cols <- ab$colname[ab$A == g]
    neg_cols <- ab$colname[ab$B == g]
    up_idx <- (rowSums(dt[, pos_cols, drop=FALSE] ==  1) > 0) |
              (rowSums(dt[, neg_cols, drop=FALSE] == -1) > 0)
    out[[g]] <- genes[which(up_idx)]
  }
  out
}

compare_go <- function(gene_list, out_dir, universe_entrez = NULL) {
  gene_list_f <- gene_list[vapply(gene_list, length, 1L) > 0]
  if (!length(gene_list_f)) return(invisible(NULL))

  cgo <- compareCluster(
    geneCluster = gene_list_f, fun = enrichGO,
    OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
    pAdjustMethod = "BH", pvalueCutoff = 0.2, qvalueCutoff = 0.05,
    universe = universe_entrez
  )
  cgo <- setReadable(cgo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.table(as.data.frame(cgo), file.path(out_dir, "all_group_enrich-GO.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(out_dir, "all_group_de_go_dotplot.pdf"), h = 10, w = 10); print(dotplot(cgo, label_format = 100)); dev.off()

  cgoBP <- compareCluster(
    geneCluster = gene_list_f, fun = enrichGO, ont = 'BP',
    OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
    pAdjustMethod = "BH", pvalueCutoff = 0.2, qvalueCutoff = 0.05,
    universe = universe_entrez
  )
  cgoBP <- setReadable(cgoBP, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.table(as.data.frame(cgoBP), file.path(out_dir, "all_group_enrich-GO_BP.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(out_dir, "all_group_de_go_BP_dotplot.pdf"), h = 15, w = 10); print(dotplot(cgoBP, label_format = 100, showCategory = 20)); dev.off()

  save(gene_list, cgo, cgoBP, file = file.path(out_dir, "compare_annotation.Rdata"))
}

make_named_contrasts <- function(desired_pairs, levels_present){
  keep <- vapply(strsplit(desired_pairs, "-", fixed = TRUE), function(x){ all(x %in% levels_present) }, TRUE)
  desired_pairs <- desired_pairs[keep]
  if (!length(desired_pairs)) return(NULL)
  nm <- gsub("-", "vs", desired_pairs)
  stats::setNames(desired_pairs, nm)
}

# -------------------- celltype  and  sex --------------------
run_celltype <- function(celltype, sex = c("F", "M"), keep_groups, lfc_threshold = 1.5) {
  sex <- match.arg(sex)
  message(sprintf("[ %s | %s ]", celltype, sex))

  expr_file <- file.path(folder_path, sprintf("%s_SV_adjusted_expression.txt", celltype))
  stopifnot(file.exists(expr_file))
  expr <- read_expr(expr_file)

  meta <- read_tsv(meta_path, show_col_types = FALSE) |>
  column_to_rownames("sample_id") |>
  dplyr::filter(sex == !!sex, sample_group %in% keep_groups) |>
  droplevels()

  #
  to_drop <- intersect(rownames(meta), EXCLUDE_SAMPLES[[sex]])
  if (length(to_drop)) {
    message(sprintf("[Exclude %s] %s", sex, paste(to_drop, collapse = ", ")))
    meta <- meta[setdiff(rownames(meta), to_drop), , drop = FALSE]
  }

  common <- intersect(rownames(meta), colnames(expr))
  if (length(common) < 3) {
    warning(sprintf("%s|%s: too few samples after filtering (n=%d). Skipped.", celltype, sex, length(common)))
    return(invisible(NULL))
  }
  expr <- expr[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  meta$sample_group <- factor(meta$sample_group, levels = intersect(keep_groups, unique(meta$sample_group)))

  base_dir <- ensure_dir(folder_path, celltype)
  de_dir   <- ensure_dir(base_dir, "DE")
  sex_dir  <- ensure_dir(de_dir, ifelse(sex == "F", "Female", "Male"), "one_vs_one")

  
  write.table(expr, file.path(de_dir, sprintf("%s_%s_SV_adjusted_expression.txt", ifelse(sex == "F", "Female", "Male"), celltype)),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


  group <- meta$sample_group
  if (nlevels(group) < 2) {
    warning(sprintf("%s|%s: <2 groups after filtering. Skipped DE.", celltype, sex))
    return(invisible(NULL))
  }
  design <- model.matrix(~ 0 + group); colnames(design) <- gsub("group", "", colnames(design))

  # 
  if (identical(sort(keep_groups), sort(c("NE","NF","NP","NT")))) {
    desired <- c("NE-NF","NE-NP","NE-NT","NF-NP","NF-NT","NP-NT")
  } else if (identical(sort(keep_groups), sort(c("NE","NF")))) {
    desired <- c("NE-NF")
  } else if (identical(sort(keep_groups), sort(c("NE","NT","PS")))) {
    desired <- c("NE-NT","NE-PS","NT-PS")
  } else desired <- character(0)

  named_contr <- make_named_contrasts(desired, levels(group))
  if (is.null(named_contr)) {
    warning(sprintf("%s|%s: No valid contrasts. Skipped.", celltype, sex))
    return(invisible(NULL))
  }
  contr.matrix <- makeContrasts(contrasts = unname(named_contr), levels = colnames(design))
  colnames(contr.matrix) <- names(named_contr)

  # 
  write.table(design,        file.path(sex_dir, "design_matrix.tsv"), sep="\t", quote=FALSE)
  write.table(contr.matrix,  file.path(sex_dir, "contrasts.tsv"),     sep="\t", quote=FALSE)

  # limma + treat
  vfit <- lmFit(as.matrix(expr), design)
  vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
  efit <- eBayes(vfit)
  pdf(file.path(sex_dir, "efit.pdf"), h = 5, w = 5); plotSA(efit, main = "Final model: Mean-variance trend"); dev.off()

  tfit <- treat(vfit, lfc = log2(lfc_threshold))
  dt   <- decideTests(tfit)
  write.table(data.frame(summary(dt)), file = file.path(sex_dir, "DE_summary.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

  genes <- rownames(expr); tfit$genes <- genes
  l2fc_cut <- log2(lfc_threshold)

  # 每个成对对比输出
  for (cn in colnames(contr.matrix)) {
    a <- sub("vs.*$", "", cn); b <- sub("^.*vs", "", cn)
    tt <- topTreat(tfit, coef = cn, n = Inf)
    tt$change <- as.factor(ifelse(tt$adj.P.Val <= 0.05 & tt$logFC >= l2fc_cut, a,
                           ifelse(tt$adj.P.Val <= 0.05 & tt$logFC <= -l2fc_cut, b, "NOT")))
    tt$gene <- rownames(tt)
    tt$contrast <- cn
    label_pair <- paste0(a, "_vs_", b)
    safe_write_table(tt, file = file.path(sex_dir, sprintf("%s_DE_%s.txt", celltype, label_pair)))

    sel <- unique(c(tt$gene[tt$change == a], tt$gene[tt$change == b]))
    if (length(sel) > 1) {
      hm <- expr[rownames(expr) %in% sel, , drop = FALSE]
      safe_write_table(hm, file = file.path(sex_dir, sprintf("%s_DE_adjusted_expression.txt", label_pair)))
    }
  }

  # union
  contrast_vec <- unname(named_contr); names(contrast_vec) <- colnames(contr.matrix)
  up_list <- union_up_genes_by_group(dt = dt, genes = genes, contrast_vec = contrast_vec)

  for (g in names(up_list)) {
    write.csv(up_list[[g]], file = file.path(sex_dir, sprintf("%s_gene.csv", g)), quote = FALSE, row.names = FALSE)
  }
  maxlen <- max(lengths(up_list))
  up_mat <- do.call(cbind, lapply(up_list, `length<-`, maxlen))
  up_mat[is.na(up_mat)] <- ""
  write.csv(up_mat, file = file.path(sex_dir, "common_gene.csv"), quote = FALSE, row.names = FALSE)

  common <- unique(unlist(up_list, use.names = FALSE))
  if (length(common) > 1) {
    hm_all <- expr[rownames(expr) %in% common, , drop = FALSE]
    safe_write_table(hm_all, file = file.path(sex_dir, "DE_adjusted_expression.txt"))

    # anno
    ann_df <- data.frame(sample_group = meta$sample_group)
    rownames(ann_df) <- rownames(meta)
    ann_cols <- if (sex == "F") part_color_female else part_color_male
    ann_cols <- ann_cols[names(ann_cols) %in% unique(as.character(meta$sample_group))]

    pdf(file.path(sex_dir, "pheatmap.pdf"), h = 6, w = 10)
    pheatmap(hm_all,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      scale = "row", show_rownames = FALSE,
      annotation_col   = ann_df,
      annotation_colors= list(sample_group = ann_cols))
    dev.off()

    pdf(file.path(sex_dir, "pheatmap_label.pdf"), h = 18, w = 10)
    pheatmap(hm_all,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      scale = "row", show_rownames = TRUE, fontsize_row = 5,
      annotation_col   = ann_df,
      annotation_colors= list(sample_group = ann_cols))
    dev.off()
  }

  cg <- unlist(up_list, use.names = TRUE)
  if (length(cg)) {
    write.table(data.frame(sampe_group = names(cg), gene = unname(cg)),
                file = file.path(sex_dir, "common_gene.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # 
  universe_entrez <- safe_bitr(rownames(expr))
  gene_list <- lapply(levels(meta$sample_group), function(g) safe_bitr(up_list[[g]]))
  names(gene_list) <- levels(meta$sample_group)
  compare_go(gene_list, out_dir = sex_dir, universe_entrez = universe_entrez)

  invisible(TRUE)
}

# no DEG, only save SV-adjusted expression
run_sv_only <- function(celltype) {
  expr_file <- file.path(folder_path, sprintf("%s_SV_adjusted_expression.txt", celltype))
  stopifnot(file.exists(expr_file))
  expr <- read_expr(expr_file)
  meta <- read_tsv(meta_path, show_col_types = FALSE) |>
    column_to_rownames("sample_id") |>
    dplyr::filter(sex == "M")
  #
  to_drop <- intersect(rownames(meta), EXCLUDE_SAMPLES[["M"]])
  if (length(to_drop)) {
    message(sprintf("[Exclude M] %s", paste(to_drop, collapse = ", ")))
    meta <- meta[setdiff(rownames(meta), to_drop), , drop = FALSE]
  }
  common <- intersect(rownames(meta), colnames(expr))
  if (length(common) < 2) return(invisible(NULL))
  expr <- expr[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE]
  base_dir <- ensure_dir(folder_path, celltype); de_dir <- ensure_dir(base_dir, "DE")
  write.table(expr, file.path(de_dir, sprintf("%s_%s_SV_adjusted_expression.txt", "Male", celltype)),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
}

# run DE
# Female： LFC=1.5
invisible(lapply(c('KC','KC_Channel','SGC','MEC','FB','Pc-vSMC','VEC','LEC','MEL','Lymphocyte','Mac-DC','LC'), function(ct)
  run_celltype(celltype = ct, sex = 'F', keep_groups = c('NE','NF','NP','NT'), lfc_threshold = 1.5)))

# Female：Schwann/Mast LFC=1.1 
invisible(lapply(c('Schwann','Mast'), function(ct)
  run_celltype(celltype = ct, sex = 'F', keep_groups = c('NE','NF','NP','NT'), lfc_threshold = 1.1)))

# Female：HFC/Sebocyte  only NE vs NF，LFC=1.1
invisible(lapply(c('HFC','Sebocyte'), function(ct)
  run_celltype(celltype = ct, sex = 'F', keep_groups = c('NE','NF'), lfc_threshold = 1.1)))

# Male：NE/NT/PS LFC=1.5
invisible(lapply(c('KC','KC_Channel','FB','Pc-vSMC','VEC','LEC','MEL','Lymphocyte','Mac-DC','LC','SGC'), function(ct)
  run_celltype(celltype = ct, sex = 'M', keep_groups = c('NE','NT','PS'), lfc_threshold = 1.5)))

# Male：Schwann/Mast/Sebocyte NE/NT/PS LFC=1.1
invisible(lapply(c('Schwann','Mast','Sebocyte'), function(ct)
  run_celltype(celltype = ct, sex = 'M', keep_groups = c('NE','NT','PS'), lfc_threshold = 1.1)))

# Male：HFC/MEC no DEG
invisible(lapply(c('HFC','MEC'), run_sv_only))

sink(file.path(folder_path, "sessionInfo.txt")); print(sessionInfo()); sink()

# =============================
# Notes
# - All outputs are written under: <folder_path>/<celltype>/DE/<Sex>/one_vs_one/
# - Can tweak lfc thresholds or the group sets above without touching the core logic.
# =============================
