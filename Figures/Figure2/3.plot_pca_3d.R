# =============================
# Plot 3D PCA for different cell types
# Author: YuYang
# =============================

library(dplyr)
library(tidyverse)
library(plot3D)
library(extrafont)
pdf.options(family = "Arial")
set.seed(1111)

# 
meta_path <- "/data/yangyu/Project/skinAtlas/70samples//sample_info.txt"

GROUP_COLOR_MAP <- c(
  'NE' = '#79a339',
  'NP' = '#24467c',
  'NT' = '#e1c7a7',
  'NF' = '#a7cee2',
  'PS' = '#9999CC'
)
GROUP_LEVELS <- c('NF','NT','NE','PS','NP')

# return list(pca = prcomp, pca_data = data.frame, axis_labels = PC)
build_pca_data <- function(expr_file,
                           meta_file = meta_path,
                           group_levels = GROUP_LEVELS,
                           group_color_map = GROUP_COLOR_MAP) {
  data_adjust <- read.table(expr_file, check.names = FALSE)
  meta_data <- readr::read_tsv(meta_file, show_col_types = FALSE) %>%
    tibble::column_to_rownames("sample_id")
  common <- intersect(colnames(data_adjust), rownames(meta_data))
  if (length(common) < 3) stop(paste0("样本重叠数不足 3：", expr_file))
  data_adjust <- data_adjust[, common, drop = FALSE]
  meta_data   <- meta_data[common, , drop = FALSE]

  if (!"sample_group" %in% colnames(meta_data)) {
    stop("meta_data 中缺少 sample_group 列。")
  }
  meta_data$sample_group <- factor(meta_data$sample_group, levels = group_levels)

  # ）
  gender_col <- if ("sex" %in% colnames(meta_data)) "sex" else if ("gender" %in% colnames(meta_data)) "gender" else NA
  gender_vec <- if (!is.na(gender_col)) as.character(meta_data[[gender_col]]) else rep(NA_character_, nrow(meta_data))
  sex_shape  <- ifelse(gender_vec == "M", 16, 17)

  # 
  pca <- prcomp(t(as.data.frame(data_adjust)), scale. = FALSE)
  pvar <- round((pca$sdev^2)/sum(pca$sdev^2)*100, 1)
  axis_labels <- paste0(colnames(pca$x), " (", pvar, "%)")

  # 
  group_chr   <- as.character(meta_data$sample_group)
  group_color <- unname(group_color_map[group_chr]); group_color[is.na(group_color)] <- "#BDBDBD"

  pca_data <- data.frame(
    Sample = rownames(pca$x),
    Classification = meta_data$sample_group,
    gender = gender_vec,
    X = pca$x[,1], Y = pca$x[,2], Z = pca$x[,3],
    group_color = group_color,
    sex_shape = sex_shape,
    stringsAsFactors = FALSE
  )
  list(pca = pca, pca_data = pca_data, axis_labels = axis_labels)
}

# ==== fun to build pca data ====
build_pca_data_list <- function(named_paths) {
  stopifnot(is.character(named_paths))
  out <- lapply(named_paths, function(p) build_pca_data(p))
  names(out) <- names(named_paths)
  out
}
expr_paths <- c(
  "whole_skin" = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/4.scVIintegrated/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "FB"         = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/5.Fibroblasts/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "KC"         = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/23.KC_VIMKC///0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "KC_Channel" = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/25.KC_gap///0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "SGC"        = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/7.Sweat_gland/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "HFC"        = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/14.Hair_follicle_cells/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "Sebocyte"   = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/21.Sebocytes//0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "MEC"        = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/20.MEC//0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "Pc"         = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/10.Pc-vSMC//0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "VEC"        = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/9.VEC/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "LEC"        = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/19.LEC/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "MEL"        = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/6.Melanocytes/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "Schwann"    = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/8.Schwann//0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "Tcell"      = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/11.T//0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "MonoMac"    = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/18.MonoMac/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "LC"         = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/17.LC/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  "Mast"       = "/data/yangyu/Project/skinAtlas/70samples/Analysis/results/16.Mast/0.pseudobulk/rawcount/DaMiRseq_adjust/SV_adjusted_expression.txt",
  # sub-celltypes
  "Pro-inflammatory.Fb" = "/home/yangyu/Project/skinAtlas/70samples/Result/2.sub_celltype_pseudobulkDE/all_subcelltype/Pro-inflammatory.Fb_SV_adjusted_expression.txt",
  "Basal.Ker.COL17A1"   = "/home/yangyu/Project/skinAtlas/70samples/Result/2.sub_celltype_pseudobulkDE/all_subcelltype/Basal.Ker.COL17A1_SV_adjusted_expression.txt",
  "MLKC"        = "/home/yangyu/Project/skinAtlas/70samples/Result/2.sub_celltype_pseudobulkDE/all_subcelltype/Hillock-Club.Epi_SV_adjusted_expression.txt"
)

# ====  pca_data list ====
res_list <- build_pca_data_list(expr_paths)


# whole skin level 
pdf('whole_skin_3d_pca.pdf', width = 8, height = 8)
pca.data <- res_list$whole_skin$pca_data
scatter3D(
  x = pca.data$X, y = pca.data$Y, z = pca.data$Z,
  colvar = as.numeric(factor(pca.data$Classification, levels = unique(pca.data$Classification))),
  col = pca.data$group_color,
  pch = pca.data$sex_shape, cex = 3,
  xlab = axis_labels[1], ylab = axis_labels[2], zlab = axis_labels[3],
  bty = "b2", colkey = FALSE, phi = 20, theta = 30, d = 10, ticktype = "simple"
)
dev.off()

# FB
pdf('FB_3d_pca.pdf', width = 8, height = 8)
pca.data <- res_list$FB$pca_data
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=20, theta=30, d=10, ticktype="simple")
dev.off()


# KC
pca.data <- res_list$KC$pca_data
pdf('KC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=30, theta=60, d=10, ticktype="simple")
dev.off()


# KC_Channel
pca.data <- res_list$KC_Channel$pca_data
pdf('KC_Channel_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=40, theta=10, d=10, ticktype="simple")
dev.off()


# SGC
pca.data <- res_list$SGC$pca_data
pdf('SGC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=20, theta=50, d=10, ticktype="simple")
dev.off()

# HFC
pca.data <- res_list$HFC$pca_data
pdf('HFC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=50, theta=45, d=10, ticktype="simple")
dev.off()


# Sebocyte
pca.data <- res_list$Sebocyte$pca_data
pdf('Sebocyte_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=50, theta=10, d=10, ticktype="simple")
dev.off()


# MEC
pca.data <- res_list$MEC$pca_data
pdf('MEC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=50, theta=15, d=5, ticktype="simple")
dev.off()


# Pc-vSMC (Pc)
pca.data <- res_list$Pc$pca_data
pdf('Pc_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=40, theta=10, d=10, ticktype="simple")
dev.off()

# VEC
pca.data <- res_list$VEC$pca_data
pdf('VEC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=40, theta=10, d=10, ticktype="simple")
dev.off()

# LEC
pca.data <- res_list$LEC$pca_data
pdf('LEC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=10, theta=45, d=30, ticktype="simple")
dev.off()

# MEL
pca.data <- res_list$MEL$pca_data
pdf('MEL_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=10, theta=30, d=10, ticktype="simple")
dev.off()

# Schwann
pca.data <- res_list$Schwann$pca_data
pdf('Schwann_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=30, theta=60, d=10, ticktype="simple")
dev.off()

# T 
pca.data <- res_list$Tcell$pca_data
pdf('Tcell_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=20, theta=60, d=10, ticktype="simple")
dev.off()

# MonoMac
pca.data <- res_list$MonoMac$pca_data
pdf('MonoMac_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=20, theta=5, d=10, ticktype="simple")
dev.off()

# LC
pca.data <- res_list$LC$pca_data
pdf('LC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=10, theta=5, d=10, ticktype="simple")
dev.off()

# Mast
pca.data <- res_list$Mast$pca_data
pdf('Mast_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=pca.data$sex_shape, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=10, theta=20, d=10, ticktype="simple")
dev.off()


# sub-celltypes
# Pro-inflammatory.Fb
pca.data <- res_list$`Pro-inflammatory.Fb`$pca_data
pdf('Pro-inflammatory.Fb_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=20, theta=40, d=10, ticktype="simple")
dev.off()

# Basal.Ker.COL17A1
pca.data <- res_list$`Basal.Ker.COL17A1`$pca_data
pdf('Basal.Ker.COL17A1_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=30, theta=75, d=10, ticktype="simple")
dev.off()

# MLKC
pca.data <- res_list$MLKC$pca_data
pdf('MLKC_3d_pca.pdf', width=8, height=8)
scatter3D(x=pca.data$X, y=pca.data$Y, z=pca.data$Z,
          colvar=as.numeric(factor(pca.data$Classification, levels=unique(pca.data$Classification))),
          col=pca.data$group_color, pch=16, cex=3,
          xlab=axis_labels[1], ylab=axis_labels[2], zlab=axis_labels[3],
          bty="b2", colkey=FALSE, phi=30, theta=10, d=20, ticktype="simple")
dev.off()
