library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(multinichenetr)
library(extrafont)
pdf.options(family = "Arial")
library(ggvenn)
setwd("/data/yangyu/Project/skinAtlas/70samples/Result/11.CCC/2.nichenet/cell_type/")

### 1. Read in NicheNet V2 ligand-target prior model ####
lr_network = readRDS('/home/yangyu/Project/pipeline/bin/nicheneter/database/lr_network_human_21122021.rds')
lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix = readRDS('/home/yangyu/Project/pipeline/bin/nicheneter/database/ligand_tf_matrix_nsga2r_final.rds')
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

weighted_networks = readRDS('/home/yangyu/Project/pipeline/bin/nicheneter/database/weighted_networks_nsga2r_final.rds')
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

### 2.Prepare potential_ligands geneset_oi#### 
# celltype  target gene####
marker_path <- "/home/yangyu/Project/skinAtlas/70samples/Result/3.Conserved_markers/"
all_marker <- read.csv('/home/yangyu/Project/skinAtlas/70samples/Result/3.Conserved_markers/combined_conserved_markers.csv')
all_marker <- all_marker  %>% filter(minimump_p_val<0.05) %>% filter(max_pval<0.05)
all_marker <- all_marker %>%
  rowwise() %>%
  mutate(avg_fc = mean(c(PS_avg_log2FC, NT_avg_log2FC, NF_avg_log2FC, NP_avg_log2FC, NE_avg_log2FC), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(cluster_id) 
all_marker <- all_marker  %>% filter(avg_fc>0.5)
#combined_data <- do.call(rbind, data_list)
# potential ligands ####
all_receptors <- unique(lr_network$to) 
expressed_receptors <- intersect(all_receptors, unique(all_marker$gene))
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
length(potential_ligands)  # 946
## module background gene set ####
all_gene <- read.delim('../all_module/var.csv',sep = ',')
background_expressed_genes <- all_gene %>% pull(X) %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)
# 22028
###4. Run the NicheNet ligand activity analysis####
library(dplyr)
calculate_ligand_activities <- function(cluster_id,  all_marker, ligand_target_matrix, background_expressed_genes, potential_ligands, save_dir) {
  print(cluster_id)
  geneset_oi <- all_marker %>% 
    filter(cluster_id == !!cluster_id) %>% 
    pull(gene) %>% 
    unique()
  
  
  print(length(geneset_oi))
  
  geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
  
  
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  
  if (is.null(ligand_activities)) {
    return(data.frame(), data.frame()) 
  }
  
  
  ligand_activities <- ligand_activities %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected)))
  
  
  best_upstream_ligands <- ligand_activities %>% arrange(-aupr_corrected) %>% top_n(50, aupr_corrected) %>% 
    pull(test_ligand) %>% unique()
  
  
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 250) %>%
    bind_rows() %>%
    drop_na()
  
  
  write.csv(ligand_activities, file = file.path(save_dir, paste0(cluster_id, "_ligand_activities.csv")), row.names = FALSE,quote = F)
  write.csv(active_ligand_target_links_df, file = file.path(save_dir, paste0(cluster_id, "_active_ligand_target_links.csv")), row.names = FALSE,quote = F)
  
  return(list(ligand_activities = ligand_activities, active_ligand_target_links_df = active_ligand_target_links_df))
}
save_dir <- "./"

results <- all_marker %>%
  group_by(cluster_id) %>%
  group_split() %>% 
  walk(~ {
    
    cluster_id <- unique(.x$cluster_id)
    print(cluster_id)
    lig_activity <- calculate_ligand_activities(cluster_id, all_marker, ligand_target_matrix, background_expressed_genes, potential_ligands, save_dir)
    
    # 
    print(paste("Saved:", cluster_id))
  })

#  all celltype  ######
marker_path <- "/home/yangyu/Project/skinAtlas/70samples/Result/3.Conserved_markers/"
broad_marker <- read.csv('/home/yangyu/Project/skinAtlas/70samples/Result/3.Conserved_markers/whole_skin_conserved_markers.csv')
broad_marker <- broad_marker  %>% filter(minimump_p_val<0.05) %>% filter(max_pval<0.05)
broad_marker <- broad_marker %>%
  rowwise() %>%
  mutate(avg_fc = mean(c(PS_avg_log2FC, NT_avg_log2FC, NF_avg_log2FC, NP_avg_log2FC, NE_avg_log2FC), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(cluster_id) 
broad_marker <- broad_marker  %>% filter(avg_fc>1)
calculate_ligand_activities <- function(cluster_id,  broad_marker, ligand_target_matrix, background_expressed_genes, potential_ligands, save_dir) {
  
  
  geneset_oi <- broad_marker %>% 
    filter(cluster_id == !!cluster_id) %>% 
    pull(gene) %>% 
    unique()
  
  print(cluster_id)
  print(length(geneset_oi))
  
  geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
  
  
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  
  if (is.null(ligand_activities)) {
    return(data.frame(), data.frame())  
  }
  
  
  ligand_activities <- ligand_activities %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected)))
  
  
  best_upstream_ligands <- ligand_activities %>% arrange(-aupr_corrected) %>% top_n(50, aupr_corrected) %>% 
    pull(test_ligand) %>% unique()
  
 
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 250) %>%
    bind_rows() %>%
    drop_na()
  
  
  write.csv(ligand_activities, file = file.path(save_dir, paste0(cluster_id, "_ligand_activities.csv")), row.names = FALSE,quote = F)
  write.csv(active_ligand_target_links_df, file = file.path(save_dir, paste0(cluster_id, "_active_ligand_target_links.csv")), row.names = FALSE,quote = F)
  
  return(list(ligand_activities = ligand_activities, active_ligand_target_links_df = active_ligand_target_links_df))
}
save_dir <- "./"
results <- broad_marker %>%
  group_by(cluster_id) %>%
  group_split() %>%  
  walk(~ {

    cluster_id <- unique(.x$cluster_id)
    
    lig_activity <- calculate_ligand_activities(cluster_id, broad_marker, ligand_target_matrix, background_expressed_genes, potential_ligands, save_dir)
    
    print(paste("Saved:", cluster_id))
  })

