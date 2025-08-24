# =============================
# generates a dendrogram of cell types based on marker gene expression
# Author: YuYang
# =============================

library(Seurat)
library(phangorn)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggsci)
library(gridExtra)
library(colorspace)
library(viridis)
library(extrafont)
library(ggdendro)
library(dendextend)

pdf.options(family = "Arial")
setwd('/data/yangyu/Project/skinAtlas/70samples/Fig/Figure1')

load('../../Analysis/results/4.scVIintegrated/seu_top10gene.rda')  
metadata <- read.csv('seu.metadata.csv', row.names = 1)  

expr_matrix <- FetchData(seu, vars = rownames(seu))
expr_matrix$categorical <- metadata$celltype_Granular

avg_expr <- aggregate(. ~ categorical, data = expr_matrix, FUN = mean)
rownames(avg_expr) <- avg_expr$categorical
avg_expr <- avg_expr[, -1]

corr_matrix <- cor(t(avg_expr), method = 'pearson')
dist_matrix <- as.dist(1 - corr_matrix)
hc <- hclust(dist_matrix, method = 'complete')

save(expr_matrix, avg_expr, dist_matrix, hc, file = 'data_for_dendrogram_binary.rda')


p <- ggtree(hc, branch.length = "none", ladderize = FALSE) +
  geom_tiplab(show.legend = FALSE)

celltype_metadata <- metadata %>%
  dplyr::select(celltype_Granular, cell_type) %>%
  distinct() %>%
  dplyr::rename(name = celltype_Granular, cell_type = cell_type)

grp <- table(celltype_metadata) %>%
  apply(2, function(x) names(x[x != 0]))



celltype_colors <- c(
  'KC' = '#1f77b4', 'KC_Channel' = '#ff7f0e', 'SGC' = '#279e68', 'HFC' = '#d62728',
  'Sebocyte' = '#aa40fc', 'MEC' = '#8c564b', 'FB' = '#e377c2', 'Pc-vSMC' = '#b5bd61',
  'VEC' = '#17becf', 'LEC' = '#aec7e8', 'MEL' = '#ffbb78', 'Schwann' = '#98df8a',
  'Lymphocyte' = '#ff9896', 'Mac-DC' = '#c5b0d5', 'LC' = '#c49c94', 'Mast' = '#f7b6d2'
)


p <- p %<+% celltype_metadata +
  geom_tiplab(show.legend = FALSE) +
  scale_color_manual(values = celltype_colors) +
  theme_tree(plot.margin = margin(2, 2, 2, 2))

p <- groupOTU(p, grp, "cell_type") + aes_(color =~ cell_type)

p +geom_label(aes(label=node))
p <- p + geom_nodepoint(aes(subset = node==68,x=x),size=3,color='grey') +
  geom_nodepoint(aes(subset = node==70,x=x),size=3,color='grey') +
  geom_nodepoint(aes(subset = node==72,x=x),size=3,color='grey')+
  geom_nodepoint(aes(subset = node==73,x=x),size=3,color='grey')+
  geom_nodepoint(aes(subset = node==76,x=x),size=3,color='grey')+
  geom_nodepoint(aes(subset = node==71,x=x),size=3,color='grey')+
  geom_nodepoint(aes(subset = node==69,x=x),size=3,color='grey')+
  geom_nodepoint(aes(subset = node==77,x=x),size=3,color='grey')
p
p <- flip(p, 68,69)
p <- flip(p, 83,84)
p <- flip(p, 72,73)
p <- flip(p, 98,99)
p <- flip(p, 90,91)
p <- flip(p, 75,76)
p <- flip(p, 81,82)
p <- flip(p, 94,95)
p <- flip(p, 33,92)
p <- flip(p, 78,77)
p <- flip(p, 96,97)
p <- flip(p, 102,101)
p <- flip(p, 88,87)
p <- flip(p, 109,108)
p <- flip(p, 92,33)
p <- flip(p, 16,1)

taxa_order <- get_taxa_name(p)
write_tsv(as.data.frame(taxa_order), file = 'celltype_order.txt')

metadata$sample_group <- factor(metadata$sample_group, levels = c('NP','PS','NE','NT','NF'))

bar_data <- metadata %>%
  count(celltype_Granular, sample_group) %>%
  group_by(celltype_Granular) %>%
  mutate(Total = sum(n), Percentage = round(n / Total * 100, 2)) %>%
  ungroup() %>%
  rename(cell_type = celltype_Granular, Group = sample_group, Freq = n)

bar_data$cell_type <- factor(bar_data$cell_type, levels = taxa_order)

p2 <- p +
  geom_facet(panel = 'bar', data = bar_data, geom = geom_bar,
             mapping = aes(x = Percentage, fill = Group),
             orientation = "y", stat = "identity", position = "stack",
             size = 0, color = "transparent") +
  scale_fill_manual(values = c(
    'NE' = '#79a339', 'NF' = '#a7cee2', 'NT' = '#e1c7a7',
    'NP' = '#24467c', 'PS' = '#9999CC'
  )) +
  geom_tiplab(show.legend = TRUE) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), legend.position = 'right')



cell_count <- metadata %>%
  count(celltype_Granular) %>%
  mutate(logncells = log10(n), len = 1) %>%
  rename(ncells = n) %>%
  filter(celltype_Granular %in% taxa_order) %>%
  mutate(celltype_Granular = factor(celltype_Granular, levels = taxa_order))

p3 <- p2 +
  geom_facet(panel = 'num', data = cell_count,
             geom = geom_text,
             mapping = aes(label = ncells, x = 5),
             orientation = "y",
             color = 'black', size = 4, fontface = 'bold') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), legend.position = 'none')


ggsave('den_dendrogram_with_bar.pdf', p2, width = 20, height = 18)
ggsave('den_count_label.pdf', p3, width = 22, height = 18)


g <- ggplot(cell_count, aes(x = celltype_Granular, y = len, fill = ncells)) +
  geom_col() + coord_flip() +
  geom_text(aes(label = ncells, y = 0.5), color = "black", size = 5) +
  scale_fill_viridis_c(trans = 'log10', option = "magma",
                       limits = c(1, max(cell_count$ncells))) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        axis.title = element_blank())

ggsave('den_count_color.pdf', g, width = 2, height = 18)
