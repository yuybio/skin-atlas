library(phangorn)
library(ggtree)
library(dendextend)
library(ggdendro)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(colorspace)
library(viridis)
library(Seurat)
library(dplyr)
library(scales)
library(ggsci)
library(gridExtra)
library(extrafont)
pdf.options(family = "Arial")

setwd('/data/yangyu/Project/skinAtlas/70samples/Fig/Figure1')
load('../../Analysis/results/4.scVIintegrated/seu_top10gene.rda')
metadata <- read.csv('seu.metadata.csv',row.names = 1)

rep_df = FetchData(seu,vars  = rownames(seu)) 

categorical = metadata$celltype_Granular
rep_df <- data.frame(rep_df, categorical = categorical) 
mean_df <- aggregate(. ~ categorical, data = rep_df, FUN = mean)
rownames(mean_df) <- mean_df$categorical
mean_df <- mean_df[,-1]

corr_matrix <- cor(t(mean_df), method = 'pearson')
corr_condensed <- as.dist(1 - corr_matrix)
hc <- hclust(corr_condensed, method = 'complete')
save(rep_df,categorical,mean_df,corr_condensed,hc,file = 'data_for_dendrogram_binary.rda')
plot(hc)
dhc <- as.dendrogram(hc)
plot(dhc)
p <- ggtree(hc, branch.length = "none",ladderize=F) +
  geom_tiplab(show.legend=FALSE)

dat <- data.frame(name = metadata$celltype_Granular,cls = metadata$cell_type)

celltype_metadata <- dat[!duplicated(dat$name),]
rownames(celltype_metadata) <- NULL


grp <- apply(table(dat), 2, function(x) names(x[x != 0])) 

p <- groupOTU(p, grp, "cell_type") + aes_(color =~ cell_type)

# extract the most recent common ancestor
noids <- lapply(grp, function(x) unlist(lapply(x, function(i) ggtree::nodeid(p, i))))
roots <- unlist(lapply(noids, function(x) ggtree::MRCA(p, x))) 
rangeX <- max(p$data$x, na.rm=TRUE) - min(p$data$x, na.rm=TRUE)


pdata <- data.frame(name = p$data$label, color2 = p$data$cell_type)
pdata <- pdata[!is.na(pdata$name), ]
cluster_color <- unique(pdata$color2)
n_color <- length(levels(cluster_color)) - length(cluster_color)


bar_data <- data.frame(table(metadata$celltype_Granular,metadata$sample_group))

metadata$sample_group <- factor(metadata$sample_group,levels=c('NP','PS','NE','NT','NF'))

bar_data$Total <- apply(bar_data,1,function(x)sum(bar_data[bar_data$Var1 == x[1],3]))
bar_data <- bar_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
colnames(bar_data) <- c('cell_type','Group','Freq','Total', 'Percentage')
bar_data <- na.omit(bar_data)


# Dendrogram ####
p <- p %<+%  celltype_metadata +
  geom_tiplab(show.legend=FALSE) + 
  scale_color_manual(values = c(
    'KC'= '#1f77b4',
    'KC_Channel'= '#ff7f0e',
    'SGC'= '#279e68',
    'HFC'= '#d62728',
    'Sebocyte'= '#aa40fc',
    'MEC'= '#8c564b',
    'FB'= '#e377c2',
    'Pc-vSMC'= '#b5bd61',
    'VEC'= '#17becf',
    'LEC'= '#aec7e8',
    'MEL'= '#ffbb78',
    'Schwann'= '#98df8a',
    'Lymphocyte'= '#ff9896',
    'Mac-DC'= '#c5b0d5',
    'LC'= '#c49c94',
    'Mast'= '#f7b6d2'
  )) +
  theme_tree(plot.margin=margin(2,2,2,2)) 
p
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

get_taxa_name(p)
aa <- get_taxa_name(p)
aa = as.data.frame(aa)
write_tsv(aa,file = 'celltype_order.txt')





p <- ggtree(hc, branch.length = "none",ladderize=F)+geom_tiplab()+theme(legend.position = 'none')+hexpand(.06)
p <- groupOTU(p, grp, "cell_type") + aes_(color =~ cell_type)
p
p <- p %<+%  celltype_metadata +
  #geom_tiplab(aes(label=celltype_Granular), size=3, hjust=.5, color='black') +
  geom_tiplab(show.legend=FALSE) + 
  scale_color_manual(values = c('KC'= '#1f77b4',
                                'KC_Channel'= '#ff7f0e',
                                'SGC'= '#279e68',
                                'HFC'= '#d62728',
                                'Sebocyte'= '#aa40fc',
                                'MEC'= '#8c564b',
                                'FB'= '#e377c2',
                                'Pc-vSMC'= '#b5bd61',
                                'VEC'= '#17becf',
                                'LEC'= '#aec7e8',
                                'MEL'= '#ffbb78',
                                'Schwann'= '#98df8a',
                                'Lymphocyte'= '#ff9896',
                                'Mac-DC'= '#c5b0d5',
                                'LC'= '#c49c94',
                                'Mast'= '#f7b6d2')) +
  theme_tree(plot.margin=margin(2,2,2,2)) 
p
p + geom_label(aes(label=node))
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

p+geom_label(aes(label=node))
p <- p+ geom_point2(aes(subset=node==68), color='steelblue', alpha = 0.6,size=7) +
  geom_point2(aes(subset=node==70), color='firebrick',alpha = 0.6,size=7)+
  geom_point2(aes(subset=node==71), color='darkgreen', alpha = 0.6,size=7)
p

df <- p$data
nodeCo <- intersect(df %>% dplyr::filter(is.na(df$x)) %>% dplyr::select(parent, 
                                                                        node) %>% unlist(), df %>% dplyr::filter(!is.na(x)) %>% 
                      dplyr::select(parent, node) %>% unlist())

labCo <- df %>% dplyr::filter(node %in% nodeCo) %>% 
  dplyr::select(label) %>% unlist()
data = bar_data
selCo <- intersect(labCo, rownames(data))
isSel <- df$label %in% selCo
df <- df[df$isTip | isSel, ]
offset = 0 
start <- max(df$x, na.rm = TRUE) + offset
dd <- as.data.frame(data)
i <- order(df$y)
i <- i[!is.na(df$y[i])]
lab <- df$label[i]

sorted_dd <- dd[order(factor(dd$cell_type, levels = lab)), ] #排序
sorted_dd$cell_type <- factor(sorted_dd$cell_type,levels = lab)
sorted_dd$cell <- sorted_dd$cell_type
sorted_dd$cell <- as.factor(sorted_dd$cell)
sorted_dd$Group <- as.factor(sorted_dd$Group)


noids <- lapply(grp, function(x) unlist(lapply(x, function(i) ggtree::nodeid(p, i))))
roots <- unlist(lapply(noids, function(x) ggtree::MRCA(p, x)))

p2 <- p +
  geom_facet(panel = 'bar',data = sorted_dd, geom = geom_bar, mapping = aes(x = Percentage, fill = Group),orientation="y",
             stat = "identity",position = "stack",size=0,color = "transparent") +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text()) + scale_fill_manual(values = c('NE'='#79a339', 
                                                                     'NF'='#a7cee2', 
                                                                     'NT'='#e1c7a7', 
                                                                     'NP'='#24467c',
                                                                     'PS'='#9999CC'))+
  geom_tiplab(show.legend=TRUE)

p2

facet_widths(p2,widths = c(0.5,0.5))


cell_count <- data.frame(table(metadata$celltype_Granular))
colnames(cell_count) <- c('celltype_Granular','ncells')
cell_count$logncells <- log10(cell_count$ncells)
cell_count <- cell_count[order(factor(cell_count$celltype_Granular, levels = lab)), ] 
cell_count$celltype_Granular <- factor(cell_count$celltype_Granular,levels = lab)
cell_count$ncells <- as.numeric(cell_count$ncells)
cell_count$logncells <- as.numeric(cell_count$logncells)
cell_count$len <- rep(1,66)
cell_count <- na.omit(cell_count)
g <- ggplot(cell_count, aes(x = celltype_Granular, y = len, fill = ncells))+
  geom_col() +
  coord_flip() +
  theme_minimal()  + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())  +
  geom_text(aes(label = ncells, y = 0.5), color = "black", size =5, fontweight = "bold") +
  scale_fill_viridis_c(trans = 'log10', option = "magma", limits = c(1, max(cell_count$ncells))) 
  
g
ggsave('den_count_color.pdf',w=2,h=18)

p3 <- p2 + geom_facet(panel = 'num',data = cell_count, geom = geom_text, mapping = aes(label=ncells,x=5),
                      orientation="y",color = 'black', size = 4,fontface = 'bold',#width = 1,
                      heigth = 1) +theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')  

facet_widths(p3,widths = c(0.4,0.1,0.15))
ggsave('den_count.pdf',w=22,h=18)
p4 <- p3 + geom_facet(panel = 'num',data = cell_count, geom = geom_col, mapping = aes(x = len,fill = ncells),
                      orientation="y") + theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())  +
  scale_fill_viridis_c(trans = 'log10', option = "magma", limits = c(1, max(cell_count$ncells))) 


p4 <- p2 + geom_facet(panel = 'num',data = cell_count, geom = geom_text, mapping = aes(label=ncells,x=5),
                      orientation="y",color = 'black', size = 4,# fontface = 'bold',width = 1,
                      heigth = 1) +theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #text = element_text(family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')  

facet_widths(p4,widths = c(0.4,0.1,0.15))
ggsave('den_count.pdf',w=22,h=18)
