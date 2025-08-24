library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggsci)
library(extrafont)
pdf.options(family = "Arial")


# Mesenchymal.Ker Normal DaMiRseq_adjust_countdata -----------------------------------------------
sampleinfo <- read_tsv("sample_info.txt")
# pseudobulk raw count matrix from 5.EMPD_scRNA.ipynb (merge Mesenchymal.Ker from EMPD and Normal(PS and Rest))
countdata <- read.delim("pseudobulk_rawCount.csv",sep = ',') 
countdata2 <- countdata%>%column_to_rownames("X")%>%as.matrix()

dim(countdata2) 
keep <- rowSums(countdata2) > 200
countdata3 <- countdata2[keep,]

dim(countdata3)

library(DaMiRseq)
order_indices <- match(colnames(countdata3), sampleinfo$sample_id_all)
class1 <- sampleinfo[order_indices, ]
class3 <- class1%>%column_to_rownames("sample_id_all")%>%dplyr::rename(class = dataset) 
class3
SE<-DaMiR.makeSE(countdata3,class3)

SE
data_norm <- DaMiR.normalization(SE, minCounts = 0, fSample = 0.7, hyper = "no") #修改fSample = 0.7

data_filt <- DaMiR.sampleFilt(data_norm, th.corr = 0.5)

dim(data_filt) #

sv<- DaMiR.SV(data_filt) 
data_filt$class<-as.factor(data_filt$class)


DaMiR.corrplot(sv, colData(data_filt),sig.level = 0.01)

# save data_norm matrix
write.table(assay(data_norm), file="DaMiRseq_filtered/vst_normalized_expression.txt", row.names = TRUE, col.names = TRUE, sep = '\t')
data_adjust <-DaMiR.SVadjust(data_filt, sv) 
write.table(assay(data_adjust), file="DaMiRseq_adjust/SV_adjusted_expression.txt", row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)

assay(data_adjust[c(1:5),c(1:5)])
dev.off()


library(bioplotr)
pdf('DaMiRseq_adjust/adjust_umap.pdf',w=5,h=4)
plot_umap(assay(data_adjust), group=data_adjust$class, size=3, label = TRUE,pal_group='d3')
plot_umap(assay(data_adjust), group=data_adjust$class, size=3,pal_group='d3')
dev.off()

pdf('DaMiRseq_adjust/adjust_pca.pdf',w=5,h=4)
plot_pca(assay(data_adjust), group=data_adjust$class, size=3, label = TRUE,pal_group='d3')
plot_pca(assay(data_adjust), group=data_adjust$class, size=3,pal_group='d3')
plot_pca(assay(data_adjust), group=data_adjust$class, size=3, pcs = c(3L,4L), label = TRUE,pal_group='d3')
dev.off()



pdf('DaMiRseq_filtered/filter_umap.pdf',w=5,h=4)
plot_umap(assay(data_filt), group=data_filt$class, size=3, label = TRUE,pal_group='d3')
plot_umap(assay(data_filt), group=data_filt$class, size=3,pal_group='d3')
dev.off()

pdf('DaMiRseq_filtered/filter_pca.pdf',w=5,h=4)
plot_pca(assay(data_filt), group=data_filt$class, size=3, label = TRUE,pal_group='d3')
plot_pca(assay(data_filt), group=data_filt$class, size=3,pal_group='d3')
plot_pca(assay(data_filt), group=data_filt$class, size=3, pcs = c(3L,4L), label = TRUE,pal_group='d3')
dev.off()
library(plotly)
library(ggplot2)
library(RColorBrewer)
pca <- prcomp(data.frame(t(assay(data_adjust))), scale=F)
pcs <- data.frame(pca$x, Group=data_adjust$class, X=pca$x[,1],Y=pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per <- paste(colnames(pcs),"(",paste(as.character(pca.var.per),"%",")",sep=""))
pca.data <- data.frame(Sample=rownames(pca$x),Classification=data_adjust$class,X=pca$x[,1],Y=pca$x[,2],Z=pca$x[,3])
fig <- plot_ly(pca.data, x = ~X, y = ~Y, z = ~Z, color = ~Classification,colors=pal_simpsons(palette = 'springfield',alpha = 1)(10))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = pca.var.per[1]),yaxis = list(title = pca.var.per[2]), zaxis = list(title = pca.var.per[3])))
fig


pca.data$Classification <- as.character(pca.data$Classification)

group_color = c('Rest_Mesenchymal'="#7C878EFF",
                'PS_Mesenchymal'='#9999CC',
                'EMPD_Mesenchymal'="#FFCD00FF",
                'EMPD_Paget'="#197EC0FF")
pca.data$group_color <- group_color[as.character(pca.data$Classification)]
pca.data$Classification <- factor(pca.data$Classification, levels = c('Rest_Mesenchymal','PS_Mesenchymal','EMPD_Mesenchymal','EMPD_Paget'))

library(plot3D)
pdf('EMPD_data_adjust_pca3d_a.pdf',h=5,w=5)
scatter3D(
  x = pca.data$X,
  y = pca.data$Y,
  z = pca.data$Z,
  colvar = as.numeric(factor(pca.data$Classification, levels = c('Rest_Mesenchymal','PS_Mesenchymal','EMPD_Mesenchymal','EMPD_Paget'))), 
  col = group_color, 
  pch = 16,
  cex = 2, 
  xlab = pca.var.per[1],
  ylab = pca.var.per[2],
  zlab = pca.var.per[3],
  bty = "b2", 
  colkey = F,
  phi = 5,theta=-15, 
  d = 50,
  ticktype = "simple" 
)
legend("topright", 
       legend = c('Rest_Mesenchymal','PS_Mesenchymal','EMPD_Mesenchymal','EMPD_Paget'), # 图例项
       fill = group_color, 
       title = "Classification", 
       cex = 1.2) 
text3D(x = pca.data$X,
       y = pca.data$Y,
       z = pca.data$Z,labels = pca.data$Sample,add = TRUE,
       col = 'black',
       cex = 0.6,
       pos = 3)
dev.off()

# PITX1 expression boxplot -------------------------------------------------
PITX1_norm <- data_adjust['PITX1',]  %>%t()%>% as.data.frame() 
colnames(PITX1_norm) <- 'PITX1'
PITX1_norm <- PITX1_norm%>%rownames_to_column('sample_id_all')
meta_data <- meta_data %>%rownames_to_column('sample_id_all')
PITX1_norm <- left_join(PITX1_norm, meta_data, by = c("sample_id_all" = "sample_id_all"))
PITX1_norm$dataset <- factor(PITX1_norm$dataset,levels=c('Rest_Mesenchymal','PS_Mesenchymal','EMPD_Mesenchymal','EMPD_Paget'))
library(ggpubr)
p <- ggboxplot(PITX1_norm, x = "dataset", y = "PITX1",color = 'dataset',add = "jitter") + 
  stat_compare_means(label = "p.signif",
                     comparisons = rev(list(c("Rest_Mesenchymal","PS_Mesenchymal"),c("Rest_Mesenchymal","EMPD_Mesenchymal"), c("Rest_Mesenchymal","EMPD_Paget"),
                                            c("PS_Mesenchymal","EMPD_Mesenchymal"), c("PS_Mesenchymal","EMPD_Paget"),
                                            c("EMPD_Mesenchymal",'EMPD_Paget'))),
                     method = "wilcox",
                     hide.ns = TRUE, 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.001), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  scale_color_manual(values = c('EMPD_Mesenchymal'="#FFCD00FF",
                                'PS_Mesenchymal'='#9999CC',
                                'Rest_Mesenchymal'="#7C878EFF",
                                'EMPD_Paget'='#197EC0FF'))+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave('PITX1_exp_boxplot.pdf',h=3,w=3)



