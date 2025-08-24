library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsci)
library(extrafont)
library(ggvenn)
library(UpSetR)
pdf.options(family = "Arial")

setwd("/data/yangyu/Project/skinAtlas/70samples/Analysis/results/26.skin_diease/BCC/6.merge_epi/hfc/0.pseudobulk")
# HFC Normal DaMiRseq_adjust_countdata -----------------------------------------------
sampleinfo <- read_tsv("sample_info.txt")
# pseudobulk_rawCount.csv from 1.BCC_scRNA.ipynb (merge HFC from BCC and Normal(Face and Rest) samples)
countdata <- read.delim("pseudobulk_rawCount.csv",sep = ',') 
countdata <- countdata %>% dplyr::select(-NP15_NP,-NP2_NP,-PS1_PS,-PS2_PS,-PS12_PS,-PS4_PS) 
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
data_norm <- DaMiR.normalization(SE, minCounts = 0, fSample = 0.7, hyper = "no") 

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
pdf.options(family = "Arial")#, encoding = "UTF-8")


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
group_color = c('Rest'="#7C878EFF",
  'NF'='#a7cee2',
             'BCC_HFC'="#CC0C00FF",
                'BCC_cancer'="#46732EFF"
                )
pca.data$group_color <- group_color[as.character(pca.data$Classification)]
pdf('BCC_data_adjust_pca3d_a.pdf',h=5,w=5)
scatter3D(
  x = pca.data$X,
  y = pca.data$Y,
  z = pca.data$Z,
  colvar = as.numeric(factor(pca.data$Classification, levels = c('Rest','NF','BCC_HFC','BCC_cancer'))), 
  col = group_color, 
  #pch = pca.data$sex_shape, 
  pch = 16,
  cex = 2, 
  xlab = pca.var.per[1],
  ylab = pca.var.per[2],
  zlab = pca.var.per[3],
  bty = "b2", 
  colkey = F,
  phi = 10,theta=-35, 
  d = 10,
  ticktype = "simple" 
)
legend("topright", 
       legend = c('Rest','NF','BCC_HFC','BCC_cancer'), 
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



SOX11_norm <- data_adjust['SOX11',]  %>%t() %>% as.data.frame()
colnames(SOX11_norm) <- 'SOX11'
SOX11_norm <- SOX11_norm%>%rownames_to_column('sample_id_all')
meta_data <- sampleinfo 
SOX11_norm <- left_join(SOX11_norm, meta_data, by = c("sample_id_all" = "sample_id_all"))
SOX11_norm$dataset <- factor(SOX11_norm$dataset,levels=c('Rest','NF','BCC_HFC','BCC_cancer'))
library(ggpubr)
p <- ggboxplot(SOX11_norm, x = "dataset", y = "SOX11",color = 'dataset',add = "jitter") + 
  stat_compare_means(label = "p.signif",
                     comparisons = rev(list(c("Rest","NF"),c("Rest","BCC_HFC"), c("Rest","BCC_cancer"),
                                            c("NF","BCC_HFC"), c("NF","BCC_cancer"),
                                            c("BCC_HFC",'BCC_cancer'))),
                     method = "wilcox",
                     hide.ns = TRUE, 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.001), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  scale_color_manual(values = c('BCC_HFC'="#CC0C00FF",
                                'BCC_cancer'="#46732EFF",
                                'NF'='#a7cee2',
                                'Rest'="#7C878EFF"))+
  theme(
    legend.position = "none",
    panel.grid = element_blank()#,
    #axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave('SOX11_exp_boxplot.pdf',h=3,w=3)






