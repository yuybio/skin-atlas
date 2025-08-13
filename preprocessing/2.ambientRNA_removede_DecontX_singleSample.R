library(celda)
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(zellkonverter)
library(ggplot2)
library(anndata)

projectdir = '/data/Project/skinAtlas/Analysis/results/2.ambiantRNAremoved_DecontX/h5adFiles/'

sampleinfo <- read.table('/data/Project/skinAtlas/70samples/sample_info.txt',sep='\t',header = T )
samples <- sampleinfo$sample_id

### import cellranger files from different data sets###
for (i in seq_along(samples)){
  assign(paste0("scs_data", samples[i]), Read10X(data.dir = paste0("/data/Project/skinAtlas/cellranger/", samples[i], "/outs/filtered_feature_bc_matrix")))
}
# Create a SingleCellExperiment object and run decontX ###
for (i in seq_along(samples)){
  assign(paste0("sceraw", samples[i]),SingleCellExperiment(list(counts = eval(parse(text = paste0("scs_data", samples[i]))))))
}

for (i in seq_along(samples)){
  assign(paste0("sce", samples[i]),decontX(eval(parse(text = paste0("sceraw", samples[i])))))
}

for  (i in seq_along(samples)){
  scename = paste0("sce", samples[i])
  save(list=scename,
       file=paste('/data/Project/skinAtlas/Analysis/results/2.ambiantRNAremoved_DecontX/sceFiles/',samples[i],'.rda',sep=''))
}

## creat AnnData from a SCE with decontX result ###
for (i in seq_along(samples)){
  adname=paste0("adata", i)
  assign(adname,SCE2AnnData(eval(parse(text = paste0("sce", samples[i]))),X_name='decontXcounts'))
  write_h5ad(eval(as.name(adname)),filename=paste(projectdir,samples[i],'.h5ad',sep=''))
}
