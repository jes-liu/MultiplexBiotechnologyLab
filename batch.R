rm(list=ls())

# library(dplyr)
# library(Seurat)
# # remotes::install_version(package = 'Seurat', version = package_version('4.0.4'))
# library(SeuratObject)
# library(patchwork)
# library(openxlsx)
# library(ggplot2)
# library(dplyr)
# 
# dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Batch_Final.xlsx'
# 
# 
# raw <- data.frame(read.xlsx(dir))
# 
# col_count = c(50, 50, 50, 50, 50, 50, 21, 50, 50, 44)
# #col_count = c(50, 50, 50, 50, 50, 50, 19)
# batch = integer()
# for (x in 1:length(col_count)) {
#   batch = append(batch, rep(x, col_count[x]))
# }
# 
# raw_fil = raw[,2:(319+1)]
# 
# #data_raw = data.matrix(raw_fil)
# data_raw = t(data.matrix(raw_fil))
# 
# data = CreateSeuratObject(data_raw)
# #data$Batch = batch
# 
# #boxplot(raw_fil)
# 
# #FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# 
# #data_norm <- NormalizeData(data)
# data_norm = data
# allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=319)
# allProteins <- rownames(allData)
# data_scale <- ScaleData(allData, features=allProteins)
# data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))
# data_nn <- FindNeighbors(data_pca, dims=1:15)
# 
# #data_cca = RunMultiCCA(data_scale, genes.use=)
# 
# 
# data_clus <- FindClusters(data_nn, resolution=.3)
# data_umap <- RunUMAP(data_clus, dims=1:15, spread=2, n.neighbors=50, min.dist=.001)
# DimPlot(data_umap, reduction='umap', label=TRUE)
# #DimPlot(data_umap, reduction='umap', label=TRUE, group.by="Batch")


###############################################
# # remotes::install_version(package = 'Seurat', version = package_version('3.2.0'))
# # use seurat 2 (RunMultiCCA) and seurat 3 (mmnCorrect)
# rm(list=ls())
# 
# library(Seurat)
# library(rliger)
# library(openxlsx)
# library(dplyr)
# 
# dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Batch_adjusted.xlsx'
# zero_dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Zero cells adjusted.xlsx'
# raw = data.frame(read.xlsx(dir))
# 
# #raw2 = t(data.matrix(raw))
# 
# batch.list = list("Batch 1" = raw[,2:51],  #50
#                   "Batch 2" = raw[,52:101],  #50
#                   "Batch 3" = raw[,102:151],  #50
#                   "Batch 4" = raw[,152:201],  #50
#                   "Batch 5" = raw[,202:251],  #50
#                   "Batch 6" = raw[,252:301],  #50
#                   "Batch 7" = raw[,302:320],  #19
#                   "Batch 8" = raw[,321:370],  #50
#                   "Batch 9" = raw[,371:420],  #50
#                   "Batch 10" = raw[,421:464],  #44
#                   "Batch 7_" = raw[,465:466]  #2
#                   )
# 
# liger = createLiger(batch.list)
# 
# #liger_norm = rliger::normalize(liger)
# liger@norm.data = liger@raw.data
# 
# 
# liger_genes = selectGenes(liger)
# 
# liger_scale = scaleNotCenter(liger_genes)
# 
# liger_als = optimizeALS(liger_scale, k=20)  #suggestK()
# # k is the inner dimension of factorization. higher k needed for datasets w/ more sub-structure
# 
# liger_quantile = quantile_norm(liger_als)  # quantile normalization
# 
# #quantileAlignSNF
# 
# liger_umap = runUMAP(liger_quantile, use.raw=T)
# p = plotByDatasetAndCluster(liger_umap, return.plots=T, pt.size=2)
# print(p[[1]])


##########################################
rm(list=ls())

library(limma)
library(Seurat)
library(openxlsx)
library(dplyr)
library(edgeR)

dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4_Final.xlsx'
zero_dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4 Zero cells adjusted.xlsx'

col_count = c(50, 50, 50, 50, 50, 50, 21, 50, 50, 44)
membrane_count = c(50, 50, 50, 50, 50, 50, 19)
intra_count = c(2, 50, 50, 44)

data = data.frame(read.xlsx(dir, sheet=1))
zero = data.frame(read.xlsx(zero_dir))
membrane_final = data[,2:(sum(membrane_count)+1)]
intra_final = data[,(sum(membrane_count)+2):(sum(membrane_count)+sum(intra_count)+1)]
zero_final = zero[,2:(sum(membrane_count)+1)]  # zero cells for membrane proteins only


batch_num = integer()
for (x in 1:(length(col_count))) {
  batch_num = append(batch_num, rep(x, col_count[x]))
}
#batch_num = append(batch_num, rep(7, 2))

#membrane_log = log2(membrane_final+1)
#intra_log = log2(intra_final+1)


membrane_DGE = DGEList(membrane_final)
membrane_norm = calcNormFactors(membrane_DGE)
membrane_cpm = cpm(membrane_norm, log=TRUE)
membrane_limma = removeBatchEffect(membrane_cpm, batch_num[1:sum(membrane_count)])
membrane_subtract = membrane_limma - rowMeans(membrane_limma)  # normalization step

intra_DGE = DGEList(intra_final)
intra_norm = calcNormFactors(intra_DGE)
intra_cpm = cpm(intra_norm, log=TRUE)
intra_limma = removeBatchEffect(intra_cpm, 
                    batch_num[(sum(membrane_count)+1):(sum(membrane_count)+sum(intra_count))])
intra_subtract = intra_limma - rowMeans(intra_limma)  # normalization step

zero_limma = removeBatchEffect(zero_final, batch_num[1:sum(membrane_count)])


membrane_limma_df = data.frame(membrane_subtract)
intra_limma_df = data.frame(intra_subtract)
zero_limma_df = data.frame(zero_limma)

#membrane_limma_df[membrane_limma_df < 0] = 0

# ncols = length(colnames(membrane_limma_df))
# for (x in 1:ncols) {
#   #zero_mean = mean(zero_limma_df[,x])
#   zero_std = sd(zero_limma_df[,x])
#   bg = log2(zero_std)  # zero_mean + 3*zero_std
#   membrane_limma_df[,x] = membrane_limma_df[,x] - bg
# }
# membrane_limma_df[membrane_limma_df < 0] = 0

# write.xlsx(membrane_limma_df,
#           'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/membrane_limma.xlsx',
#           overwrite=T)
# 
# 
# membrane_limma_df$Cell = c(1:nrow(data))
# intra_limma_df$Cell = c(1:nrow(data))
# 
# combined = merge(membrane_limma_df, intra_limma_df, by="Cell")
# 
# membrane_limma_df = membrane_limma_df[-length(membrane_limma_df)]
# intra_limma_df = intra_limma_df[-length(intra_limma_df)]
# combined = combined[-1]
# 
# list_df = list("membrane" = membrane_limma_df, "intracellular" = intra_limma_df)
# 
# write.xlsx(list_df,'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Batch_limma.xlsx', 
#            overwrite=TRUE)
# write.xlsx(combined,
#           'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Batch_limma_combined.xlsx', 
#            overwrite=TRUE)


##################### seurat after limma
#data_raw= data.matrix(membrane_limma_df)
data_raw = t(data.matrix(membrane_limma_df))

#row.names(data_raw) = row.names(membrane_limma_df)  # for combined

data = CreateSeuratObject(data_raw)
#data$Batch = batch_num[1:sum(membrane_count)]  # for membrane df
#data$Batch = batch_num[(sum(membrane_count)+1):(sum(intra_count)+sum(membrane_count))]  
# for intra df
#data$Batch = batch_num

#data_norm = NormalizeData(data)
data_norm = data
allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=120)
allProteins <- rownames(allData)
data_scale <- ScaleData(allData, features=allProteins)
data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))
data_nn <- FindNeighbors(data_pca, dims=1:15)
data_clus <- FindClusters(data_nn, resolution=.5)
data_umap <- RunUMAP(data_clus, dims=1:15, n.neighbors=50, spread=2, min.dist=.001)
DimPlot(data_umap, reduction='umap', label=TRUE)
#DimPlot(data_umap, reduction='umap', label=TRUE, group.by="Batch")

# membrane_limma_df$Cell = c(1:nrow(membrane_limma_df))
# membrane_limma_df$seurat_clusters = as.numeric(data_umap$seurat_clusters) + 1
# write.xlsx(membrane_limma_df,
#            'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/membrane_limma.xlsx',
#            overwrite=T)


markers <- FindAllMarkers(data_umap)
top <- markers %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
DoHeatmap(data_umap, features=unique(top$gene)) + NoLegend()

#boxplot(membrane_limma_df)

#######################################
#library(harmony)


