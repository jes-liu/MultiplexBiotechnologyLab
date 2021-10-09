rm(list=ls())

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)
library(ggplot2)

dir = "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Intensities_minBG.xlsx"


raw <- data.frame(read.xlsx(dir))
data_raw <- raw[,1:length(raw)]

#boxplot(data_raw)


data_raw <- t(data.matrix(data_raw))

data <- CreateSeuratObject(data_raw)

# ----

data_norm <- data
#data_norm <- NormalizeData(data)
allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=29)
allProteins <- rownames(allData)
data_scale <- ScaleData(allData, features=allProteins)
data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))

# data_jackstraw <- JackStraw(data_pca)
# data_score <- ScoreJackStraw(data_jackstraw, dims=1:15)
# ElbowPlot(data_score)

data_nn <- FindNeighbors(data_pca, dims=1:28)
data_clus <- FindClusters(data_nn)
data_umap <- RunUMAP(data_clus, dims=1:28)
# min.dist, spread, n.neighbors

DimPlot(data_umap, label=TRUE)

FeaturePlot(data_umap, features=allProteins)  # 1400x2200
# min.cutoff=


# # get hist of num of values per cell in cluster 0
# # get hist of the value distribution in cluster 0
# index = which(data_tsne$seurat_clusters == 0)
# data_re = data.frame(t(data_raw))
# subset = data_re[index,]
# nrow = nrow(subset)
# ncol = ncol(subset)
# 
# array = vector()
# values = vector()
# for (x in 1:nrow) {
#   num = 0
#   for (y in 1:ncol) {
#     if (subset[x, y] > 0) {
#       num = num + 1
#       values = append(values, subset[x, y])
#     }
#   }
#   array = append(array, num)
# }
# hist(array)
# hist(values)
# 
# # get hist of num of values per cell of all else
# # get hist of the value distribution of all else
# not_index = which(data_tsne$seurat_clusters != 0)
# not_subset = data_re[not_index,]
# not_nrow = nrow(not_subset)
# 
# not_array = vector()
# not_values = vector()
# for (x in 1:not_nrow) {
#   num = 0
#   for (y in 1:ncol) {
#     if (not_subset[x, y] > 0) {
#       num = num + 1
#       not_values = append(not_values, not_subset[x, y])
#     }
#   }
#   not_array = append(not_array, num)
# }
# hist(not_array)
# hist(not_values)


coords <- data.frame("cell" = 1:length(data_umap$orig.ident),
                     "seurat_clusters" = data_umap$seurat_clusters)

write.xlsx(coords, "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/UMAP_clusters.xlsx",
           overwrite=TRUE)
