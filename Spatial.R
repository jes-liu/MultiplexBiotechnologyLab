library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)
library(ggplot2)

dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Intensities_75CellLogNorm.xlsx"
#dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Intensities_Filtered50.xlsx"

raw <- data.frame(read.xlsx(dir, sheet=1))
data_raw <- raw[,2:length(raw)]

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

data_nn <- FindNeighbors(data_pca, dims=1:9)
data_clus <- FindClusters(data_nn, resolution=.5)
data_tsne <- RunTSNE(data_clus, dims=1:9, check_duplicates = FALSE)

DimPlot(data_tsne, label=TRUE)
#DimPlot(Embeddings(data_pca)[1:100,1:3], reduction='pca', label=TRUE, group.by='markers')

FeaturePlot(data_tsne, features=allProteins, min.cutoff=.75)  # 1400x2200

if(FALSE) {
  # protein 1 - cluster 0
  # protein 2 - cluster 17
  # protein 3 - cluster 28
  # protein 4 - cluster 16
  # protein 5 - cluster 20
  # protein 6 - cluster 19
  # protein 7 - cluster 11
  # protein 8 - cluster 0
  # protein 9 - cluster 24
  # protein 10 - cluster 6
  # protein 11 - cluster 22
  # protein 12 - cluster 18
  # protein 13 - cluster 29
  # protein 14 - cluster 14
  # protein 15 - cluster 27
  # protein 16 - cluster 23
  # protein 17 - cluster 3
  # protein 18 - cluster 2+5
  # protein 21 - cluster 21
  # protein 22 - cluster 7 + parts of 1
  # protein 23 - cluster 13 + parts of 1
  # protein 24 - cluster 15
  # protein 25 - cluster 4
  # protein 26 - cluster 12
  # protein 27 - cluster 8
  # protein 28 - cluster 25
  # protein 29 - cluster 10
  # protein 30 - cluster 9
  # protein 31 - cluster 26
  }



coords <- data.frame("x" = data_tsne[["tsne"]]@cell.embeddings[,1],
                     "y" = data_tsne[["tsne"]]@cell.embeddings[,2],
                     "cluster" = data_tsne$seurat_clusters)

write.xlsx(coords, "TSNE_coords.xlsx")
