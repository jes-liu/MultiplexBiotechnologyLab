rm(list=ls())

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)

print("Reading Data")
dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unifiedB.xlsx"


# preprocessing ----
# read the data into a dataframe
raw <- data.frame(read.xlsx(dir))

# get data groups
type <- raw$type  #ad - 1:721, wt - 722:1396
label <- raw$Label

# read only the columns with expression data and type
#data_raw <- raw[,c(3,5:length(raw))]
data_raw <- raw[,5:length(raw)]

# changing dataframe into matrix
data_raw <- data.matrix(data_raw)

# transpose matrix
data_raw <- t(data_raw)

# initialize the seurat object
data <- CreateSeuratObject(counts=data_raw)
print("Data Read")
data$Condition <- type
data$Sample <- label
head(data, 5)

# plot of features
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC
print("Preprocessing Data")
data_qc <- subset(data, subset = nFeature_RNA>181 & nCount_RNA<600000)

# normalization
data_norm <- NormalizeData(data_qc)
print("Data Preprocessed")

# identification of highly variable features ----
# i,e. they are highly expressed in some cells, and lowly expressed in others
print("Finding Variable Features")
allData <- FindVariableFeatures(data_norm, selection.method="vst")

# scaling the data
# only add in features that will be used in PCA
allProteins <- rownames(allData)
print("Scaling Data")
data_scale <- ScaleData(allData, features=allProteins)

# running linear dimension reduction
print("Running PCA")
data_pca <- RunPCA(data_scale, features=VariableFeatures(data_scale))
DimPlot(data_pca)

# clustering ----

# filter out the technical noise of the dataset from its dimensions of pca
print("Filtering Technical Noise")
data_jackstraw <- JackStraw(data_pca)
data_score <- ScoreJackStraw(data_jackstraw, dims=1:15)
JackStrawPlot(data_score)
ElbowPlot(data_score)
print("Dimensions Found")

# find neighbors and cluster
print("Starting Cluster")
data_nn <- FindNeighbors(data_score, dims=1:10)
data_clus <- FindClusters(data_nn, resolution=.5)

# run umap
data_umap <- RunUMAP(data_clus, dims=1:10)
print("Cluster Finished")
DimPlot(data_umap, reduction='umap', label=TRUE,
        group.by=c("seurat_clusters", "Condition", "Sample"), ncol=3)


# heatmap ----

# finding markers
print("Finding Markers")
AD_markers <- FindMarkers(data_umap, ident.1="AD", only.pos=TRUE,
                          group.by="Condition", logfc.threshold=.1)
WT_markers <- FindMarkers(data_umap, ident.1="WT", only.pos=TRUE,
                          group.by="Condition", logfc.threshold=.1)
topAD <- AD_markers %>% top_n(n=30, wt=avg_log2FC)
topWT <- WT_markers %>% top_n(n=30, wt=avg_log2FC)
total_top <- c(rownames(topAD), rownames(topWT))
print("Markers Found")
DoHeatmap(data_umap, features=total_top, group.by="Condition") + NoLegend()

