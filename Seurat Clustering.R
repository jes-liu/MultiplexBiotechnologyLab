rm(list=ls())

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)

dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unifiedB.xlsx"


# preprocessing ----
# read the data into a dataframe
raw <- data.frame(read.xlsx(dir))

# get data labels
label <- raw$type

# read only the columns with expression data
data_raw <- raw[,5:length(raw)]

# changing dataframe into matrix
data_raw <- data.matrix(data_raw)

# transpose matrix
data_raw <- t(data_raw)

# initialize the seurat object
data <- CreateSeuratObject(counts=data_raw)
#data2 <- GetAssayData(data, slot="data")
head(data, 5)

# plot of features
VlnPlot(data2, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC
data_qc <- subset(data, subset = nFeature_RNA>181 & nCount_RNA<600000)

# normalization
data_norm <- NormalizeData(data_qc)

# identification of highly variable features
# i,e. they are highly expressed in some cells, and lowly expressed in others
data_id <- FindVariableFeatures(data_norm)
top10 <- head(VariableFeatures(data_id), 10)
plot1 <- VariableFeaturePlot(data_id)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# scaling the data
row.names <- rownames(data_id)  # only add in features that will be used in PCA
data_scale <- ScaleData(data_id, features=row.names)


# performing ----
# running linear dimension reduction
data_pca <- RunPCA(data_scale)  # features = VariableFeatures(object = subset)

# visualization
VizDimLoadings(data_pca, dims = 3, reduction = "pca")  # dims = 1:5
# graphs output of dim reduc technique
# each point is a cell
DimPlot(data_pca)
# Heat map using pca components
DimHeatmap(data_pca, dims=1, cells=250, balanced=TRUE)
