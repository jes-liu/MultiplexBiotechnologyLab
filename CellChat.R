library(patchwork)
library(CellChat)
library(openxlsx)
library(Seurat)
library(SeuratObject)


# Seurat
dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unified_minusBackG_Split.xlsx"

ad1_raw <- data.frame(read.xlsx(dir, sheet=1))
ad1_raw <- ad1_raw[,4:length(ad1_raw)]
ad1_raw <- t(ad1_raw)

ad2_raw <- data.frame(read.xlsx(dir, sheet=2))
ad2_raw <- ad2_raw[,4:length(ad2_raw)]
ad2_raw <- t(ad2_raw)

wt1_raw <- data.frame(read.xlsx(dir, sheet=3))
wt1_raw <- wt1_raw[,4:length(wt1_raw)]
wt1_raw <- t(wt1_raw)

wt2_raw <- data.frame(read.xlsx(dir, sheet=4))
wt2_raw <- wt2_raw[,4:length(wt2_raw)]
wt2_raw <- t(wt2_raw)

proteins <- rownames(ad1_raw)

AD1 <- CreateSeuratObject(counts=ad1_raw)
AD1 <- subset(AD1, subset = nCount_RNA<600000)
AD1 <- NormalizeData(AD1)
AD1 <- FindVariableFeatures(AD1)
AD1 <- ScaleData(AD1, features=proteins)
AD1 <- RunPCA(AD1, features=VariableFeatures(AD1))
AD1 <- FindNeighbors(AD1, dims=1:10)
AD1 <- FindClusters(AD1, resolution=.4)
AD1 <- RunUMAP(AD1, dims=1:10)
DimPlot(AD1, reduction='umap', label=T)

AD2 <- CreateSeuratObject(counts=ad2_raw)
AD2 <- subset(AD2, subset = nCount_RNA<600000)
AD2 <- NormalizeData(AD2)
AD2 <- FindVariableFeatures(AD2)
AD2 <- ScaleData(AD2, features=proteins)
AD2 <- RunPCA(AD2, features=VariableFeatures(AD2))
AD2 <- FindNeighbors(AD2, dims=1:10)
AD2 <- FindClusters(AD2, resolution=.6)
AD2 <- RunUMAP(AD2, dims=1:10)
DimPlot(AD2, reduction='umap', label=T)

WT1 <- CreateSeuratObject(counts=wt1_raw)
WT1 <- subset(WT1, subset = nCount_RNA<600000)
WT1 <- NormalizeData(WT1)
WT1 <- FindVariableFeatures(WT1)
WT1 <- ScaleData(WT1, features=proteins)
WT1 <- RunPCA(WT1, features=VariableFeatures(WT1))
WT1 <- FindNeighbors(WT1, dims=1:10)
WT1 <- FindClusters(WT1, resolution=.3)
WT1 <- RunUMAP(WT1, dims=1:10)
DimPlot(WT1, reduction='umap', label=T)

WT2 <- CreateSeuratObject(counts=wt2_raw)
WT2 <- subset(WT2, subset = nCount_RNA<600000)
WT2 <- NormalizeData(WT2)
WT2 <- FindVariableFeatures(WT2)
WT2 <- ScaleData(WT2, features=proteins)
WT2 <- RunPCA(WT2, features=VariableFeatures(WT2))
WT2 <- FindNeighbors(WT2, dims=1:10)
WT2 <- FindClusters(WT2, resolution=.5)
WT2 <- RunUMAP(WT2, dims=1:10)
DimPlot(WT2, reduction='umap', label=T)


# CellChat
ad1.input <- GetAssayData(data_umap, assay="RNA", slot="data")
labels <- data_umap$Sample
meta <- data.frame(group=labels, row.names=names(labels))


