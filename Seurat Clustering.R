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
label <- raw$type  #ad - 1:721, wt - 722:1396

# read only the columns with expression data and type
#data_raw <- raw[,c(3,5:length(raw))]
data_raw <- raw[,5:length(raw)]

# changing dataframe into matrix
data_raw <- data.matrix(data_raw)

# transpose matrix
data_raw <- t(data_raw)

# initialize the seurat object
data <- CreateSeuratObject(counts=data_raw)
data$groups <- label
head(data, 5)

# plot of features
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC
data_qc <- subset(data, subset = nFeature_RNA>181 & nCount_RNA<600000)

# normalization
data_norm <- NormalizeData(data_qc)

# identification of highly variable features ----
# separating the AD and WT for normalization then variable features
data_AD <- subset(data_norm, subset = groups=="AD")
data_WT <- subset(data_norm, subset = groups=="WT")

# i,e. they are highly expressed in some cells, and lowly expressed in others
allData <- FindVariableFeatures(data_norm, selection.method="vst")
features <- VariableFeatures(allData)
data_featAD <- FindVariableFeatures(data_AD, selection.method="vst", nfeatures=30)
data_featWT <- FindVariableFeatures(data_WT, selection.method="vst", nfeatures=30)
top30AD <- VariableFeatures(data_featAD)
top30WT <- VariableFeatures(data_featWT)

# combining top features of both groups
top <- unique(c(top30AD, top30WT))

# visualization of top features
plot1 <- VariableFeaturePlot(data_featAD)
plot2 <- LabelPoints(plot = plot1, points = top30AD, repel = TRUE)
plot2

# scaling the data
# only add in features that will be used in PCA
data_scale <- ScaleData(allData)

# running linear dimension reduction
data_pca <- RunPCA(data_scale, features=features)  # features=top?

# visualization
# Heat map using pca components
DoHeatmap(data_pca, features=top, group.by = 'groups')


# clustering ----

# filter out the technical noise of the dataset from its dimensions of pca
# data_jackstraw <- JackStraw(data_pca)
# data_score <- ScoreJackStraw(data_jackstraw)
# JackStrawPlot(data_jackstraw)

# find neighbors and cluster
data_nn <- FindNeighbors(data_pca)
data_clus <- FindClusters(data_nn, resolution=.5)

# run umap
feat_test <- c("TrkA", "Arp2", "p.SREBP1", "BID", "PDCD4", "LEF1", "K.Ras", "CHRFAM7A",
               "MEK1", "APOE", "GSK.3a.ß", "Lamin.A.C", "COX2", "CDK5", "PAK1",
               "Hsp90a.ß", "pAPP", "HMGCS1", "pSHP2", "TREM2", "NLK", "PLP1",
               "GLUT4", "cleaved.Caspase.9", "ERAB", "AXIN1", "p.JNK1.JNK2",
               "Acetylcholinesterase", "pß.Catenin", "ATG5", "GNAI3",
               "HSP60", "Protein.APC", "UFD1", "Notch.3", "TSC1", "Src", "PDK.1",
               "pEGFR", "p.PLC.U.03B3.1", "GRB2", "X4EBP1", "Cyclin.D1.D2",
               "GAPDH", "TNFR1", "ß.actin", "HSF1", "Notch.1", "CD68", "pAkt",
               "p53", "JAK3", "p.INSR", "PAK3", "Rac1", "Akt", "Rheb", "LC3",
               "Active.Caspase.3", "GSK.3a.ß")
data_umap <- RunUMAP(data_clus, features=top)
DimPlot(data_umap, reduction='umap')
