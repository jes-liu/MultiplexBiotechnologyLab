rm(list=ls())

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)
library(ggplot2)

print("Reading Data")
dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unified_minusBackG.xlsx"


# preprocessing ----
# read the data into a dataframe
raw <- data.frame(read.xlsx(dir))

# get data groups
type <- raw$type
label <- raw$Label

# read only the columns with expression data
data_raw <- raw[,4:length(raw)]

# remove housekeeping genes
data_raw <- subset(data_raw, select= -c(GAPDH, Lamin.A.C, Lamin.B1, X.U.0251..Adaptin))

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

# beta-actin, gapdh, laminb1, lamin a/c

# plot of features
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC
print("Preprocessing Data")
data_qc <- subset(data, subset = nCount_RNA<600000)

# normalization
data_norm <- NormalizeData(data_qc)
print("Data Preprocessed")

# identification of highly variable features ----
# i,e. they are highly expressed in some cells, and lowly expressed in others
# print("Finding Variable Features")
allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=100)

# scaling the data
# only add in features that will be used in PCA
allProteins <- rownames(allData)
print("Scaling Data")
data_scale <- ScaleData(allData, features=allProteins)

# running linear dimension reduction
print("Running PCA")
data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))

# clustering ----

# filter out the technical noise of the dataset from its dimensions of pca
# print("Filtering Technical Noise")
# data_jackstraw <- JackStraw(data_pca)
# data_score <- ScoreJackStraw(data_jackstraw, dims=1:15)
# ElbowPlot(data_score)
# print("Dimensions Found")

# find neighbors and cluster
print("Starting Cluster")
data_nn <- FindNeighbors(data_pca, dims=1:10)
data_clus <- FindClusters(data_nn, resolution=.24)

# run umap
data_umap <- RunUMAP(data_clus, dims=1:10)
print("Cluster Finished")
DimPlot(data_umap, reduction='umap', label=TRUE,
        group.by=c("seurat_clusters", "Condition", "Sample"), ncol=3)
DimPlot(data_umap, reduction='umap', label=TRUE, split.by="Condition")


# heatmap ----

# finding markers
print("Finding Markers")

# AD-WT heatmap
AD_markers <- FindMarkers(data_umap, ident.1="AD", only.pos=TRUE,
                          group.by="Condition", logfc.threshold=.1)
AD_markers
WT_markers <- FindMarkers(data_umap, ident.1="WT", only.pos=TRUE,
                          group.by="Condition", logfc.threshold=.1)
topAD <- AD_markers %>% top_n(n=30, wt=avg_log2FC)
topWT <- WT_markers %>% top_n(n=30, wt=avg_log2FC)
total_top <- c(rownames(topAD), rownames(topWT))
DoHeatmap(data_umap, features=total_top, group.by="Condition") + NoLegend()

# cluster heatmap
cluster_markers <- FindAllMarkers(data_umap, only.pos=TRUE, logfc.threshold=.1,
                                  min.pct=.1)
top <- cluster_markers %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
DoHeatmap(data_umap, features=top$gene) + NoLegend()


# dotplot ----
dot_markers <- FindAllMarkers(data_umap, test.use = "wilcox", 
                              logfc.threshold=.1, only.pos=TRUE)
#dot_markers_fil <- subset(dot_markers, subset= pct.2<1)
top_dot <- dot_markers %>% group_by(cluster) %>% top_n(10, wt=pct.2)
DotPlot(data_umap, features=unique(top_dot$gene)) + RotatedAxis() + coord_flip()
print("Plots Created")

FeaturePlot(data_umap, features=unique(head(top_dot$gene, 9)))
VlnPlot(data_umap, features=head(top_dot$gene, 9), sort=T, log=T, pt.size=0)




# Seurat - separate AD and WT ----

dir<-"C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unified_minusBackG_SplitCondition.xlsx"

ad_raw <- data.frame(read.xlsx(dir, sheet=1))
wt_raw <- data.frame(read.xlsx(dir, sheet=2))

ad_raw <- ad_raw[,4:length(ad_raw)]
wt_raw <- wt_raw[,4:length(wt_raw)]

# remove housekeeping genes
ad_filt <- subset(ad_raw, select= -c(GAPDH, Lamin.A.C, Lamin.B1, X.U.0251..Adaptin))
wt_filt <- subset(wt_raw, select= -c(GAPDH, Lamin.A.C, Lamin.B1, X.U.0251..Adaptin))

# preprocessing
ad_filt <- t(ad_filt)
wt_filt <- t(wt_filt)

ad_data <- CreateSeuratObject(counts=ad_filt)
wt_data <- CreateSeuratObject(counts=wt_filt)

# qc
ad_qc <- subset(ad_data, subset = nCount_RNA<600000)
wt_qc <- subset(wt_data, subset = nCount_RNA<600000)

ad_norm <- NormalizeData(ad_qc)
wt_norm <- NormalizeData(wt_qc)


# scaling
ad_var <- FindVariableFeatures(ad_norm, selection.method="vst", nfeatures=100)
wt_var <- FindVariableFeatures(wt_norm, selection.method="vst", nfeatures=100)

proteins_filt <- rownames(ad_var)

ad_scale <- ScaleData(ad_var, features=proteins_filt)
wt_scale <- ScaleData(wt_var, features=proteins_filt)

ad_pca <- RunPCA(ad_scale, features=VariableFeatures(ad_var))
wt_pca <- RunPCA(wt_scale, features=VariableFeatures(wt_var))

ad_nn <- FindNeighbors(ad_pca, dims=1:10)
ad_clus <- FindClusters(ad_nn, resolution=.35)
wt_nn <- FindNeighbors(wt_pca, dims=1:10)
wt_clus <- FindClusters(wt_nn, resolution=.4)

# umap
ad_umap <- RunUMAP(ad_clus, dims=1:10)
wt_umap <- RunUMAP(wt_clus, dims=1:10)
DimPlot(ad_umap, reduction='umap', cols=c(rgb(1,0,0), rgb(0,0,1), rgb(0,.4,0), rgb(.5,.5,0)),
        label=TRUE) + xlim(-6.5, 6.5) + ylim(-4.5, 4.5)
DimPlot(wt_umap, reduction='umap', cols=c(rgb(1,0,.4), rgb(.4,0,1), rgb(.1,.5,.4), rgb(1,.5,0)),
        label=TRUE) + xlim(-6.5, 6.5) + ylim(4.5, -4.5)

# heatmap
ad_clusters <- FindAllMarkers(ad_umap, only.pos=TRUE, logfc.threshold=.1, min.pct=.1)
top_ad <- ad_clusters %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
DoHeatmap(ad_umap, features=top_ad$gene) + NoLegend()

wt_clusters <- FindAllMarkers(wt_umap, only.pos=TRUE, logfc.threshold=.1, min.pct=.1)
top_wt <- wt_clusters %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
DoHeatmap(wt_umap, features=top_wt$gene) + NoLegend()
