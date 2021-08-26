rm(list=ls())

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)
library(ggplot2)
library(CellChat)

print("Reading Data")
dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unified_minusBackG.xlsx"
#dir <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Data_unifiedB.xlsx"


# preprocessing ----
# read the data into a dataframe
raw <- data.frame(read.xlsx(dir))

# get data groups
type <- raw$type
label <- raw$Label

# read only the columns with expression data
data_raw <- raw[,4:length(raw)]

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
data_clus <- FindClusters(data_nn, resolution=.3)

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
AD_markers
WT_markers <- FindMarkers(data_umap, ident.1="WT", only.pos=TRUE,
                          group.by="Condition", logfc.threshold=.1)
topAD <- AD_markers %>% top_n(n=30, wt=avg_log2FC)
topWT <- WT_markers %>% top_n(n=30, wt=avg_log2FC)
total_top <- c(rownames(topAD), rownames(topWT))
DoHeatmap(data_umap, features=total_top, group.by="Condition") + NoLegend()


# dotplot ----
dot_markers <- FindAllMarkers(data_umap, test.use = "wilcox", 
                              logfc.threshold=.1, only.pos=TRUE)
#dot_markers_fil <- subset(dot_markers, subset= pct.2<1)
top_dot <- dot_markers %>% group_by(cluster) %>% top_n(10, wt=pct.2)
DotPlot(data_umap, features=unique(top_dot$gene)) + RotatedAxis() + coord_flip()
print("Plots Created")

FeaturePlot(data_umap, features=unique(head(top_dot$gene, 9)))
VlnPlot(data_umap, features=head(top_dot$gene, 9), sort=T, log=T, pt.size=0)


# CellChat ----

data.input <- GetAssayData(data_umap, assay="RNA", slot="data")
labels <- data_umap$Sample
meta <- data.frame(group=labels, row.names=names(labels))

# creating cellchat object from normalized data
cellchat <- createCellChat(object=data.input, meta=meta, group.by="group")

# running on mouse data
DB <- CellChatDB.mouse
showDatabaseCategory(DB)

DB.use <- DB  # use entire database
#DB.use <- subsetDB(DB, search = "Secreted Signaling")

# set the used database in the object
cellchat@DB <- DB.use

# feature names
proteins <- rownames(cellchat@data)

# subset the expression data of signaling genes for saving computation cost
cellchat_pp <- subsetData(cellchat, features=proteins)

cellchat_expr <- identifyOverExpressedGenes(cellchat_pp)
cellchat_int <- identifyOverExpressedInteractions(cellchat_expr)
cellchat_int@data.signaling[1:5, 1:5]

#cellchat_ppi <- projectData(cellchat_pp, PPI.mouse)  # optional

cellchat_prob <- computeCommunProb(cellchat_int)

compareInteractions(cellchat_expr)



