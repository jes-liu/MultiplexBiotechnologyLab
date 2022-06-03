rm(list=ls())

"
This program takes in a preprocessed excel file of single-cell protein expression values.
The format should have cells as the rows and proteins as the columns. The aim of this script 
is to analyze the single-cells to find upregulated proteins. Analysis includes:
Clustering (UMAP)
Heatmap
FeaturePlot
Dotplot

For this specific script, the input file is a single mouse-data type.
This allows for  comparisons betweens proteins only within that single dataset.
"

## Imports and Reading Files ----

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)
library(ggplot2)
library(limma)
library(edgeR)


# Reading the excel file into a data frame
input_path = "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Batch_Final.xlsx"
data = data.frame(read.xlsx(input_path))

# Subset the data to only include protein expression values
data = data[,2:length(data)]
print("Data Read")



## Clustering(UMAP) ----

"
Clustering groups the cells together dependent on the protein expression values of each cell.
This allows us to visualize which cells can be identified by a particular cell type or 
how similar the functionality of the cell types are (data-dependent)
"

data = t(data.matrix(data))

# Initialize Seurat Object
data = CreateSeuratObject(data)

# Visualize the data values
# data@assays$RNA@data  # <- run this code if needed

# Visualize Features
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC for outliers and normalization
#data = subset(data, subset = nCount_RNA<600000)  # need to change nCount_RNA to appropriate value
data = NormalizeData(data, normalization.method = "RC")  # can change method between "RC" and "LogNormalize"

# Find Highly Variable Features:
# Pick out proteins that exhibit high cell-to-cell variation
# i.e. they are highly expressed in some cells, and lowly expressed in others
# nfeatures = the number of proteins to choose
allData = FindVariableFeatures(data, selection.method="vst", nfeatures=ncol(data))


# Scaling the Data and running dimensional reduction (PCA)
allProteins = rownames(allData)
data = ScaleData(allData, features=allProteins)
data = RunPCA(data, features=VariableFeatures(allData))

# Run UMAP clustering
data = FindNeighbors(data, dims=1:15) # dim = number of PCA dimensions to include
data = FindClusters(data, resolution=.5) # resolution = how large the community is

data_umap = RunUMAP(data, dims=1:15)
print("Clustering Completed")
DimPlot(data_umap, reduction='umap', label=TRUE)



## Heatmap ----

"
Heatmaps visualize the top x number of proteins per group (i.e. seurat_clusters, 
conditions, samples, etc.) by representing the protein's logFC (log fold change; 
main weight parameter) value on the heatmap. Parameters can be change to find
downregulated proteins, number of proteins per group, type of group, different
weight parameters)

Colors: yellow = highly expressed; purple = poorly expressed
"

# Find top 10 markers of the seurat_clusters labeled groups using logFC
heat_markers = FindAllMarkers(data_umap, only.pos=TRUE)
top = heat_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(data_umap, features=unique(top$gene)) + NoLegend()




## FeaturePlot ----

"
A plot that shows the location of the expression of a single protein for all proteins.
This plot is useful for determining which protein belong within which cluster and helping 
with distinguishing cell types per cluster. The input data must be subtracted by the
background values to ensure no excessive expression representation in the plots. There must
be 0's within the data to show zero expression, otherwise all plots will be over saturated.
"

# Features= can be changed to a list of specific proteins for convenience
# This function will take a very long time to load if there are >50 proteins chosen
# The resulting image must also be saved with large dimensions (i.e. 3000x6000 for 30 proteins)
FeaturePlot(data_umap, features=allProteins)


## DotPlot ----

"
A plot that shows the percent expression within all clusters.
This plot is useful for determining which protein belong within which cluster and helping 
with distinguishing cell types per cluster. The input data must be subtracted by the
background values to ensure no excessive expression representation in the plots. There must
be 0's within the data to show zero expression, otherwise all plots will be over saturated.
"
  
dot_markers <- FindAllMarkers(data_umap, only.pos=TRUE)

# pct.1 = percent of expression within its own cluster
# pct.2 = percent of expression within global clusters
top_dot <- dot_markers %>% group_by(cluster) %>% top_n(10, wt=pct.2)
DotPlot(data_umap, features=unique(top_dot$gene)) + RotatedAxis() + coord_flip()

