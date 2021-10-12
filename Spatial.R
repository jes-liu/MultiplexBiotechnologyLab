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

data_raw <- t(data.matrix(data_raw))

data <- CreateSeuratObject(data_raw)

# Creating Clusters ----

data_norm <- data
allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=29)
allProteins <- rownames(allData)
data_scale <- ScaleData(allData, features=allProteins)
data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))

data_nn <- FindNeighbors(data_pca, dims=1:25)
data_clus <- FindClusters(data_nn)
data_umap <- RunUMAP(data_clus, dims=1:25, min.dist=.1, spread=2, n.neighbors=50)

DimPlot(data_umap, label=TRUE)

FeaturePlot(data_umap, features=allProteins)  # 1400x2200


coords <- data.frame("cell" = 1:length(data_umap$orig.ident),
                     "seurat_clusters" = data_umap$seurat_clusters)
umap_dir = "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/UMAP_clusters.xlsx"
write.xlsx(coords, umap_dir, overwrite=TRUE)


# Average Nearest Neighbor Distance ----
location = "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/cell_spatial.xlsx"
df <- data.frame(read.xlsx(location))

category <- sort(unique(df$cluster_num))

cluster1 = integer()
cluster2 = integer()
euclidean_nn = integer()

for (x in 1:(length(category)-1)) { # filtering the data for first cluster
  first = filter(df, cluster_num == category[x])

  for (v in (x+1):length(category)) {  # filtering the data for second cluster
    second = filter(df, cluster_num == category[v])
    print(x)
    print(v)
    avg_pair = integer()
    all_cells = integer()
  
    for (y in 1:nrow(first)) { # getting first point's position and creating empty vector
      x0 = first$x.position[y]
      y0 = first$y.position[y]
      one_cell = integer()
      
      for (z in 1:nrow(second)) { # finding distance from first point to all second points
        x1 = second$x.position[z]
        y1 = second$y.position[z]
        dist = sqrt((x1-x0)^2 + (y1-y0)^2)
        one_cell = append(one_cell, dist)
      }
      min_dist = min(one_cell)  # nearest neighbor for each first's cell
      all_cells = append(all_cells, min_dist)  # add a cell's nn to a compilation
    }
    avg_pair = mean(all_cells)  # average of all nearest neighbors for both clusters
    euclidean_nn = append(euclidean_nn, avg_pair)  # compiling all nn to a vector
    cluster1 = append(cluster1, x)
    cluster2 = append(cluster2, v)
  }
}
cluster1 = cluster1 - 1  # cluster nomenclature starts at 0
cluster2 = cluster2 - 1

spatial_distance = data.frame('first_cluster' = cluster1,
                              'second_cluster' = cluster2,
                              'avg_nn_distance' = euclidean_nn)

spatial_dir = "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/spatial_distance.xlsx"
write.xlsx(spatial_distance, spatial_dir, overwrite=TRUE)
