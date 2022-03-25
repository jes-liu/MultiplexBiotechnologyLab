##########################################
rm(list=ls())

library(limma)
library(Seurat)
library(openxlsx)
library(dplyr)
library(edgeR)
library(ggplot2)

dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4_Final1.xlsx'
#zero_dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4 Zero cells adjusted.xlsx'

col_count = c(50, 50, 50, 50, 50, 50, 21, 50, 50, 44)
membrane_count = c(50, 50, 50, 50, 50, 50, 19)
intra_count = c(2, 50, 50, 44)

data = data.frame(read.xlsx(dir, sheet=1))
#zero = data.frame(read.xlsx(zero_dir))
membrane_final = data[,2:(sum(membrane_count)+1)]
intra_final = data[,(sum(membrane_count)+2):(sum(membrane_count)+sum(intra_count)+1)]
#zero_final = zero[,2:(sum(membrane_count)+1)]  # zero cells for membrane proteins only


batch_num = integer()
for (x in 1:(length(col_count))) {
  batch_num = append(batch_num, rep(x, col_count[x]))
}
#batch_num = append(batch_num, rep(7, 2))

#membrane_log = log2(membrane_final+1)
#intra_log = log2(intra_final+1)

membrane_DGE = DGEList(membrane_final)
#membrane_DGE$samples$lib.size = membrane_DGE$samples$lib.size + 1
membrane_norm = calcNormFactors(membrane_DGE)
membrane_cpm = cpm(membrane_norm, log=TRUE)
membrane_limma = removeBatchEffect(membrane_cpm, batch_num[1:sum(membrane_count)])
membrane_subtract = membrane_limma - rowMeans(membrane_limma)  # normalization step

intra_DGE = DGEList(intra_final)
intra_norm = calcNormFactors(intra_DGE)
intra_cpm = cpm(intra_norm, log=TRUE)
intra_limma = removeBatchEffect(intra_cpm,
                    batch_num[(sum(membrane_count)+1):(sum(membrane_count)+sum(intra_count))])
intra_subtract = intra_limma - rowMeans(intra_limma)  # normalization step
# 
# zero_limma = removeBatchEffect(zero_final, batch_num[1:sum(membrane_count)])


intra_limma_df = data.frame(intra_subtract)
#zero_limma_df = data.frame(zero_limma)

membrane_limma_df = data.frame(membrane_subtract)


membrane_limma_df$Cell = c(1:nrow(data))
intra_limma_df$Cell = c(1:nrow(data))

combined = merge(membrane_limma_df, intra_limma_df, by="Cell")

membrane_limma_df = membrane_limma_df[-length(membrane_limma_df)]
intra_limma_df = intra_limma_df[-length(intra_limma_df)]
combined = combined[-1]

list_df = list("membrane" = membrane_limma_df, "intracellular" = intra_limma_df)

# write.xlsx(list_df,'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4_limma_separated.xlsx',
#             overwrite=TRUE)
# write.xlsx(combined,
#           'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4_limma_combined.xlsx',
#            overwrite=TRUE)



##################### seurat after limma
#data_raw= data.matrix(membrane_limma_df)
data_raw = t(data.matrix(intra_limma_df))

#row.names(data_raw) = row.names(membrane_limma_df)  # for combined

data = CreateSeuratObject(data_raw)
#data$Batch = batch_num[1:sum(membrane_count)]  # for membrane df
#data$Batch = batch_num[(sum(membrane_count)+1):(sum(intra_count)+sum(membrane_count))]  
# for intra df
#data$Batch = batch_num


#VlnPlot(data, features=c("nCount_RNA", "nFeature_RNA"))
#FeatureScatter(data, feature1="nCount_RNA", feature2="nFeature_RNA")
#data_qc <- subset(data, subset= nCount_RNA<27000)


#data_norm = NormalizeData(data)
data_norm = data
allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=100)
allProteins <- rownames(allData)
data_scale <- ScaleData(allData, features=allProteins)
data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))
data_nn <- FindNeighbors(data_pca, dims=1:15)
data_clus <- FindClusters(data_nn, resolution=.5)
data_umap <- RunUMAP(data_clus, dims=1:15, n.neighbors=50, spread=2, min.dist=.001)
DimPlot(data_umap, reduction='umap', label=TRUE, pt.size=1, label.size=6) + FontSize(
  x.text=20, y.text=20, x.title=20, y.title=20) + theme(legend.text=element_text(size=20))
#DimPlot(data_umap, reduction='umap', label=TRUE, group.by="Batch")

# membrane_limma_df$Cell = c(1:nrow(membrane_limma_df))
# membrane_limma_df$seurat_clusters = as.numeric(data_umap$seurat_clusters) + 1
# write.xlsx(membrane_limma_df,
#            'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/membrane_limma.xlsx',
#            overwrite=T)


markers <- FindAllMarkers(data_umap, test.use="wilcox")
sort(markers$avg_log2FC, decreasing=T)
#markers[markers$gene=="CD45RA" | markers$gene=="CD45RO" | markers$gene=="CD25",]
top <- markers %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
#bottom <- markers %>% group_by(cluster) %>% top_n(n=-5, wt=avg_log2FC)
#genes = unique(c(top$gene, bottom$gene))
DoHeatmap(data_umap, features=top$gene) + NoLegend() + theme(
  text = element_text(size = 15, face = "bold"))
 
#boxplot(membrane_limma_df)

