rm(list=ls())

library(limma)
library(Seurat)
library(openxlsx)
library(dplyr)
library(edgeR)
library(ggplot2)

dir = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/CD4_limma_combined.xlsx'
dir2 = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/Batch_Final_separate.xlsx'


combined_df = data.frame(read.xlsx(dir))
membrane_limma_df = data.frame(read.xlsx(dir2, sheet=1))
intra_limma_df = data.frame(read.xlsx(dir2, sheet=2))
intra_limma_df = intra_limma_df[-1]  # fixes max(nu, nv) error




membrane_raw = data.matrix(membrane_limma_df)
intra_raw = data.matrix(intra_limma_df)


#membrane_raw = t(data.matrix(membrane_limma_df))
#intra_raw = t(data.matrix(intra_limma_df))
#row.names(data_raw) = row.names(membrane_limma_df)  # for combined



#### first two proteins are too similar and conflict for intra
#intra_raw = intra_raw[,-1]



seurat_membrane = CreateSeuratObject(membrane_raw)
seurat_membrane$type = "membrane"
#membrane_sub = subset(seurat_membrane, subset= nCount_RNA>-4000)  # protein
membraneVF <- FindVariableFeatures(seurat_membrane, selection.method="vst", nfeatures=300)

seurat_intra = CreateSeuratObject(intra_raw)
seurat_intra$type = "intracellular"
intraVF = FindVariableFeatures(seurat_intra, selection.method="vst", nfeatures=100)

# combined = merge(seurat_membrane, y=seurat_intra)
# allData = FindVariableFeatures(combined, nfeatures=300)


anchors <- FindIntegrationAnchors(object.list = list(membraneVF, intraVF), dims = 1:20)
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

integrate_scale = ScaleData(integrated)
integrate_pca = RunPCA(integrate_scale)
integrate_nn <- FindNeighbors(integrate_pca, dims=1:15)
integrate_clus <- FindClusters(integrate_nn, resolution=.5)
integrate_umap <- RunUMAP(integrate_clus, dims=1:15, n.neighbors=50, spread=1, min.dist=.001)
DimPlot(integrate_umap, reduction='umap', label=TRUE, group.by="type")

#### making cluster for DEGs using integration
membrane_raw_t = t(data.matrix(membrane_limma_df))
intra_raw_t = t(data.matrix(intra_limma_df))
seurat_membrane_t = CreateSeuratObject(membrane_raw_t)
membraneVF <- FindVariableFeatures(seurat_membrane_t, selection.method="vst", nfeatures=400)
seurat_intra_t = CreateSeuratObject(intra_raw_t)
intraVF = FindVariableFeatures(seurat_intra_t, selection.method="vst", nfeatures=400)
anchors <- FindIntegrationAnchors(object.list = list(membraneVF, intraVF), dims = 1:20)
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
integrate_scale = ScaleData(integrated)
integrate_pca = RunPCA(integrate_scale)
integrate_nn <- FindNeighbors(integrate_pca, dims=1:15)
integrate_clus <- FindClusters(integrate_nn, resolution=.5)
integrate_umap <- RunUMAP(integrate_clus, dims=1:15, n.neighbors=50, spread=1, min.dist=.001)
DimPlot(integrate_umap, reduction='umap', label=TRUE)

markers <- FindAllMarkers(integrate_umap, test.use="wilcox")
sort(markers$avg_log2FC, decreasing=T)
top <- markers %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
DoHeatmap(integrate_umap, features=top$gene) + NoLegend() + theme(
  text = element_text(size = 15, face = "bold"))





############## normal analysis
combined_raw = t(data.matrix(combined_df))
data = CreateSeuratObject(combined_raw)
data_norm = data
allData <- FindVariableFeatures(data_norm, selection.method="vst", nfeatures=300)
allProteins <- rownames(allData)
data_scale <- ScaleData(allData, features=allProteins)
data_pca <- RunPCA(data_scale, features=VariableFeatures(allData))
data_nn <- FindNeighbors(data_pca, dims=1:15)
data_clus <- FindClusters(data_nn, resolution=.7)
data_umap <- RunUMAP(data_clus, dims=1:15, n.neighbors=50, spread=1, min.dist=.001)
DimPlot(data_umap, reduction='umap', label=TRUE, pt.size=1, label.size=6) + FontSize(
  x.text=20, y.text=20, x.title=20, y.title=20) + theme(legend.text=element_text(size=20))
#DimPlot(data_umap, reduction='umap', group.by="type", label=TRUE)

markers <- FindAllMarkers(data_umap, test.use="wilcox")
sort(markers$avg_log2FC, decreasing=T)
top <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top$gene
DoHeatmap(data_umap, features=top$gene) + NoLegend() + theme(
  text = element_text(size = 15, face = "bold"))



######## correlations and covariance
deg = unique(top$gene)
correlation = cor(combined_df[, deg])
length(abs(correlation[correlation > .4]))
covariance = cov(combined_df[, deg])

write.xlsx(correlation,
           'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/correlation.xlsx',
           overwrite=T)
write.xlsx(correlation,
           'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/Batch/correlation.xlsx',
           overwrite=T)



######################
# CD95
# CD147
# CD8a
# CD97
# CD212
# CD253
# CD44
# CD57
# CD196
# CD23
# CD161
# NKB1
# JNK
# CD33
# CD231
# Arp2
# CD45
# CD146
# MEK1
# CD4
# HLA-DR-DP-DQ
# CD90
# CD142
# CD94
# CD150
# CD156c
# Cyclin-D1-D2
# Prkca
# Ube1
# CD226
# CD28
# CD223
# CD141
# CD155
# CD89
# HLA-DR
# EIF2S1
# CD83
# CD148
# CD84
# CD160
# PTEN
# TSC2
# CLTC
# ITPR1
# CD77
###############################

library(ggplot2)
library(dplyr)
library(gridExtra)

path = 'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/ToppGeneMatrix.csv'

pathways = data.frame(read.csv(path))  # 6-51 columns

category = unique(pathways$Category)

# created a sum column with total counts for that specific pathway/function
for (x in 1:nrow(pathways)) {
  pathways$sum[x] = sum(pathways[x,6:(ncol(pathways)-1)])
}

# divide the sum by the number of proteins (42) and get percent of presence
pathways$percent = pathways$sum/42 * 100

pathways = pathways[,c('Category', 'Name', 'sum', 'percent')]

# filter for top 15 of proteins present in each category
pathways = pathways %>% group_by(Category) %>% top_n(n=12, wt=percent)

# control the num of characters shown on x labels
for (x in 1:nrow(pathways)) {
  pathways$Name[x] = substr(pathways$Name[x], 1, 40)
}


plots = vector("list", 6)
for (x in 1:6) {
  print(x)
  category1 = pathways[pathways$Category == category[x],]
  category1 = category1 %>%
              arrange(percent) %>%
              dplyr::mutate(func=factor(Name, levels=unique(Name)))
  
  plots[[x]] = ggplot(category1, 
             aes(x = func, y = percent)) + 
             geom_bar(stat = "identity") +
             coord_flip() + scale_y_continuous(name="Percent") +
             scale_x_discrete(name="Function") +
             theme(axis.text.x = element_text(face="bold"), 
                   axis.text.y = element_text(face="bold")) +
             ggtitle(category[x])

}
marrangeGrob(grobs=plots[1:6], nrow=3, ncol=2)

