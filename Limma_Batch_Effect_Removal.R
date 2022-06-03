rm(list=ls())

"
This program takes in an excel file of protein expression numerical data only (no strings).
This program is typically used after preprocessing and removing the artifical correlation.
Visualization plots have been added to valdiate that there is indeed batch effect present
in the data.
"

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(openxlsx)
library(limma)
library(edgeR)


## Limma Batch Effect Removal ----


# Read the file and convert to data frame
input_path = "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/NewADvsWT/Data Unified 2SD.xlsx"
data = data.frame(read.xlsx(input_path))


# Visualization of the raw data
new_data = data
new_data[new_data == 0] = 1
new_data = apply(log2(new_data), 2, mean)
plot(new_data)  # check to make sure there is visible batch effect before proceeding

par(mfrow=c(2,1))  # changes the number of plots to 2 rows and 1 column for comparison later on
boxplot(data)  # simple boxplot


# Get the number of proteins (columns) per batch
# This needs to be change depending on how many batches there are and how many proteins are in
# each batch. c(47, 45, 48, 42) means that there are 4 batches with 47, 45, 48, and 42 proteins
batches = c(47, 45, 48, 42)

# Create and assign a batch number for each protein column based on its batch number
batch_num = integer()
for (x in 1:(length(batches))) {
  batch_num = append(batch_num, rep(x, batches[x]))
}


# Go through the limma processing by log normalizing the data then run it through the function
limma = DGEList(data)
limma = calcNormFactors(limma)
limma = cpm(limma, log=TRUE)
limma = removeBatchEffect(limma, batch_num)  # input here must be log-normalized

# Visualize the limma-processed data
boxplot(limma)
par(mfrow=c(1,1))


# Run Seurat visualization of the cluster of the proteins
data_by_proteins = data.matrix(limma)

dataProtein = CreateSeuratObject(data_by_proteins)

dataProtein$Batch = batch_num
dataProtein = NormalizeData(dataProtein, normalization.method = "RC")
dataProtein_allData <- FindVariableFeatures(dataProtein, selection.method="vst", 
                                            nfeatures=nrow(data_by_proteins))
dataProtein_allProteins <- rownames(dataProtein_allData)
dataProtein <- ScaleData(dataProtein_allData, features=dataProtein_allProteins)
dataProtein <- RunPCA(dataProtein, features=VariableFeatures(dataProtein_allData))
dataProtein <- FindNeighbors(dataProtein, dims=1:15)
dataProtein <- FindClusters(dataProtein, resolution=.5)
dataProtein_umap <- RunUMAP(dataProtein, dims=1:15)
DimPlot(dataProtein_umap, group.by='Batch')

# Output UMAP should have all proteins in those batches to be well mixed


# Writing the limma-processed excel file for the next analysis
write.xlsx(limma,'Data_limma_processed.xlsx', overwrite=TRUE)
