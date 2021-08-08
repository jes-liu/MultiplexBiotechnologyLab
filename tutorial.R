rm(list=ls())
# what is r? ----
# r is a programming language mainly used for statistical analysis and graphing

# what is r-studio? ----
# when you download r, it only downloads the console where all the function
# are executed. console on bottom
# r-studio, however, is an IDE (integrated development environment) to be able
# to create scripts and longer functions for r

# naming objects ----
# assignment convention & executing that line
a <- 2
a <- 1 + 1
b <- 5*a

# case-sensitive (a =/= A) ----
A <- 3

# variable names ----
12obj <- 5
_obj <- 6
(obj) <- 9
obj.1 <- 1
obj_2 <- 2


# data classes ----
84.9  # numeric
"char"  # character
T/F/TRUE/FALSE  #boolean


# data structures ----
c()  # combines inputted values into vector or list
?c 
x <- c(1, 2, 3, 5, 10)
y <- c(1:10)
z <- c("red", 12, TRUE, a)
x + 2
y*2
y^y


# vector ----
# 1 column/row of data of 1 data type
vector <- c("red", "blue", "green")
vector


# list ----
# 1 column/row of data of >1 data type
l <- list("red", 12, TRUE, a)
l
num <- list(1, 2, 3, 4)
num*2
is.vector()
is.list()
list(vector, a)
c(list, 36)
c(vector,36)


# subsetting ----
vector <- c("red", "blue", "green", "yellow", "purple", "grey")
vector[1]
vector[1:4]
vector[c(1, 3, 6)]
vector[-5]


# matrix ----
# 2+ columns/rows of data of 1 data type
mat <- matrix(1:9, nrow=3, ncol=3)
mat[1, 2]

# dataframe ----
# 2+ columns/rows of data of >1 data type
df <- data.frame("Age" = c(18:29),
                 "Gender" = c("Male", "Female", "Male", "Female", "Male", "Female",
                              "Male", "Female", "Male", "Female", "Male", "Female"))
head(df)
tail(df)
df$Age
df$Gender
df[,1]
df[,1:2]
df[1,]
df[1:2,]
df$Eyes <- c("Blue", "Brown", "Green", "Blue", "Brown", "Green",
             "Blue", "Brown", "Green", "Blue", "Brown", "Green")
df$Twenties <- df$Age > 19
df$BlueEyes <- df$Eyes == "Blue"
df[df$Age >= 25 | df$Gender == "Male",]
df[df$Age >= 25 & df$Gender != "Male",]
df


# functions ----
values <- 1:1000
mean(values)
m <- mean(values)
# median, max, min, sum, sd, log, sqrt, ?/help()
class(values)
length(values)
rep("Apples", 100)
summary(values)

# statistical tests
before <- sample(150:200, 50, replace=TRUE)
after <- sample(175:225, 50, replace=TRUE)
t.test(before, after)  # var, paired
cor(before, after)

# plots ----

# x-y plots
# arrow keys
points <- data.frame("x" = c(1:8),
                     "y" = c(2:9))
plot(points$x, points[,2])  #type


# rand hist plots
norm <- rnorm(100000, mean=0)
hist(norm)
unif <- runif(100000)
hist(unif)

# overlap
data1 <- hist(rnorm(500, mean=4))
data2 <- hist(rnorm(500, mean=7))
plot(data1)
plot(data2, add=TRUE)
#xlim=c(), col="red",rgb()
# delete later:
#plot(data1, xlim=c(0,12), col=rgb(1,0,0,.25))
#plot(data2, xlim=c(0,12), col=rgb(0,0,1,.25), add=TRUE)

# side by side
par(mfrow=c(2,2))
plot(data1)
plot(data2)
plot(data1)
plot(data2)


# loops ----
x <- 1
if (x==1) {
  # code goes here
} else if (x>1) {
  # code goes here
} else if (x<1) {
  # code goes here
} else {
  # code goes here
}

for (x in 1:10) {
  # code goes here
}


# creating functions ----
# function(inputs)
our_mean <- function(numbers) {
  # code goes here
  avg <- sum(numbers)/length(numbers)
  return(avg)
}

vector <- c(1:20)
our_mean(vector)


# packages ----
# Bioconductor
# Biocmanager
# reticulate

library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)


# reading data in r ----
dir <- "C:/Users/jesse/OneDrive/Documents/R/data/"
list.files(dir)

pbmc_data <- Read10X(dir)

pbmc <- CreateSeuratObject(pbmc_data)

head(pbmc)  # explain data

# visualization of data ----
VlnPlot(pbmc, features=c("nCount_RNA", "nFeature_RNA"), ncol=2)
FeatureScatter(pbmc, feature1="nCount_RNA", feature2="nFeature_RNA")


# preprocessing ----
# QC
head(sort(pbmc$nCount_RNA, decreasing=TRUE), 20)
pbmc_qc <- subset(pbmc, subset= nCount_RNA<10000 & nFeature_RNA<2100)
# visualize - replot

# normalization
pbmc_norm <- NormalizeData(pbmc_qc)
# visualize - pbmc_norm[["RNA"]]@data


# finding variable features ----
pbmc_var <- FindVariableFeatures(pbmc_norm, selection.method="vst")
pbmc_var

top <- head(VariableFeatures(pbmc_var), 30)

plot1 <- VariableFeaturePlot(pbmc_var)
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
plot2


# dimensional reduction PCA ----
allGenes <- rownames(pbmc_var)
pbmc_scale <- ScaleData(pbmc_var, features=allGenes)  # explain features
# mean 0, std 1
# visualize - pbmc_scale[["RNA"]]@scale.data

pbmc_pca <- RunPCA(pbmc_scale, features=VariableFeatures(object=pbmc_scale))
# visualize - pbmc_pca[["pca"]]


# reduce technical noise ----
pbmc_jack <- JackStraw(pbmc_pca)
pbmc_score <- ScoreJackStraw(pbmc_jack, dims = 1:20)
ElbowPlot(pbmc_score)  # true signal


# clustering ----
pbmc_nn <- FindNeighbors(pbmc_score, dims=1:10)
pbmc_clus <- FindClusters(pbmc_nn, resolution=0.5)

pbmc_umap <- RunUMAP(pbmc_clus, dim=1:10)
DimPlot(pbmc_umap, reduction='umap')  # label


# finding markers ----
pbmc_markers <- FindAllMarkers(pbmc_umap)
FeaturePlot(pbmc_umap, 
            features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                         "FCGR3A", "LYZ", "PPBP","CD8A"))
top10 <- pbmc_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(pbmc_umap, features=top10$gene) + NoLegend()
