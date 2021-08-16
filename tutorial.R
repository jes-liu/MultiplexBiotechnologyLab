rm(list=ls())

# downloads ----
# install r - https://cran.r-project.org/bin/windows/base/
# install rstudio - https://www.rstudio.com/products/rstudio/
# cocalc.com

# what is r? ----
# r is a programming language mainly used for statistical analysis and graphing

# what is r-studio? ----
# r-studio, however, is an IDE (integrated development environment) to be able
# to create scripts and longer functions for r
# when you download r, it only downloads the console where all the function
# are executed. console on bottom


# basic controls ----
# this is a comment
# ctrl+enter, run, source


# data classes ----
84.9  # numeric
"char"  # character
T/F/TRUE/FALSE  #boolean


# naming objects ----
# assignment convention & executing that line
a <- 2
a
a <- 1 + 1
b <- 5*a
b <- 5*"a"  # cannot use mathematical operations for a numeric and a non-numeric


# case-sensitive ----
A <- 3  # (a =/= A)


# variable names ----
12obj <- 5  # cannot have a number in the front
_obj <- 6  # cannot have special characters in the front
obj) <- 9  # cannot have some special characters at the end
# only two characters
obj.1 <- 1
obj_2 <- 2


# data structures ----
c()  # combines inputted values into vector or list
?c   # opens the function in the help menu
x <- c(1, 2, 3, 5, 10)  # add multiple values to a variable
y <- c(1:10)  # add consecutive values to a variable
z <- c("red", 12, TRUE, a)  # try to add multiple values of different data
                            # classes to a variable
x + 2  # computes mathematical operation to each individual value
y*2
y^y
z + 2  # cannot add to character class


# vector ----
# 1 column/row of data of 1 data type
vector <- c("red", "blue", "green")  # only one data type
vector
vector <- c(vector, "orange")  # add another value to an existing vector


# list ----
# 1 column/row of data of >1 data type
l <- list("red", 12, TRUE, a)  # >1 data type
l
num <- list(1, 2, 3, 4)
num*2  # will not perform mathematical operation on lists
list(vector, a)  # add another value to an existing list


# subsetting ----
vector <- c("red", "blue", "green", "yellow", "purple", "grey")
vector[5]  # calls the 5th index of the vector
vector[1:4]  # calls the first 4 indicies of the vector
vector[c(1, 3, 6)]  # calls the 1st, 3rd, and 6th indicies of the vector
vector[-5]  # calls every index except the 5th
vector[-c(2, 4)]  # calls every index except the 2nd and the 4th


# matrix ----
# 2+ columns/rows of data of 1 data type
mat <- matrix(1:9, nrow=3, ncol=3)  # creating a simple matrix
mat[1, 2]  # gets the 1st row, 2nd column value
mat[1,]  # gets every column in the 1st row
mat[,2]  # gets every row in the 2nd column


# dataframe ----
# 2+ columns/rows of data of >1 data type
df <- data.frame("Age" = c(18:29),
                 "Gender" = c("Male", "Female", "Male", "Female", "Male", "Female",
                              "Male", "Female", "Male", "Female", "Male", "Female"))
                # creating a simple data frame with "Age" and "Gender" columns
head(df)  # gets the first x rows of the data frame
tail(df)  # gets the last x rows of the data frame
df$Age  # gets the "Age" column
df$Gender  # gets the "Gender" column
df[,1]  # gets every row of the first column
df[,1:2]  # gets every row of the first and second columns
df[1,]  # gets every column of the first row
df[1:2,]  # gets every column of the first and second rows
df$Eyes <- c("Blue", "Brown", "Green", "Blue", "Brown", "Green",
             "Blue", "Brown", "Green", "Blue", "Brown", "Green")
           # creating a new column with data onto an existing data frame
df$Twenties <- df$Age > 19  # creating a new data frame column from an existing
                            # column with conditions
df$BlueEyes <- df$Eyes == "Blue"  # same as previous
df[df$Age >= 25 | df$Gender == "Male",]  # gets every row where the conditions
                                         # for each column values are met
df[df$Age >= 25 & df$Gender != "Male",]  # same as previous
df


# functions ----
values <- 1:1000
mean(values)  # returns the result of the function
m <- mean(values)  # stores the function result into a variable
# median, max, min, sum, sd, log, sqrt, ?/help()
class(values)  # returns the class of the variable
length(values)  # returns how many values are in the variable
rep("Apples", 100)  # repeats the inputted value an x amount of times
summary(values)  # a summary of basic statistical analyses

# statistical tests
before <- sample(150:200, 50, replace=TRUE)  # chooses 50 values within the 
                                             # range given, with replacement
after <- sample(175:225, 50, replace=TRUE)
t.test(before, after, paired=TRUE)  # basic t-test
cor(before, after)  # correlation test

# plots ----

# x-y plots
# arrow keys on "Plots" tab lets you go back and forth between plots
points <- data.frame("x" = c(1:8),
                     "y" = c(2:9))
plot(points$x, points[,2])  # labels
abline(1, 1)  # creates a line on the previous graph with slope and intercept

# other plots
norm <- rnorm(100000, mean=0)  # generates 100000 random points of a normal
                               # distribution around a mean of 0
hist(norm)  # makes a histogram plot
unif <- runif(100000)  # generates 100000 random points of a uniform distribution
hist(unif)

data1 <- rnorm(500, mean=4)
data2 <- rnorm(500, mean=7)
box <- data.frame("one" = data1,
                "two" = data2)
boxplot(box)  # makes a box plot

# overlap
hist1 <- hist(data1)
hist2 <- hist(data2)
plot(hist1)
plot(hist2, add=TRUE)
# add - if TRUE, then plots the current graph onto the previous graph
# xlim - change the x axis limit from a range c(start, end)
# col - color of the plot. this can be character value (i.e. "red") or a
#       rgb value (i.e. rgb(256, 256, 256))


# side by side
par(mfrow=c(2,2))  # changes the number of possible plots on one graph 
plot(hist1)
plot(hist2)
plot(hist1)
plot(hist2)
par(mfrow=c(1,1))  # changes the number of possible plots back to 1

# if/else statments ----
x <- 2
if (x==1) {  # used to run a set of steps if the conditional statement is met
  # code goes here
} else if (x>1) {
  # code goes here
} else if (x<1) {
  # code goes here
} else {
  # code goes here
}


# loops ----
y <- 1
for (x in 1:10) {  # used to repeat a set of steps an x number of times
  # code goes here
}


# creating functions ----
# function(inputs)
our_mean <- function(numbers) {
  # code goes here
  avg <- sum(numbers)/length(numbers)
}

vector <- c(1:20)
our_mean(vector)


# packages ----
# Bioconductor
# Biocmanager

library(dplyr)  # call the packages that you will be using
library(Seurat)
library(SeuratObject)
library(patchwork)


# reading data in r ----
dir <- "C:/Users/jesse/OneDrive/Documents/R/data/"  # location of where the file
                                                    # is stored
list.files(dir)  # show a list of files inside the directory

pbmc_data <- Read10X(dir)  # read the files within the directory onto R
                           # must contain the following: 
                           # "barcodes.tsv", "genes.tsv", "matrix.mtx"

dir_1 <- "C:/Users/jesse/OneDrive/Documents/R/data/matrix.mtx"
dir_2 <- "C:/Users/jesse/OneDrive/Documents/R/data/barcodes.tsv"
dir_3 <- "C:/Users/jesse/OneDrive/Documents/R/data/genes.tsv"
ReadMtx(dir_1, dir_2, dir_3)  # read the files and show the count data


pbmc <- CreateSeuratObject(pbmc_data)  # create a Seurat Object from the data
                                  # so it can be processed using Seurat

head(pbmc)  # data explanation:
            # orig.ident - cell names
            # nCountRNA - total amount of expression for each cell
            # nFeature_RNA - number of genes that cell expresses


# visualization of data ----
VlnPlot(pbmc, features=c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(pbmc, feature1="nCount_RNA", feature2="nFeature_RNA")


# preprocessing ----
# QC
head(sort(pbmc$nCount_RNA, decreasing=TRUE), 20)  # used the sort function to
    # sort the data by descending RNA counts, then took the top 20 using head
pbmc_qc <- subset(pbmc, subset= nCount_RNA<10000 & nFeature_RNA<2200)
    # subset - a function that filters data by creating conditions to be met
# visualize - replot

# normalization
pbmc_norm <- NormalizeData(pbmc_qc)
# visualize - pbmc_norm[["RNA"]]@data


# finding variable features ----
pbmc_var <- FindVariableFeatures(pbmc_norm, nfeatures = 2000)
    # finds the features that has the most variable expressions (i.e. it is
    # highly expression in some cells but lowly expressed in other cells)
    # nfeatures - number of top variable features to find 
pbmc_var  # now introduces variable features which is what we just found

top <- head(VariableFeatures(pbmc_var), 30)

plot1 <- VariableFeaturePlot(pbmc_var)  # plots both the variable and
                                        # non-variable features (color coded)
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)  # label points
    # repel - nicely label the gene names to avoid cluttering
plot2


# dimensional reduction PCA ----
allGenes <- rownames(pbmc_var)
pbmc_scale <- ScaleData(pbmc_var, features=allGenes)  # ScaleData:
    # changes the dataset to mean 0, std 1
    # features - the genes that we want to scale for
# visualize - pbmc_scale[["RNA"]]@scale.data

pbmc_pca <- RunPCA(pbmc_scale, features=VariableFeatures(object=pbmc_scale))
    # PCA - prinicpal component analysis, a dimensional reduction technique
    # reduce dimensionality while minimizing information loss - helps to
    # make data easier to interpret
    # returns a set of PC dimensions
# visualize - pbmc_pca[["pca"]]


# reduce technical noise ----
pbmc_jack <- JackStraw(pbmc_pca)  # finds the statistical significance (p-value)
                                  # for each gene to a principal component
                                  # returns scores
pbmc_score <- ScoreJackStraw(pbmc_jack, dims = 1:20)
    # computes the significance of each jackstraw score
ElbowPlot(pbmc_score)  # plots each principal component's standard deviation
                       # true signal - shown by a high deviation for the PC


# clustering ----
pbmc_nn <- FindNeighbors(pbmc_score, dims=1:10)  # runs the k-nearest neighbor
    # algorithm to find the closest neighbors for each cell in a graph
pbmc_clus <- FindClusters(pbmc_nn, resolution=0.5)  # identify the clusters
    # within the graph based on the k-nearest neighbors

pbmc_umap <- RunUMAP(pbmc_clus, dim=1:10)  # run a UMAP (Uniform Manifold
    # Approximation and Projection) dimensional reduction technique used to 
    # help visualize the clusters
    # constructs a high dimensional graph from kNN and FindClusters, then
    # optimizes it to a low dimensional graph
DimPlot(pbmc_umap, reduction='umap')  # plots the umap


# finding markers ----
pbmc_markers <- FindAllMarkers(pbmc_umap)  # find differentially expressed
    # genes per group (usually clusters)
top10 <- pbmc_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
    # take the top 10 differentially expressed genes in each cluster
DoHeatmap(pbmc_umap, features=top10$gene) + NoLegend()
    # plots a heatmap with no legend
