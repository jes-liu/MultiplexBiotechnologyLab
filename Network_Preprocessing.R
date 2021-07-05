

#  Train Data Preprocessing ----
#train_data <- ""
# create 'interaction' column in train_data
#train_df <- data.frame(read_csv(train_data))



# Test Data Preprocessing ----
rm(list=ls())

library(openxlsx)
library(preprocessCore)

# RMA
# step 1: background correction
ad_bgcorrected <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/AD_BGcorrected.xlsx"
wt_bgcorrected <- "C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/WT_BGcorrected.xlsx"

# reading the excel file into a dataframe
ad_bg <- data.frame(read.xlsx(ad_bgcorrected))
wt_bg <- data.frame(read.xlsx(wt_bgcorrected))

# changing the dataframe into a matrix for the functions
ad_matrix <- data.matrix(ad_bg)
wt_matrix <- data.matrix(wt_bg)


# step 2: quantile normalization
ad_norm <- normalize.quantiles(ad_matrix)
wt_norm <- normalize.quantiles(wt_matrix)

# quantile normalization visualization
ad_matrix_log <- log(ad_matrix+1)
boxplot(ad_matrix_log[,1:10], las=2, xlab='Protein', ylab='Intensity (log)',
        main='Without normalization')
boxplot(ad_norm[,1:10], las=2, xlab='Protein', ylab='Intensity')
ad_norm_log <- log(ad_norm+1)
boxplot(ad_norm_log[,1:10], las=2, xlab='Protein', ylab='Intensity (log)',
        main='With normalization')


# step 3: median polish
ad_log2 <- log2(ad_norm)
wt_log2 <- log2(wt_norm)
ad_polish <- rcModelMedianPolish(ad_log2)
wt_polish <- rcModelMedianPolish(wt_log2)

# median polish visualization
boxplot(ad_norm[,1:10], las=2, ylim=c(1000, 10000),
        xlab='Protein', ylab='Intensity', main='No Median Polish')
boxplot(ad_polish$Residuals[,1:10], las=2, 
        xlab='Protein', ylab='Intensity', main='With Median Polish and log2')


# writing df to excel
ad_preprocessed <- data.frame(ad_polish$Residuals)
wt_preprocessed <- data.frame(wt_polish$Residuals)

write.xlsx(ad_preprocessed, 'protein_preprocessed_AD.xlsx',
           col.names=TRUE, row.names=TRUE, append=FALSE)
write.xlsx(wt_preprocessed, 'protein_preprocessed_WT.xlsx',
           col.names=TRUE, row.names=TRUE, append=FALSE)

print('Preprocessing: Done')




# Protein pair analysis ----
rm(list=ls())

# mean, p-value, correlation, distance, std?

# p-value - fisher's exact test: test <- fisher.test(dat); dat in matrix
# FDR <- p.adjust(p, method="BH")


library(openxlsx)


ad_path<-'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/protein_preprocessed_AD.xlsx'
wt_path<-'C:/Users/jesse/OneDrive/Documents/Multiplex Lab/Data/protein_preprocessed_WT.xlsx'

ad_preprocess <- read.xlsx(ad_path)
wt_preprocess <- read.xlsx(wt_path)

ad_matrix <- as.matrix(ad_preprocess)
wt_matrix <- as.matrix(wt_preprocess) # wt pairs have not been made yet


# making a list of protein 1 and protein 2 pairs
ncol <- ncol(ad_matrix)  # 183
protein1 <- integer()
protein2 <- integer()
col.names <- colnames(ad_matrix)
for (x in 1:(ncol-1)) {
  protein1 <- append(protein1, rep(col.names[x], (ncol-x)))
  protein2 <- append(protein2, rep(col.names[(x+1):ncol], 1))
  print(length(protein1))
  print(length(protein2))
}


# microarray analysis
corr_final <- integer()
for (x in 1:(ncol-1)) {
  for (y in (x+1):(ncol)) {
    # pearson correlation
    # 'complete.obs' only uses complete pairs of data in both categories (no NA)
    corr <- cor(ad_matrix[,x], ad_matrix[,y], method='pearson', use='complete.obs')
    corr_final <- append(corr_final, corr)
  }
}


library(limma)
library(edgeR)



ad_pairs <- data.frame(protein1 = protein1,
                       protein2 = protein2,
                       # p.value = pvalue,
                       pearson_corr = corr_final
                       )

write.xlsx(ad_pairs, 'protein_pairs_AD.xlsx', overwrite=TRUE)



