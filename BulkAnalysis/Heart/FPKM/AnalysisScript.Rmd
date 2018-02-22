---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
# http://genomicsclass.github.io/book/pages/intro_to_highthroughput_data.html
# https://www.bioconductor.org/help/workflows/RNAseq123/
# RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR

# Creating the three table model
```{r}
# Read in the sample_sheet
sample_sheet = read.delim('~/Desktop/YoungCECvsOldCEC/FPKM/ExpressionMatrix/samples.table', row.names=1)
# Read in the gene annotations
gene_annotations = read.delim('~/Desktop/YoungCECvsOldCEC/FPKM/ExpressionMatrix/genes.attr_table', row.names=1)
# Read in the gene expression matrix
fpkm_matrix <- read.delim("~/Desktop/YoungCECvsOldCEC/FPKM/ExpressionMatrix/genes.fpkm_table", row.names=1)

# Use the sample_sheet, gene annotations and gene expression matrix to build the three tables for modelling

# I noticed that the order of the gene expression table is out of sync with the sample table.
fpkm_matrix[1:3,1:7]
rownames(sample_sheet)

# fpkm_matrix column order is YH_0,YH_2,YH_3,YH_1,OH_1,OH_0 and OH_2 while the sample_sheet table rownames order is "YH_0" "YH_1" "YH_2" "YH_3" "OH_0" "OH_1" "OH_2". Fix that by shuffling the gene expression matrix.
# Also,YH3 tunred out to be a duplicate, so we drop YH3
# Note do this for kidney but dropping OK1
# Use dplyr library
library(dplyr)
fpkm_matrix1 = fpkm_matrix %>% select(YH_0 ,YH_1 ,YH_2,OH_0 ,OH_1 ,OH_2)
# Sanity check
fpkm_matrix1[1:3,1:6]

# The sample_sheet contains a column called file which is a link to the pwd of the gene expression quantification.
# Lets create a different dataframe called sampleInfo without that file column
sampleInfo = sample_sheet %>% select(-file)
# Sanity check
head(sampleInfo)
# Now add another column called sample_names and this column should allow us to connect the rows of this table to the
# columns of the gene expression matrix. 
sampleInfo$sample_names = rownames(sampleInfo)
# Add another column called the phenotype column that will describe each sample as either young or aged
sampleInfo$phenotype = c('Young','Young','Young','Young','Aged','Aged','Aged')
# Cool now we have a table describing the samples in the gene expression matrix called the sampleInfo
# Note that YH_3 was dropped, so drop the YH_3 in the sampleInfo
sampleInfo_2 = sampleInfo %>% filter(sample_names!='YH_3')
# Add back the rownames
rownames(sampleInfo_2) = sampleInfo_2$sample_names
# Sanity check
sampleInfo_2
# Lets make sure the connecting column matches
match(sampleInfo_2$sample_names,colnames(fpkm_matrix1))
all(sampleInfo_2$sample_names == colnames(fpkm_matrix1))

# Get the table that describes the features of the gene expression matrix
# That is in the table called gene_annotations. The column that connects this table to the row names of the gene expression table is the gene_id. Lets test that
# Not going to use the match function since the vectors involved are very long. So lets use this:
all(rownames(fpkm_matrix1)==gene_annotations$gene_id)
# Everything checks out as this boolean returned True

# However, the gene_id are not informative to use unless we look them up. Switch them to actual gene names which we already have as gene_short_name in the gene_annotation table. Make rownames of gene expression matrix to be gene_short_name
# It turns out that there some duplicates in the gene name suggesting isoforms? 
# That was a surprise as I thought I was analysing gene level expression and not isoform.  Went back to the documentation and genes_fpkm.table seems to be a gene level analysis.
# Will look into this more for clarification (Question for Andrew later)
# rownames(fpkm_matrix1) = gene_annotations$gene_short_name
# So lets leave the table as it is for now. I can import the gene_short names as a new column in the expression_matrix down the line.
fpkm_matrix2 = as.data.frame(fpkm_matrix1)
fpkm_matrix2$gene_names = gene_annotations$gene_short_name
# Filter out genes that are lowly expressed. I define lowly expressed as an gene whose mean expression in the data set is greater than 1
filt_edata = subset(fpkm_matrix2,rowMeans(fpkm_matrix2[,1:6]) > 1)
# This operation downsized the number of genes from ~40 k genes to 13 k genes
# Log2 transform the data
LogData = log2(filt_edata[,1:6]+1)
# Make sure this is in matrix form
LogData2 = as.matrix(LogData)


# Full data set
LogData3 = log2(fpkm_matrix1+1)
LogData4 = as.matrix(LogData3)
```

# Exploratory Data Analysis
```{r}
# Create a directory to save the results of the EDA
dir.create('EDA_results')
# Use the table function to summarize the phenotype
table(sampleInfo_2$phenotype)
summary(LogData2)
# Plot the data using a boxplot
boxplot(LogData2,col=2,range=0)
# Also plot using histogram
hist(LogData2,col=2)
# And density plot
plot(density(LogData2),col=2)
# These plots didnt raise any flags
# Use highly expressed genes to create heatmap (Highly expressed is defined as fpkm of greater than 32 ie log2 to the power of 5)
e_heatmap = LogData2[rowMeans(LogData2) >5,]
dim(e_heatmap)
# This downsampled genes from ~9K to ~3K 
# Now plot the heat map. This clusters the data
heatmap(as.matrix(e_heatmap))
# I am having issues with YH_0, it does cluster away. That is something we need to think about
# The default color sucks. So we can specify our own custom color 
colramp = colorRampPalette(c(3,"white",2))(9)
# Then redraw the heat map but now setting the col to the custom color we made. 
heatmap(as.matrix(e_heatmap),col=colramp)
# Note that once we draw the heat map, it clusters the samples based on similarity
# However, if you are not really interested in the clustering based on similarity and you
# want to maintain the order in which the samples were rendered, you can turn off the clustering 
# To do that, set Rowv equals NA and Colv equals NA
heatmap(as.matrix(e_heatmap),col=colramp,Rowv=NA,Colv=NA)

# Perform hierichical clustering
library(rafalib)
distance_matrix <- dist( t(LogData2) )
hc <- hclust(distance_matrix)
plot(hc)
# This plot seem to indicate a batch effect. 
# Note that we filtered the genes from ~40K to ~9K
# Lets look at the clustering with the entire gene set
distance_matrix1 <- dist( t(LogData4) )
hc1 <- hclust(distance_matrix1)
plot(hc1)

# Distribution of genes across all samples as shown by a density plot
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(LogData2[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:6){lines(density(LogData2[,i]),lwd=3,col=colramp[i])}
# All the samples seems to have a pretty similar distribution of genes

# Create a model matrix for the age variable
mod = model.matrix(~ sampleInfo_2$phenotype)
# Then we fit a linear model 
library(limma)
fit_limma = lmFit(LogData2,mod)
# We moderate those linear models
ebayes_limma = eBayes(fit_limma)
topTable(ebayes_limma)
head(ebayes_limma$t)

# ebayes_limma$coefficients has 2 columns. The first column is reflective of the average
# expression level in samples. The second column is the fold change (Young - Old) in log2 
# Then we have the p-value in ebayes_limma$p.value[,2]
# Save these three as a dataframe
results_limma = cbind(ebayes_limma$coefficients,ebayes_limma$p.value[,2])
# Change column names
colnames(results_limma) = c('Avg_Expr','Fold_Change','p_value')

# Then add this to the Log2Data
Results = as.data.frame(cbind(LogData2,results_limma))

library(qvalue)
qobj <- qvalue(p = Results$p_value)
Results$q_value = qobj$qvalues

# Add the gene names
Results$gene_names = filt_edata$gene_names
# Make the gene_names as the first column
Results2 = Results %>% select(gene_names,YH_0 ,YH_1 ,YH_2,OH_0 ,OH_1 ,OH_2,Avg_Expr,Fold_Change,p_value,q_value)
# Save file
write.csv(Results2,'Heart_Fpkm_p_values.csv')

```





