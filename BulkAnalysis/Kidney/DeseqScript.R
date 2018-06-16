---
title: "Differential Gene Expression Analysis: Young KEC vs Aged KEC"
author: "Uchenna Emechebe"
output: html_document
---

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

# REFERENCES
# http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2
# http://www.bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
```{r}
setwd("~/Desktop/Kidney")
# read in count data
CountData = read.csv('CombinedKidneyCounts.csv',header=TRUE, row.names=1)
head(CountData)
library(dplyr)
# Reshuffle the columns 
CountData1 = CountData %>% select(YK1,YK3,YK4,OK1,OK3,OK4)
# Rename the sample names
colnames(CountData1) = c('YK1','YK2','YK3','OK1','OK2','OK3')
# Adding the Metadata. This is the phenotype table describing the samples in the 
# expression data set
row.names = colnames(CountData1)
condition =c('Young','Young','Young','Old','Old','Old')
libType = c('single-end','single-end','single-end','single-end','single-end','single-end')
KidneyMetadata = data.frame(row.names,condition,libType, row.names=1)
# Note the condition variable is supposed to be a factor variable indicating that the experiment has 2 levels: Young and Old. Lets factorize that

# condition = factor(c('Young','Young','Young','Old','Old','Old'))
# Now this has a factor of 2 levels: Old and Young
# Using the HeartCountdata and the condition data we can now
# create the CountDataSet which is the central data structure in Deseq package.

```

# DeSeq Analysis
```{r}
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library( "DESeq2" )

# Build the DeSeq Object
# There are three functions used to create the DeSeq object depending on how the count data was acquired
# 1) DESeqDataSet : Use this when the count data was when you have a SummarizedExperiment object
# 2) DESeqDataSetFromMatrix : Use this when you just have a simple matrix of count values
# 3) DESeqDataSetFromHTSeq : Use this when the output is from the HTSeq python package
# Here I will be using DESeqDataSetFromMatrix for creation of the DeSeq object

# Getting the countdata set
# cds = newCountDataSet(HeartCountData, condition )
# This didnt work because I was using the Deseq tutorial
# Well this is Deseq2...so I got myself a Deseq2 tutorial

# Using Count matrix input for Deseq2 package
head(CountData1)
head(KidneyMetadata)
# Make sure that the row names in the meta data match the col names in
# the count matrix data
all(rownames(KidneyMetadata) %in% colnames(CountData1))
# This should be True
all(rownames(KidneyMetadata)==colnames(CountData1))
# This should be true too
# Now make sure the condition have two levels (Young and Old)
KidneyMetadata$condition
# 2 levels Old and Young
# Now construct a DeSeq2 Data set
dds = DESeqDataSetFromMatrix(countData=CountData1, colData=KidneyMetadata, design = ~ condition)
# Pre-filtering to remove genes with 1 or 0 read
dds = dds[rowSums(counts(dds))>1,]
# I want to explicitly tell Deseq to set Young as the control group
# Thus, set my levels
dds$condition = factor(dds$condition,levels=c('Young','Old'))

```

# Differential analysis
```{r}
# Analyse the DeSeq object with the DESeq function
dds = DESeq(dds)
# Perform independent filtering based on mean of normalized counts.
# The purpose of this is to filter out tests that have little or no chance of
# showing significant evidence. Type 1 error increases with increase in the
# size of data set and so removing tests that are more likely not to be significantly
# different before applying a Type 1 error control will lead to decrease in the size of
# dataset controlled for Type 1 error. Typically, this results in increased detection power at the same experiment-wide type I error.
# Genes with very low counts are not likely to see significant differences typically due to high dispersion.Default is 0.1 alpha. 
res = results(dds, alpha=0.1)
# Reorder res by smallest adjusted p value
resOrdered = res[order(res$padj),]

# How many were less than 0.1 for adjusted p-values
sum(res$padj<0.1,na.rm=TRUE)
# It says 628 genes at an alpha level of 0.1
mcols(res,use.names = TRUE)

# Identify differential expression that is the most statistically significant
# ie one with the smallest p value
idx = which.min(res$pvalue)
counts(dds[idx,])
# That was the unnormalized , raw counts. To get this result but with the normalized count
# set the normalized option to True. 
counts(dds,normalized=TRUE)[idx,]
# res = results(dds,contrast=c("condition","Old","Young"))

# Gene names are ensemble id's. Lets use a custom script from http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2 to add 
# the gene names
# Load the annotation package for Mus Musculus package
library("org.Mm.eg.db")
# View the contents of the annotation database
columns(org.Mm.eg.db)

#ids = list of IDS
#fromKey = key type; toKey = key type we want to convert to
#db = the AnnotationDb object to use.
#ifMultiple = the argument specifies what to do if one source ID maps to several target IDs:
  #should the function return an NA or simply the first of the multiple IDs?
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
   stopifnot( inherits( db, "AnnotationDb" ) )
   ifMultiple <- match.arg( ifMultiple )
   suppressWarnings( selRes <- AnnotationDbi::select( 
      db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
   if( ifMultiple == "putNA" ) {
      duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
      selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
   return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

res$hgnc_symbol <- convertIDs( row.names(res), "ENSEMBL", "SYMBOL", org.Mm.eg.db )
res$entrezid <- convertIDs( row.names(res), "ENSEMBL", "ENTREZID", org.Mm.eg.db )

# DeSeq2 uses the so-called BH adjustment for multiple testing problem. This adjusted pvalue
# will be in padj
sum( res$padj < 0.05, na.rm=TRUE )
# 380 genes were significant after multiple test correction

# Analyze the significant results a little bit more
# We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation:
resSig <- subset(res, res$padj < 0.05 )
# Sort by decreasing fold change to identify genes with the maximum downregulation in the old CEC
# Select just 4 
head( resSig[ order( resSig$log2FoldChange ), ], 4)
# Top 4 genes that are downregulated in aging CEC relative to young CEC include Agt,
# Pcdhb11, Lif
# Sort by increasing fold change to identify genes with the maximum upregulation in the old CEC
# Select just 4
head( resSig[ order( -resSig$log2FoldChange ), ], 4)
 

# Data diagnostics
# MA plot
plotMA( res, ylim = c(-3, 3) )
# Dispersion plot
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
# Histogram of the p value
hist( res$pvalue, breaks=20, col="grey" )
# write.csv( as.data.frame(res), file="results.csv" )
write.csv( as.data.frame(res), file="results1.csv" )
# results.csv only has the results but not the counts for each sample.
# The counts for each sample are in the dataframe called CountData1
# The row names are the common variables between the two tables.
# So I will make the rownames a column in each data set 
# And since Results is a filtered data set and so a subset of CountData1
# I will do an inner join with the Result_data as the left table and the 
# CountData1 as the right table

# Read back in the Results csv data
Result_data = read.csv("results1.csv",row.names = 1)

# For inner join to work, I need a common column. I proposed using Gene id names as the common column
# Create a new column called Gene_id which will be just the row names
Result_data$Gene_id = rownames(Result_data)
# Create the same new column also in the count data
CountData1$Gene_id = rownames(CountData1)
# Now do the inner join to create a new table that has the results as well as the
# count data. Note that the count data is the raw counts
New = inner_join(Result_data,CountData1)


# For a more comprehensive dataset, I wanted to add the normalized count values from DeSeq2
# Get the normalized count data from the DeSeq object
NormalizedCountData <- counts(dds, normalized = TRUE)
# Convert to a dataframe
NormalizedCountData = as.data.frame(NormalizedCountData)
# Add the Gene_id column to NormalizedCountData to create a common column for a join operation
NormalizedCountData$Gene_id = rownames(NormalizedCountData)
# Rename the columns to reflect that these are normalized counts
colnames(NormalizedCountData) = c("YK1_norm","YK2_norm","YK3_norm","OK1_norm","OK2_norm","OK3_norm","Gene_id")
# Use inner join to add the normalized count data to the dataset
New1 = inner_join(New,NormalizedCountData)
# I have the raw counts and the normalized counts in the same data frame. Now time to add the transformed normalized data
# Get the regularized transformed data with the rld function in Deseq
rld = rlog(dds,blind=FALSE)
Log2_transformed = assay(rld)
Log2_transformed = as.data.frame(Log2_transformed)
# Rename the columns to reflect that these are regularized log transformed data
colnames(Log2_transformed) = c("YK1_lognorm","YK2_lognorm","YK3_lognorm","OK1_lognorm","OK2_lognorm","OK3_lognorm")
# I just realized that the data frame that I have now all matched; so I didnt really need to use inner join. A simple
# cbind would suffice
# Confirm that the data I am about to merge together matches 
all(rownames(Log2_transformed) %in% New1$Gene_id)
# This returned True so the two data frames match
# Merge column wise
New2 = cbind(New1,Log2_transformed)
#Select relevant columns for the final result data frame
New3 = New2 %>% select_('Gene_id','hgnc_symbol',"YK1","YK2","YK3","OK1","OK2","OK3","YK1_norm","YK2_norm","YK3_norm","OK1_norm","OK2_norm","OK3_norm","YK1_lognorm","YK2_lognorm","YK3_lognorm","OK1_lognorm","OK2_lognorm","OK3_lognorm", "baseMean","log2FoldChange",'pvalue',"padj")
# Rename the hgnc_symbol as just the Gene_symbol
colnames(New3)[2] = 'Gene_symbol'
# Now save this result
write.csv( New3, file="Kidney_Counts_p_values.csv" )


```





```{r, echo=FALSE}

```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
