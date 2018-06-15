# Purpose: Differential gene expression analysis between young brain endothelial cells (BEC) and aged brain endothelial cells
## Software used : DeSeq2  https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## Input: Counts of reads generated from HTSeq python script. This is contained in CombinedBrainCounts.csv

## Normalization/Transformation of data: Raw counts were normalized using the median to ratio method whil;e transformation of normalized values was done using a regularized log transformation. Both of these operations are implemented in the DeSeq2 package. 

## An explanation of normalization methods implemented by DeSeq2 package as well as other packages like EdgeR can be seen here https://github.com/hbc/knowledgebase/wiki/Count-normalization-methods

## In addition, these 2 papers (in addition to the DeSeq2 paper "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 2014"), do a pretty good job in comparing the different analysis pipelines: 1) Comparison of normalization anddifferential expression analyses usingRNA-Seq data from 726 individualDrosophila melanogaster 2016 2) A Comparison of Methods: Normalizing High-Throughput RNA Sequencing Data 2015 

## Result of analysis: The script for this analysis is called DeSeqScript.Rmd. The result of this analysis is a csv file called Brain_Counts_p_values.csv.

## Data Dictionary for the csv file
### Gene_id : Ensemble numbers for Genes
### Gene_symbol : Abbreviated Gene names
### YB1,YB2,YB3,OB1,OB2,OB3 : Raw counts of genes for each sample where YB represents endothelial cells from young brain and OB represents endothelial cells from ageing brain.
### YB1_norm,YB2_norm,YB3_norm,OB1_norm,OB2_norm,OB3_norm : Normalized counts
### YB1_lognorm,YB2_lognorm,YB3_lognorm,OB1_lognorm,OB2_lognorm,OB3_lognorm : Regularized Log transformation of normalized counts
### baseMean : Mean expression of gene across all samples
### log2foldchange : Fold chnage in log2 
### pvalue 
### padj : Multiple test correction of pvalues

## Figure summary of young and aged brain endothelial cell analysis. 
