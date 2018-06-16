# Purpose: Differential gene expression analysis between young cardiac endothelial cells and aged cardiac endothelial cells
## Software used : DeSeq2  https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
##Input: Counts of reads generated from HTSeq python script. This is contained in CombinedHeartCounts.csv

## Normalization/Transformation of data: Raw counts were normalized using the median to ratio method whil;e transformation of normalized values was done using a regularized log transformation. Both of these operations are implemented in the DeSeq2 package. 

## An explanation of normalization methods implemented by DeSeq2 package as well as other packages like EdgeR can be seen here https://github.com/hbc/knowledgebase/wiki/Count-normalization-methods

## In addition, these 2 papers (in addition to the DeSeq2 paper "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 2014"), do a pretty good job in comparing the different analysis pipelines:1) Comparison of normalization anddifferential expression analyses usingRNA-Seq data from 726 individualDrosophila melanogaster 2016 2) A Comparison of Methods: Normalizing High-Throughput RNA Sequencing Data 2015 

## Result of analysis: The script for this analysis is called DeSeqScript.Rmd. The result of this analysis is a csv file called Heart_Counts_p_values.csv.

## Data Dictionary for the csv file

### Gene_id : Ensemble numbers for Genes
### Gene_symbol : Abbreviated Gene names
### YH_0,YH_1,YH_2,OH_0,OH_1,OH_2 : Raw counts of genes for each sample where YH represents endothelial cells from young heart and OH represents endothelial cells from ageing heart.
### YH_0_norm,YH_1_norm,YH_2_norm,OH_0_norm,OH_1_norm,OH_2_norm : Normalized counts
### YH_0_lognorm,YH_1_lognorm,YH_2_lognorm,OH_0_lognorm,OH_1_lognorm,OH_2_lognorm : Regularized Log transformation of normalized counts
### baseMean : Mean expression of gene across all samples
### log2foldchange : Fold change in log2 
### pvalue 
### padj : Multiple test correction of pvalues

### Identification of transcription motifs enriched in genes dysregulated due to aging of cardiac endothelial cells. Used this technique and software from this paper: Identify transcription factor binding motifs enriched on a gene list Aibar et al 2016

