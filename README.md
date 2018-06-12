# Endothelial Cell Transcriptome Analysis

Project aimed at understanding age-induced endothelial dysfunction through transcriptional reprogramming using both conventional bulk sequencing as well as single cell sequencing of young and aged endothelial cells. 


## Bulk Sequencing

Endothelial cells from either young and aged heart, brain and kidney were isolated using conventional lab protocols. RNA was extracted from all tissues in replicates of 3 and sequenced on the illumina HiSeq. Age specific differences were then identified by comparing each young vasuclar bed with its corresponding aging vascular bed. Differential analysis was done using DeSeq2. Analysis can be found in the folders as listed below:
### BulkAnalysis/Brain
### BulkAnalysis/Heart
### BulkAnalysis/Kidney


## Single Cell Sequencing

Endothelial cells from either young and aged heart and brain were isolated into single cell suspensions using conventional lab protocols. Single cell cDNA generation and tagging was performed on the 10X genomic platform. cDNA from all single cells were sequenced on the Nextera Sequencer and sequences were demultiplexed and assigned to original cell of origin using CellRanger. Clustering of single cells from each tissue was done using the Seurat Software. Combined clustering and analysis of single cells from both young and aged endothelial cells was also done using Seurat. Analysis can be found in the folders listed below:
### SingleCellAnalysis/Heart
### SingleCellAnalysis/Brain



