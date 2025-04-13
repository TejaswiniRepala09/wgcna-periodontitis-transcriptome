# HP_PROJECT
### Title: Analysis of Gene Expression Patterns in Periodontitis Using Weighted Gene Co-expression Network Analysis (WGCNA)

### Author: Tejaswini Repala

### Date: May 02, 2024

### Programming Language & Version: R 4.1.2

### Description:

This project employs RNA sequencing and Weighted Gene Co-expression Network Analysis (WGCNA) to study periodontitis, aiming to identify gene expression patterns and modules associated with the disease. The analysis includes quality control, normalization, network construction, module identification, and correlation with clinical traits to discover potential biomarkers and therapeutic targets.

### Dependencies:

WGCNA: For network analysis and module detection.

DESeq2: For RNA-seq data normalization and analysis.

tidyverse: For data manipulation and visualization.

GEOquery: To fetch datasets from the GEO database.

gplots, Hmisc: Additional visualization tools.

Installation of Dependencies:
These packages can be installed from CRAN and Bioconductor using:


                  install.packages(c("WGCNA", "gplots", "Hmisc", "impute", "preprocessCore", "tidyverse"))
                      BiocManager::install(c("GEOquery", "DESeq2"))
                      
### Input:

GSE173078_rnaseq_raw_counts.txt: Raw counts of RNA-seq data.

### Output:

Graphical outputs including PCA plots, dendrograms, and heatmaps.
Statistical summaries of gene modules correlated with periodontitis traits.
List of key driver genes potentially relevant for further biomedical research.
