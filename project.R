library(data.table)
install.packages("WGCNA")
install.packages("gplots")
install.packages("Hmisc")
install.packages("impute")
install.packages("preprocessCore")

# Load the libraries
library(WGCNA)
library(gplots)
library(Hmisc)
library(impute)
library(preprocessCore)
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(gridExtra)

getwd()

setwd("C:/Users/18136/Desktop/spring 2024/High Throughput/class project/GSE173082_RAW")

#allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------

data <- read.delim("C:/Users/18136/Desktop/spring 2024/High Throughput/class project/GSE173078_rnaseq_raw_counts.txt/GSE173078_rnaseq_raw_counts.txt", header = T)

# get metadata
geo_id <- "GSE173078"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,43:46)]

# prepare data
data[1:10,1:10]

library(dplyr)
library(tidyr)

# Make the row names into a column called 'Gene'
data <- tibble::rownames_to_column(data, var = "Gene")

# Check the result to ensure the new column has been created
head(data)


library(tidyverse)


# Reshape the data to long format
long_data <- data %>%
  gather(key = 'SampleID', value = 'counts', -Gene)

# Prepare phenoData by creating a 'samples' column that matches the 'SampleID' format in long_data
phenoData <- phenoData %>%
  mutate(samples = gsub("_.*", "", title)) # Strips everything after the underscore

# Join the long format data with phenoData
joined_data <- long_data %>%
  inner_join(phenoData, by = c('SampleID' = 'samples'))

# Spread the joined data to wide format
wide_data <- joined_data %>%
  select(Gene, geo_accession, counts) %>%
  spread(key = geo_accession, value = counts)

# Convert 'Gene' column back to row names if required
final_data <- column_to_rownames(wide_data, var = "Gene")
data <- final_data
# Check the result
print(head(data))

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead

# exclude outlier samples
samples.to.be.excluded <- c('GSM5259457', 'GSM5259456', 'GSM5259462')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) # 13284 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 12
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# Relate modules to traits --------------------------------------------------
# module trait associations

# create traits file - binarize categorical variables
traits <- colData %>%
  mutate(
    periodontitis_bin = ifelse(phenotype == 'Periodontitis', 1, 0),
    healthy_bin = ifelse(phenotype == 'Healthy', 1, 0),
    gingivitis_bin = ifelse(phenotype == 'Gingivitis', 1, 0)
  )
select(traits, 8:10)



# binarize categorical variables
colData$severity <- factor(colData$phenotype, 
                           levels = c("Healthy", "Gingivitis", "Periodontitis"),
                           labels = c("None", "Moderate", "Severe"))


severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)


traits <- cbind(traits, severity.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Select only the numeric columns for correlation analysis from traits
numeric_traits <- traits %>% select(periodontitis_bin, healthy_bin, gingivitis_bin)

# Computing the correlation between module eigengenes and numeric trait data
module.trait.corr <- cor(module_eigengenes, numeric_traits, use = "p")

module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

library(CorLevelPlot)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[23:27],
             y = names(heatmap.data)[1:15],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

#Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvals[1:10,1:10]

# Calculate the driver gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)