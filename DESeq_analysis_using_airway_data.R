# script to perform differential gene expression analysis using DESeq2 package

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# Step 1: preparing count data 

# counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)


# sample info
sample_data <- read.csv('sample_info.csv')
head(sample_data)


# making sure the row names in sample_data matches to column names in counts_data
all(colnames(counts_data) %in% rownames(sample_data))

# same order
all(colnames(counts_data) == rownames(sample_data))


# Step 2: constructing a DESeqDataSet object 

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = sample_data,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts

# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")


# Step 3: Run DESeq 
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results 

summary(res)

res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

# contrasts
resultsNames(dds)

#Dispersion Estimates 
png("Dispersion plot.png", width = 800, height = 600)
plotDispEsts(dds)
dev.off()

# MA plot
png("MA plot.png", width = 800, height = 600)
plotMA(res)
dev.off()

