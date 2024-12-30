# Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(airway)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load count data and sample info
counts_data <- read.csv('counts_data.csv', row.names = 1)
sample_data <- read.csv('sample_info.csv', row.names = 1)

# Ensure column names match row names
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = sample_data, 
                              design = ~ dexamethasone)

# Pre-filtering: Remove rows with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Set factor levels for the condition
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# Run DESeq analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05) # Using a significance level of 0.05

# Annotate results with gene names and symbols
anno <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), 
                              columns = c("ENSEMBL", "SYMBOL", "GENENAME"), 
                              keytype = "ENSEMBL")
res <- as.data.frame(res)
res <- cbind(ENSEMBL = rownames(res), res)
res <- merge(res, anno, by = "ENSEMBL", all.x = TRUE)
res

# Save results to a CSV file
write.csv(res, "DESeq2_results_annotated.csv")

# Dispersion plot
png("Dispersion_plot.png", width = 800, height = 600)
plotDispEsts(dds)
dev.off()

# MA plot
png("MA_plot.png", width = 800, height = 600)
plotMA(res)
dev.off()

# PCA Plot
vsd <- vst(dds, blind = FALSE) # Variance stabilizing transformation
pcaData <- plotPCA(vsd, intgroup = "dexamethasone", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("PCA_plot.png", width = 800, height = 600)
ggplot(pcaData, aes(PC1, PC2, color = dexamethasone)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()
dev.off()

# Heatmap of top variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(vsd)),
         main = "Heatmap of Top Variable Genes",
         filename = "Heatmap_top_genes.png")

# Volcano plot
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
res$pvalue[is.na(res$pvalue)] <- 1
res$padj[is.na(res$padj)] <- 1

# Define significance thresholds
threshold_padj <- 0.05
threshold_logFC <- 1


#Volcano plot
png("Volcano_plot.png", width = 800, height = 600)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < threshold_padj & abs(log2FoldChange) > threshold_logFC), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +  # Color significant genes red, others gray
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted p-value", title = "Volcano Plot") +
    theme_minimal()
dev.off()

# Inspect the top results
head(res,10)
write.csv(head(res,10), "DESeq2_results_annotated_top10_genes_final.csv")
