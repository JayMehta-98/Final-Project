if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("tidyverse")
library(tidyverse)
install.packages("readr")
library(readr)
BiocManager::install("DESeq2")
library(DESeq2)
cts <- read.csv('ctssort2.csv', check.names=FALSE , row.names=1)
coldata <- read.csv("coldatasort4.csv", row.names=1)
coldata[coldata=="Not"] <- "No"

all(colnames(cts) %in% rownames(coldata))

coldata$alcohol_history <- factor(coldata$alcohol_history)
coldata$gender <- factor(coldata$gender)


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ gender + alcohol_history)
dds

# Pre-Filtering

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# Note on factor levels

dds$alcohol_history <- factor(dds$alcohol_history, levels = c("No","Yes"))
dds$alcohol_history <- droplevels(dds$alcohol_history)

# Differential expression analysis

dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, name="alcohol_history_Yes_vs_No")

# Log fold change shrinkage for visualization and ranking

resultsNames(dds)

BiocManager::install("apeglm")
library("apeglm")
resLFC <- lfcShrink(dds, coef="alcohol_history_Yes_vs_No", type="apeglm")
resLFC

# p-values and adjusted p-values

resOrdered <- res[order(res$pvalue),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE) # How many adjusted p-values were less than 0.1?

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

# Independent hypothesis weighting
BiocManager::install("IHW")
library(IHW)
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
# Exploring and exporting results

# MA-plot

plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-5,5))


# Plot counts

plotCounts(dds, gene=which.min(res$padj), intgroup="alcohol_history")

# For customized plotting

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="alcohol_history", 
                returnData=TRUE)
install.packages("ggplot2")
library("ggplot2")
ggplot(d, aes(x=alcohol_history, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Extracting transformed values

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

# Effects of transformations on the variance

ntd <- normTransform(dds)
library(BiocManager)
BiocManager::install("vsn")
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

# Heat Map
BiocManager::install("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("alcohol_history","gender")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


# sample to sample distance
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$alcohol_history, vsd$gender, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("alcohol_history", "gender"))
