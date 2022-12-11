# Final-Project Outline

## Title
Differential gene expression in TCGA within pancreatic cancer between drinkers and non drinkers with gender as control using DESeq2.

## Author
Jay Mehta

## Overview of project
I will identify diffrentially expressed genes for pancreatic cancer between people who drink vs those who don't. This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
For this analysis, I'll use the TCGA cohort and have identified 143 star counts files for tumors that fit within my cohort with 82 drinkers and 61 non drinkers. Within the analysis, I will control for gender.

## Data
I will use the data from https://portal.gdc.cancer.gov/repository. Examining clinical data, there are 143 samples out of which 82 are defined as drinkers and 61 are identified as non-drinkers. The specific files are available are here https://github.com/JayMehta-98/Final-Project/tree/master/Clinical%20Cases

## Milestone 1
**Due Date**: Thursday November 22th
I have downloaded 143 samples from the repository and I have written a script which will read all the files and create a count matrix containing the gene expression for all the samples along with the gene id. The star count files can be found over here

https://github.com/JayMehta-98/Final-Project/tree/master/Star%20Count%20Genes

The code below will create a txt file with all the gene names in it.
```
cat *any file name*.tsv | awk '{print $1}' | tail -n+7 > test.txt
```


Run the shell script below on the terminal. This script will read all the tsv files you have downloaded, take the unstranded column and create a new txt file with the file name as the header. It will then paste all the columns with the test file created by the code above thus giving a complete matrix.

```
#!/bin/bash

for f in *.tsv; do
cat $f | awk '{print $4}' | tail -n+7 > $f"_1".txt;
done;

for i in *_1.txt; do
echo -e "$i" | cut -d "." -f1 | cat - $i > $i"_2".txt;
done;

paste test.txt *_2.txt > testmatrix.tsv

for j in *_*.txt; do
rm $j;
done;
```

**Known Issues**: The file name for each sample had to be manually changed to their corresponding case id and I had to fish out two files which were extra since I had a 143 cases but had 145 files. Also, I had to add the gender of each sample manually in the dataframe.

## Milestone 2
After preparing the dataset, the major issue I faced was getting the columns of count matrix and rows of column data to match. I had to get rid of the index columns in both the datasets for the vignette to work and to load the data in DESeq 2. I had to modify the vignette to load the data and make it compliant to DESeq2.

You can find the prepared data sets over here - https://github.com/JayMehta-98/Final-Project/tree/master/Matrix%20and%20coldata

### To load the dataset on R Studio

```
cts <- read.csv('ctssort2.csv', check.names=FALSE , row.names=1)
coldata <- read.csv("coldatasort3.csv", row.names=1)
```

To match the columns of cts to the rows of col data, I used python to wrangle the data, Please use the script below on Google Colab or any python based notebook to get the rows and columns in the same order.

### Run this script on Python

```
import pandas as pd 
df1 = pd.read_csv('/path/to/coldata.tsv' , sep='\t')
df2=df1.sort_values(by=['case_submitter_id'] , ascending=True)
df3 = df2.set_index("case_submitter_id")
df3.to_csv("/path/to/coldatasort3.csv")
```

Check to make sure that columns and rows in cts and coldata are identical. If not, use the script above for cts as well. Alternatively, you can run the code given below on R Studio once you have loaded the data to match the columns and rows of the datasets.

```
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
```
To make sure the data you have loaded is formatted according to DESeq2, run the codes given below and make sure both of them give the output ```TRUE```

```
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
```

## Deliverable
**Due Date**: December 3rd

**Link to csv file with HUGO names** - https://drive.google.com/drive/folders/1Sx3c0lasWT9LqI1duLsf06RjpS0RKyDi?usp=sharing

### Loading the Libraries

Before begining to run the vignette on R Studio, we will have to load a few libraries.

```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("tidyverse")
library(tidyverse)
install.packages("readr")
library(readr)
library(DESeq2)
```
For any package not installed for DESeq2, Please use the following syntax to install them.

```
BiocManager::install("Package")
```

### Creating DESeq data object using Count matrix input method

Please use the vignette used in Milestone 2 to load the count matrix and column data. Once the data is loaded, use the code provided below to make the DESeq2 dataset.

```
coldata[coldata=="Not"] <- "No"
coldata$condition <- factor(coldata$alcohol_history)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ gender + alcohol_history)
                              
dds
```

**Output**
```
class: DESeqDataSet 
dim: 60660 143 
metadata(1): version
assays(1): counts
rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1
  ENSG00000288675.1
rowData names(0):
colnames(143): TCGA-2J-AAB6 TCGA-2J-AAB8 ... TCGA-YY-A8LH TCGA-Z5-AAPL
colData names(2): alcohol_history gender
```

### Pre filtering
This step is done to reduce the memory size of the dds data object, and increase the speed of the transformation and testing functions within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted. Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total. 

```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
```
**Output**
```
class: DESeqDataSet 
dim: 46086 143 
metadata(1): version
assays(1): counts
rownames(46086): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1
  ENSG00000288675.1
rowData names(0):
colnames(143): TCGA-2J-AAB6 TCGA-2J-AAB8 ... TCGA-YY-A8LH TCGA-Z5-AAPL
colData names(2): alcohol_history gender
```
As you notice, the rownames reduce from 60660 to 46086.

###  Note on factor levels
 To explicitly set the factors levels, follow the code below rather than letting R choose the factor level alphabetically.
 ```
 dds$alcohol_history <- factor(dds$alcohol_history, levels = c("No","Yes"))
 ```
 
 ### Differential expression analysis
 Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
```
dds <- DESeq(dds)
res <- results(dds)
res
```
**Output**
```
log2 fold change (MLE): alcohol history Yes vs No 
Wald test p-value: alcohol history Yes vs No 
DataFrame with 46086 rows and 6 columns
                     baseMean log2FoldChange     lfcSE       stat      pvalue
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
ENSG00000000003.15  1855.2253     0.07624090 0.0858040  0.8885473 3.74246e-01
ENSG00000000005.6     13.4887    -1.85070266 0.4406193 -4.2002304 2.66644e-05
ENSG00000000419.13  1426.9202     0.00543541 0.0657004  0.0827303 9.34066e-01
ENSG00000000457.14   717.0735    -0.08154333 0.0781599 -1.0432891 2.96814e-01
ENSG00000000460.17   222.9388     0.08214938 0.0873735  0.9402095 3.47110e-01
...                       ...            ...       ...        ...         ...
ENSG00000288667.1    0.198414       0.247456 0.8333935   0.296925   0.7665236
ENSG00000288669.1    0.141307      -0.362868 0.8765822  -0.413957   0.6789053
ENSG00000288670.1  219.037515       0.144392 0.0912022   1.583208   0.1133740
ENSG00000288674.1    5.198035       0.306895 0.1700604   1.804621   0.0711340
ENSG00000288675.1   29.090789       0.299778 0.1401010   2.139731   0.0323765
                        padj
                   <numeric>
ENSG00000000003.15 0.7096810
ENSG00000000005.6  0.0137437
ENSG00000000419.13 0.9795046
ENSG00000000457.14 0.6540050
ENSG00000000460.17 0.6903789
...                      ...
ENSG00000288667.1         NA
ENSG00000288669.1         NA
ENSG00000288670.1   0.452076
ENSG00000288674.1   0.376282
ENSG00000288675.1   0.283922
```
**Note** - We can specify the coefficient or contrast we want to build a results table for, by using the following command.
```
res <- results(dds, name="alcohol_history_Yes_vs_No")
```
## Log fold change shrinkage for visualization and ranking
Shrinkage of effect size is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function lfcShrink using apeglm method.
```
resultsNames(dds)
BiocManager::install("apeglm")
library("apeglm")
resLFC <- lfcShrink(dds, coef="alcohol_history_Yes_vs_No", type="apeglm")
resLFC
```
**Output**
```
log2 fold change (MAP): alcohol history Yes vs No 
Wald test p-value: alcohol history Yes vs No 
DataFrame with 46086 rows and 5 columns
                     baseMean log2FoldChange      lfcSE      pvalue      padj
                    <numeric>      <numeric>  <numeric>   <numeric> <numeric>
ENSG00000000003.15  1855.2253    5.07141e-05 0.00144294 3.74246e-01 0.7096810
ENSG00000000005.6     13.4887   -1.61255e+00 0.47153732 2.66644e-05 0.0137437
ENSG00000000419.13  1426.9202    1.32527e-06 0.00144235 9.34066e-01 0.9795046
ENSG00000000457.14   717.0735   -1.39188e-05 0.00144248 2.96814e-01 0.6540050
ENSG00000000460.17   222.9388    5.76564e-05 0.00144307 3.47110e-01 0.6903789
...                       ...            ...        ...         ...       ...
ENSG00000288667.1    0.198414    9.56677e-07 0.00144269   0.7665236        NA
ENSG00000288669.1    0.141307   -1.40886e-06 0.00144269   0.6789053        NA
ENSG00000288670.1  219.037515    1.80254e-05 0.00144257   0.1133740  0.452076
ENSG00000288674.1    5.198035    1.05446e-05 0.00144266   0.0711340  0.376282
ENSG00000288675.1   29.090789    1.24697e-05 0.00144265   0.0323765  0.283922
```
### p-values and adjusted p-values

```
resOrdered <- res[order(res$pvalue),]
summary(res)
```
**Output**
```
out of 46079 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 406, 0.88%
LFC < 0 (down)     : 127, 0.28%
outliers [1]       : 0, 0%
low counts [2]     : 18768, 41%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```
To check how many adjusted p-values were less than 0.1, run
```
sum(res$padj < 0.1, na.rm=TRUE) 
```
**Output**
``` 
# 533
```
## Results

### MA Plot
```
plotMA(res, ylim=c(-2,2))
```
![MA plot](https://user-images.githubusercontent.com/112113115/204729536-f33d1789-70aa-4197-84f2-efd618330fae.png)

```
# Shrink data
plotMA(resLFC, ylim=c(-5,5))
```
![MA Plot 2](https://user-images.githubusercontent.com/112113115/205462312-0a6e4f0d-0f81-4d28-b5d3-4fc599676b06.png)

### Plot Counts

To examine the counts of reads for a single gene across the groups, we use plot counts.
```
plotCounts(dds, gene=which.min(res$padj), intgroup="alcohol_history")
```
![plotcount gender](https://user-images.githubusercontent.com/112113115/205469874-b3385e74-bbe9-448a-b3c2-8ccc9afefae7.png)
**For customized plotting**- an argument returnData specifies that the function should only return a data.frame for plotting with ggplot.
```
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="alcohol_history", 
                returnData=TRUE)
install.packages("ggplot2")
library("ggplot2")
ggplot(d, aes(x=alcohol_history, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
  ```
![gg plot gender](https://user-images.githubusercontent.com/112113115/205469905-8fc1ca8c-cafd-482f-847d-d8bfdc660da8.png)
  
### Extracting transformed values
```
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
```

### Effects of transformations on the variance
```
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
```
![Rplot](https://user-images.githubusercontent.com/112113115/205465131-cd5718aa-684e-4c95-aa1e-f1ba39f32580.png)
```
meanSdPlot(assay(vsd))
```
![meansd vsd](https://user-images.githubusercontent.com/112113115/205467035-64fd7125-66fa-471c-b4d0-44995d4e6291.png)
```
meanSdPlot(assay(rld))
```
![meansd rld](https://user-images.githubusercontent.com/112113115/205467052-d69b5747-ffec-4df7-89ee-670f372976cd.png)

### Data quality assessment by sample clustering and visualization
```
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("alcohol_history","gender")])
```

### Heatmap of the count matrix
```
BiocManager::install("pheatmap")
library("pheatmap")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![heatmap1](https://user-images.githubusercontent.com/112113115/205470266-c4295782-c94f-4cd8-b3af-d83f7dcfbec1.png)
```
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![heatmap vsd](https://user-images.githubusercontent.com/112113115/205470360-cee8f8ea-2616-492d-93bb-bd9adb052caf.png)
```
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![heatmap rlt](https://user-images.githubusercontent.com/112113115/205479289-4cf6ecb6-a877-4cb7-92ee-3995f4c7bb54.png)

**Heatmap of the sample-to-sample distances**
```
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
```
![sts heatmap](https://user-images.githubusercontent.com/112113115/205470746-06543abb-80ca-46cc-adf2-fc9dfb9fbda5.png)

## Principal Component Analysis Plot
```
plotPCA(vsd, intgroup=c("alcohol_history", "gender"))
```
![PCA](https://user-images.githubusercontent.com/112113115/205470941-b0ac4465-f02b-43d1-a3a4-765bbc3b2315.png)

