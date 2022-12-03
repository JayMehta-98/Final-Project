# Final-Project Outline

## Title
Differential gene expression in TCGA within pancreatic cancer between drinkers and non drinkers using DESeq2.

## Author
Jay Mehta

## Overview of project
I will identify diffrentially expressed genes for pancreatic cancer between people who drink vs those who don't. This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
For this analysis, I'll use the TCGA cohort and have identified 185 ht-seq counts files for tumors that fit within my cohort with 143 drinkers and 121 non drinkers and 17 not reported. Within the analysis, I will control for gender.

## Data
I will use the data from https://portal.gdc.cancer.gov/repository. Examining clinical data, there are 143 samples out of which 135 are defined as drinkers, 121 are identified as non-drinkers and 17 are not reported. The specific files are available are here https://github.com/JayMehta-98/Final-Project/tree/master/Clinical%20Cases

## Milestone 1
**Due Date**: Thursday November 22th
I have downloaded 134 samples from the repository and I have written a script which will read all the files and create a count matrix containing the gene expression for all the samples along with the gene id. The star count files can be found over here

https://github.com/JayMehta-98/Final-Project/tree/master/Star%20Count%20Genes

Once the files have been downloaded, run the shell script below on the terminal.

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

**Known Issues**: The file name for each sample had to be manually changed to their corresponding case id and I had to fish out two files which were extra since I had a 143 cases but had 145 files.

## Milestone 2
After preparing the dataset, the major issue I faced was getting the columns of count matrix and rows of column data to match. I had to get rid of the index columns in both the datasets for the vignette to work and to load the data in DESeq 2. I had to modify the vignette to load the data and make it compliant to DESeq2.

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
                              design = ~ condition)
                              
dds
```

**Output**
```
class: DESeqDataSet 
dim: 60660 143 
metadata(1): version
assays(1): counts
rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ...
  ENSG00000288674.1 ENSG00000288675.1
rowData names(0):
colnames(143): TCGA-2J-AAB6 TCGA-2J-AAB8 ... TCGA-YY-A8LH TCGA-Z5-AAPL
colData names(2): alcohol_history condition
```

### Pre filtering
This step is done to reduce the memory size of the dds data object, and increase the speed of the transformation and testing functions within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted. Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total. 

```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
**Output**
```
class: DESeqDataSet 
dim: 46086 143 
metadata(1): version
assays(6): counts mu ... replaceCounts replaceCooks
rownames(46086): ENSG00000000003.15 ENSG00000000005.6 ...
  ENSG00000288674.1 ENSG00000288675.1
rowData names(23): baseMean baseVar ... maxCooks replace
colnames(143): TCGA-2J-AAB6 TCGA-2J-AAB8 ... TCGA-YY-A8LH TCGA-Z5-AAPL
colData names(3): alcohol_history sizeFactor replaceable
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
                     baseMean log2FoldChange     lfcSE       stat      pvalue      padj
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric> <numeric>
ENSG00000000003.15  1855.2253     0.07236326 0.0853062  0.8482770 0.396283720  0.737160
ENSG00000000005.6     13.4887    -1.89301102 0.4471450 -4.2335508 0.000023003  0.012106
ENSG00000000419.13  1426.9202     0.00235198 0.0653884  0.0359694 0.971306743  0.992380
ENSG00000000457.14   717.0735    -0.06781630 0.0786889 -0.8618277 0.388782335  0.731776
ENSG00000000460.17   222.9388     0.08262278 0.0868017  0.9518566 0.341169698  0.696882
...                       ...            ...       ...        ...         ...       ...
ENSG00000288667.1    0.198414       0.248522  0.811166   0.306377   0.7593177        NA
ENSG00000288669.1    0.141307      -0.350842  0.852374  -0.411606   0.6806281        NA
ENSG00000288670.1  219.037515       0.154180  0.092819   1.661087   0.0966959  0.432746
ENSG00000288674.1    5.198035       0.326996  0.169878   1.924884   0.0542439  0.353417
ENSG00000288675.1   29.090789       0.283183  0.139729   2.026655   0.0426977  0.324511
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
ENSG00000000003.15  1855.2253    1.48507e-05 0.00144253 0.396283720  0.737160
ENSG00000000005.6     13.4887   -9.50690e-06 0.00144270 0.000023003  0.012106
ENSG00000000419.13  1426.9202    1.52006e-04 0.00144638 0.971306743  0.992380
ENSG00000000457.14   717.0735    2.29001e-05 0.00144254 0.388782335  0.731776
ENSG00000000460.17   222.9388   -2.25623e-05 0.00144259 0.341169698  0.696882
...                       ...            ...        ...         ...       ...
ENSG00000288667.1    0.198414    9.37736e-07 0.00144269   0.7593177        NA
ENSG00000288669.1    0.141307   -1.29785e-06 0.00144269   0.6806281        NA
ENSG00000288670.1  219.037515    1.82277e-05 0.00144258   0.0966959  0.432746
ENSG00000288674.1    5.198035    1.15593e-05 0.00144267   0.0542439  0.353417
ENSG00000288675.1   29.090789    2.03537e-05 0.00144269   0.0426977  0.324511
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
LFC > 0 (up)       : 432, 0.94%
LFC < 0 (down)     : 136, 0.3%
outliers [1]       : 0, 0%
low counts [2]     : 16088, 35%
(mean count < 1)
```
To check how many adjusted p-values were less than 0.1, run
```
sum(res$padj < 0.1, na.rm=TRUE) 
```
**Output**
``` 
# 568
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
![Plot count](https://user-images.githubusercontent.com/112113115/205466740-589b8a13-2b87-4527-b2c0-d4a3cd84953f.png)
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
  ![GGplot cust](https://user-images.githubusercontent.com/112113115/205464809-04ba625f-8707-4c33-bba6-76cffa729858.png)
  
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

