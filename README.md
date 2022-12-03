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

I have ran the analysis and have posted a few plots below.

![MA plot](https://user-images.githubusercontent.com/112113115/204729536-f33d1789-70aa-4197-84f2-efd618330fae.png)


![MA plot 2](https://user-images.githubusercontent.com/112113115/204729700-6a4c301e-1f8c-4246-b5ed-a1c512b6b9b7.png)



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
 
