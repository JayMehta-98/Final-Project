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

```
cts <- read.csv('ctssort2.csv', check.names=FALSE , row.names=1)
coldata <- read.csv("coldatasort3.csv", row.names=1)
```
I have ran the analysis and have posted a few plots below.

![MA plot](https://user-images.githubusercontent.com/112113115/204729536-f33d1789-70aa-4197-84f2-efd618330fae.png)





## Deliverable
**Due Date**: December 3rd

A complete repository with clear documentation and description of your analysis and results.
