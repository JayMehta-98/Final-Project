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
I have downloaded 134 samples from the repository and I have written a script which will read all the files and create a count matrix containing the gene expression for all the samples along with the gene id.

**Known Issues**: The file name for each sample had to be manually changed to their corresponding case id and I had to fish out two files which were extra since I had a 143 cases but had 145 files.

## Milestone 2
**An initial completion of vignette.** I will complete an entire first draft of analysis analyzed through the vignette.Data loaded into vignette (through htseq), for seeking feedback. Not all sections in the writing will be completed, but will be final project.

## Deliverable
**Due Date**: December 3rd

A complete repository with clear documentation and description of your analysis and results.
