## Name: Meiheng Liang
## Programming Language: R v4.3.1
## Date: Apr 11, 2024
## Description:

This script is designed to search E.coli Promoter region based on the scoring of sequence using PSSM. This 20-mer PSSM was computed from psudo-frequency which derived from real frequency of regions of E.coli. The weighing matrix was then used to compute the scoring of sequence which was extracted from 400bp upstream and 50bp within the promoter region. Completed scores was further ranked and extracted where top 30 genes are ones with high likelihood as being promoter region.

### Required files: 

argR-counts-matrix.txt
E_coli_K12_MG1655_.400_50

### Required packages:

BiocManager
BiocManager::install("Biostrings")
library(Biostrings)
rentrez # fetch sequence from the sequence file
 
### Execution:
following function

[step by step information to run the script]
---
title: "528Assignment4_PSSM"
author: "Meiheng Liang"
date: "2024-04-10"
output: html_document
---

```{r setup, include=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("rentrez", quietly = TRUE))
    install.packages("rentrez")
BiocManager::install("Biostrings")
library(Biostrings)

```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r prep data}
# Reading a comma-separated .txt file
data <- read.csv("C:/path/to/your/file/argR-counts-matrix.txt", header = FALSE, sep = "\t", fill = TRUE, strip.white = TRUE)
rownames(data) <- data[,1]  # Set the first column as row names
rownames(data) <- toupper(rownames(data))
matrix <- data[,-c(1,2)]

```

## Including Plots
```{r frequency matrix }
matrix.freq <- (matrix ) / apply(matrix, 2, sum)
matrix.freq 
```

```{r pseudo frequency matrix with augmentation of 1}
k <- 1 
Pseudo.freq <- (matrix + k) / apply(matrix + k, 2, sum) #adding 1 to each position, and frequency calculation as well
Pseudo.freq
weights <- log(Pseudo.freq / 0.25) #calculate log-odds ratio of a nucleotides presence at given position
weights #printing scoring matrix
```

```{r data prep: top 30 gene}
inputFilePath <- ("C:/path/to/your/file/E_coli_K12_MG1655_.400_50") # input file with 450 bp length sequence
outputFilePath <- "C:/path/to/your/file/E_coli_K12.fasta" #output fasta file
# Read the raw data
rawfile <- readLines(inputFilePath) 
# Open a connection to the output file
outputConnection <- file(outputFilePath, "w")
# Process each line to remove the back slash and replace with">" for later analysis 
for(line in rawfile) {
  # Split the line into ID and sequence parts
  parts <- strsplit(line, " \\\\ ")[[1]]  # Note the double escape for backslash
  
  # Write the formatted entry to the output file
  cat(sprintf(">%s\n%s\n", parts[1], parts[2]), file = outputConnection)
}

# Close the connection
close(outputConnection)

E.coli <- readDNAStringSet ("C:/path/to/your/file/E_coli_K12.fasta")
head(E.coli)
```

```{r screening: top 30 gene}
scoreSequenceWithPSWM <- function(sequence, pswm) {
  sequence <- as.character(sequence)  # Ensure the sequence is in character format
  scores <- numeric(length = max(0, length(sequence) - ncol(pswm) + 1)) # Corrected to prevent negative length
  for (pos in 1:length(scores)) {
    subseq <- substring(sequence, pos, pos + ncol(pswm) - 1) #indexing sequence at given position
    tempScore <- 0  # Initialize a temporary score variable for this position
    for (i in 1:ncol(pswm)) {
      nucleotide <- substr(subseq, i, i) 
      if(nucleotide %in% rownames(pswm)) { # Add the score for this nucleotide at this position
        tempScore <- tempScore + pswm[nucleotide, i]
      } else {
        tempScore <- NA  # Assign NA if nucleotide is not in PSWM
        break  # Exit the loop for this position
      }
    }
    scores[pos] <- tempScore  # Corrected to use tempScore
  }
  return(scores)
}

scoresList <- t(data.frame(lapply(E.coli, scoreSequenceWithPSWM, pswm=weights))) #connverting 
scoresList <- data.frame(score = as.numeric(scoresList), row.names = rownames(scoresList))
GeneRank <- scoresList[order(scoresList$score, decreasing = TRUE), , drop = FALSE]
top30 <- GeneRank[1:30, , drop = FALSE]
colnames(top30)<-'score'


```

### output file: 
E_coli_K12.fasta
528Assignment4.htm