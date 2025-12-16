---
title: "Differential expression using DESeq2"
teaching: 40
exercises: 45
author:
  - Arun Seetharam
  - Michael Gribskov (contributed material)

---

:::::::::::::::::::::::::::::::::::::: questions

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

This script should work for many simple designs ranging from two samples (i.e., control/treatment)
to designs with large numbers of samples, for instance a time course. The first sample group is 
used as the denominator in the fold-change calculations. This version is set up for the simple two
sample case, but you can switch to the code in comments for more complicated setups.

Required packages:
DESeq2 (bioconductor)
ggplot2
reshape2
ashr or apelgm (for shrunken LFC calculation)
scales (volcano plot)

For bioconductor packages such as DESeq2, you need to load using BiocManager.

## Attribution

This section is adapted from materials developed by **Michael Gribskov**, Professor of Computational Genomics & Systems Biology at Purdue University.

Original and related materials are available via the **CGSB Wiki**:
<https://cgsb.miraheze.org/wiki/Main_Page>

Michael Gribskov  v 2.4 updated: 2025 December 12

::::::::::::::::::::::::::::::::::::::: prereq

## What you need for this episode

:::::::::::::::::::::::::::::::::::::::

## Working directory and setup

It's easiest to work in the directory where your data is.

```r
getwd()
setwd("A:/mrg/Dropbox/projects/25plantago")
dir()
````



## Package installation and loading

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
```



## Reading raw read count data

Read in raw read count data from merged salmon output (NumReads column). You may want to adjust
row and column labels to be shorter or more informative. It will be convenient if you make the
first sample group the one that will be the denominator in the log fold change calculation.

Normally this will correctly read the rownames and column names from the counts file.

```r
salmon <- read.delim("JM13_18_2.numreads.out", row.names=1, header=T)

row.names(salmon)
colnames(salmon)
summary(salmon)
```

If the rownames and colnames are OK, skip ahead to raw counts.

Possible adjustments to rows and columns include reordering columns, selecting subsets, manually
moving row or column names, and editing names to remove redundant strings.

Remove suffix `.salmon2` from column names.

```r
names  <- colnames(salmon)
names <- sub(".salmon2","",names)
colnames(salmon) <- names

head(salmon)
```



## Raw counts quality control

The count matrix from the previous section is named `salmon`.

Total counts per sample can be checked to verify agreement with read counts and expected dimensions.

```r
head(salmon)
colSums(salmon)
dim(salmon)
```

Boxplot and density plots of raw counts provide a visual summary of count distributions.

```r
melted <- melt(round(salmon), value.name="count", variable="sample")
ggplot(melted,aes(x=sample, y=count, fill=sample)) + 
    geom_boxplot(show.legend=F ) + 
    scale_y_continuous(trans='log10') +
    ggtitle("Un-normalized Counts") + 
    ylab("log( Count )") + xlab("Sample") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
```

```r
ggplot(melted, aes(x=count,fill=sample)) +
  geom_density(alpha=0.3) +
  scale_x_continuous(trans='log10', limits=c(0.1,100000)) +
  annotation_logticks() +
  ggtitle("Un-normalized Counts")
```

```r
ggplot(melted, aes(x=count,fill=sample)) +
    geom_density(alpha=0.3, show.legend=F) +
    scale_x_continuous(trans='log10', limits=c(0.1,100000)) +
    ggtitle("Un-normalized Counts")
```



## Optional low-count filtering

Optional removal of very low count transcripts. A threshold of approximately 5 reads per sample
is typical, but this is experiment specific. Be careful not to remove genes expressed in only one
group.

```r
cutoff <- 20
table( rowSums(salmon[,]) > cutoff )
cut <- salmon[rowSums(salmon[,]) > cutoff,]
```

```r
melted  <- melt(cut, variable="sample")
title <- sprintf("Totals counts > %d", cutoff)
ggplot(melted, aes(x=value,fill=sample)) +
    geom_density(alpha=0.3) +
    scale_x_continuous(trans='log10', limits=c(0.01,100000000)) +
    annotation_logticks() +
    ggtitle(title)
```

Rename filtered data for downstream analysis.

```r
salmon <- cut
```



## Preliminary normalization with DESeq2

Preliminary normalization uses one pass of DESeq2 to estimate size factors and enable prefiltering.

```r
replicates <- as.factor(c(rep("Sminus", 3), rep("Splus", 3)))

metadata <- data.frame(sampleNames=colnames(salmon), sample=replicates)

count_init <- DESeqDataSetFromMatrix( round(salmon), colData=metadata, design=~sample )
init <- DESeq(count_init)
plotDispEsts(init)
```

```r
summary(counts(init,normalized=T))
```

```r
melted <- melt(round(counts(init,normalized=T)), variable="sample")
ggplot(melted, aes(x=value,fill=Var2)) +
  geom_density(alpha=0.3) +
  scale_x_continuous(trans='log10', limits=c(0.1,100000)) +
  annotation_logticks() +
  ggtitle("Preliminary normalization")
```



## Pre-filtering based on normalized counts

Prefiltering removes poorly measured genes unlikely to show significant differential expression.
Multiple cutoffs are evaluated to guide selection.

```r
trial =c()
svector = c()
cutoff <- c(0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 100)
min_over <- 3
nrep <- 3
nsample <- length(levels(replicates))
for( mincount in cutoff) {
    select <- rowSums(counts(init,normalized=T)[,1:nrep] >= mincount) >= min_over
    for (c in 2:nsample) {
        first <- (c - 1) * nrep + 1
        last <- first + nrep - 1
        select <- select | rowSums(counts(init,normalized=T)[,first:last] >= mincount) >= min_over
    }
    trial[sprintf("%d",mincount)]=list(table(select))
    svector[sprintf("%d",mincount)]=list(select)
}
```

```r
for (cutoff in names(trial)) {
    s <- sprintf("cutoff:%4s     discard:%6d   keep:%6d", cutoff, trial[[cutoff]]["FALSE"], trial[[cutoff]]["TRUE"])
    print(s)
}
```

```r
cutoff_selected <- 20
min_selected = svector[[sprintf("%d",cutoff_selected)]]
table(min_selected)
```

```r
melted <- melt(round(counts(init,normalized=T)[min_selected,]), variable="sample", varnames=c("gene", "sample"))
title <- sprintf("Normalized counts selected for %d samples over %d counts", min_over, cutoff_selected)
ggplot(melted, aes(x=value,fill=sample)) +
  geom_density(alpha=0.3) +
  scale_x_continuous(trans='log10', limits=c(0.1,100000)) +
  annotation_logticks() +
  ggtitle(title)
```



## Final differential expression analysis

Final DESeq2 analysis is performed using the prefiltered counts.

```r
selected <- salmon[min_selected,]

count_filtered <- DESeqDataSetFromMatrix( round(selected), colData=metadata, design=~sample )
filtered <- DESeq(count_filtered)
plotDispEsts(filtered)
summary(counts(filtered,normalized=T))
```



## PCA with variance stabilization

Variance stabilized PCA provides improved separation of biological signal.

```r
se <- SummarizedExperiment(vst(round(counts(filtered, normalized=TRUE))),
                           colData=colData(filtered))

se$sampleNames <- c(rep("C16",12), rep("C31",12))
times <- c(rep("T00",3),rep("T48",3),rep("T06",3),rep("T72",3))
se$time <-  c(rep(times,2))

pcaData <- plotPCA(DESeqTransform(se), intgroup=c("sampleNames", "time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=se$sampleNames, shape=se$time)) +
  ggtitle("Normalized/Filtered Counts - variance stabilizing transform") + 
  scale_shape_manual(values = c(15,16,17,18)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
```



## Differential expression contrasts

```r
nrep <- 3
foreground <- levels(filtered$sample)
```

```r
contrasts <- c()
background = foreground[1]
for (f in foreground[-1]) {
    contrast_name <- paste(f, background, sep="_vs_")
    contrasts <- c(contrasts, contrast_name)
    contrast <- results(filtered, contrast=c("sample", f, background))
    assign(contrast_name,contrast)
    volcano(contrast, 5e-1, log2(2), stattype="padj", symmetric=T )
}
```



## Saving results

```r
output_file = "final.tsv"
write.table(statout, file=output_file, sep='\t',col.names=NA,row.names=T,quote=F)
```



## Summary

::::::::::::::::::::::::::::::::::::: keypoints

::::::::::::::::::::::::::::::::::::::::::::::::

## Appendix: functions

```r
volcano <- function( contrast, p_cut=0.1, lfc_cut=1.0, stattype="svalue",
                     title="Volcano Plot", xspace=5, yspace=1, symmetric=FALSE ) {
    sig <- contrast[[stattype]]<p_cut & !is.na(contrast[[stattype]]<p_cut)
    lfc <- abs(contrast$log2FoldChange)>lfc_cut & !is.na(contrast$log2FoldChange>lfc_cut)
    both <- sig & lfc
    neither <- !sig & !lfc 
    other <- xor(sig,lfc)
    
    if (title == "Volcano Plot") {
        title <- sprintf("Volcano Plot -  %s", contrast_name)
    }
    
    low <- min(contrast$log2FoldChange)
    high <- max(contrast$log2FoldChange)
    if (symmetric) {
        l <- max(abs(high), abs(low))
        limits_x <- c(-l, l)
    } else {
        limits_x <- c(low, high)
    }
    
    p <- ggplot(contrast, aes(x=log2FoldChange, y=log10(.data[[stattype]]))) + 
        ggtitle(title) +
        scale_x_continuous(limits=limits_x) +
        geom_hline(yintercept=log10(p_cut), linewidth=0.25, color="red") +
        geom_vline(xintercept=lfc_cut, linewidth=0.25, color="red") +
        geom_vline(xintercept=-lfc_cut, linewidth=0.25, color="red") +
        geom_point(data=contrast[neither,], color="black", size=1.2) +
        geom_point(data=contrast[other,], shape=21, fill="yellow", size=2, stroke=0.6) +
        geom_point(data=contrast[both,], shape=21, color="black", fill="red", size=3, stroke=0.6)
    
    print(p)
}
```

