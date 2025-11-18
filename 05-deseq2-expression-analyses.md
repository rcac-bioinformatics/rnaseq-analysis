---
source: Rmd
title: "Gene-level QC and differential expression (DESeq2)"
teaching: 40
exercises: 45
---

:::::::::::::::::::::::::::::::::::::: questions

- How do we import a featureCounts matrix into R.
- How do we attach sample metadata.
- How do we perform basic QC such as library size checks and PCA.
- How do we construct and run a DESeq2 model.
- How do we extract and visualize differential expression results.

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Load a cleaned gene count matrix from featureCounts.
- Add sample information and prepare a DESeq2 dataset.
- Perform QC using variance stabilizing transforms, PCA, and clustering.
- Run a simple two-group differential expression analysis.
- Create standard output plots such as MA, volcano, and heatmap.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

In the previous episode, we mapped reads using STAR and counted reads per gene using featureCounts.  
In this episode, we take the resulting **gene-level count matrix** and perform:

1. Import and metadata setup  
2. Quality control  
3. Differential expression with DESeq2  
4. Visualization and reporting

This workflow is the standard entry point for downstream RNA-seq interpretation.

::::::::::::::::::::::::::::::::::::::: callout

## What you need for this episode

- `gene_counts_clean.txt` generated from featureCounts  
- a `samples.csv` file describing the experimental groups  
- the `r-rnaseq` module on RCAC

:::::::::::::::::::::::::::::::::::::::

## Step 1: Start Open OnDemand R session


We will be using the OOD to start an interactive session on Scholar:
[https://gateway.scholar.rcac.purdue.edu](https://gateway.scholar.rcac.purdue.edu/)

1. Login using your Purdue credentials after clicking the above link
2. Click on "Interactive Apps" in the top menu, and select "RStudio (Bioconductor)"
3. Fill the job submission form as follows:
   - queue: `scholar`
   - Walltime: `4`
   - Number of cores: `4`
4. Click "Launch" and wait for the RStudio session to start
5. Once the session starts, you'll will be able to click on "Connect to RStudio server" which will open the RStudio interface.




<div class="figure" style="text-align: center">
<img src="fig/05_deseq/open-on-demand.png" alt="Open OnDemand interface"  />
<p class="caption">Open OnDemand interface</p>
</div>


In R, load the required packages:

```r
.localPaths()
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(readr)
library(dplyr)
```

## Step 2: Import gene-level counts

Load the cleaned featureCounts matrix:

```r
setwd("/scratch/negishi/aseethar/rnaseq-workshop")
counts_file <- "results/counts/gene_counts_clean.txt"
counts <- read_tsv(counts_file)
```

Check the first few rows:

```r
head(counts)
```

Extract gene IDs and numeric count matrix:

```r
gene_ids <- counts$Geneid
count_matrix <- counts[, -1]
count_matrix <- as.matrix(count_matrix)
rownames(count_matrix) <- gene_ids
```

## Step 3: Load sample metadata

Assume you prepared a file `samples.csv` containing:

```
sample,condition
SRR33253285,Nat10-CKO
SRR33253286,Nat10-CKO
SRR33253287,Nat10-CKO
SRR33253288,Nat10flox/flox
SRR33253289,Nat10flox/flox
SRR33253290,Nat10flox/flox
```

Load it:

```r
samples <- read_csv("samples.csv", show_col_types = FALSE)
samples
```

::::::::::::::::::::::::::::::::::::::: callout

## Sample order must match the columns of the count matrix

If the metadata rows are not in the same order as the columns in the matrix,  
your model will be incorrect. Always verify column names.

:::::::::::::::::::::::::::::::::::::::

Check consistency:

```r
colnames(count_matrix)
samples$sample
```

If mismatch occurs, reorder:

```r
count_matrix <- count_matrix[, samples$sample]
```

Convert condition to a factor and set the reference level:

```r
samples$condition <- factor(samples$condition)
samples$condition <- relevel(samples$condition, ref = "Nat10flox/flox")
```

## Step 4: Create DESeq2 dataset

```r
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = samples,
                              design = ~ condition)
```

Filter low-count genes:

```r
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
```

## Step 5: Variance stabilizing transform and QC

Transform:

```r
vsd <- vst(dds, blind = FALSE)
```

### Library size plot

```r
lib_sizes <- colSums(counts(dds))
barplot(lib_sizes,
        las = 2,
        main = "Library sizes",
        ylab = "Total counts")
```

### PCA

```r
p <- plotPCA(vsd, intgroup = "condition")
p + ggplot2::ggtitle("PCA of variance stabilized data")
```

### Sample distance heatmap

```r
sample_dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dists),
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample distance heatmap")
```

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: examine QC plots

Using the PCA and distance heatmap:

1. Which samples cluster together.
2. Are the two conditions separated.
3. Is any sample unusually distant from its group.

::::::::::::::::::::::::::::::::::: solution

- Nat10-CKO replicates cluster together.  
- Nat10flox/flox replicates cluster together.  
- No sample appears as an outlier.  

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::


## Step 6: Differential expression analysis

Run DESeq2:

```r
dds <- DESeq(dds)
```

Extract results comparing Nat10-CKO vs Nat10flox/flox:

```r
res <- results(dds, contrast = c("condition", "Nat10-CKO", "Nat10flox/flox"))
summary(res)
```

Save results:

```r
res_tbl <- as.data.frame(res)
write_csv(res_tbl, "DESeq2_results_raw.csv")
```

## Step 7: Shrink log fold changes

Shrink LFC values:

```r
res_shrink <- lfcShrink(dds,
                        coef = "condition_Nat10.CKO_vs_Nat10flox.flox",
                        type = "apeglm")
```

Save:

```r
write_csv(as.data.frame(res_shrink), "DESeq2_results_shrunken.csv")
```

## Step 8: Visualization

### MA plot

```r
plotMA(res_shrink, main = "MA plot (shrunken LFC)")
```

### Volcano plot

```r
res_df <- as.data.frame(res_shrink) %>% mutate(gene = rownames(res_shrink))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), color = "red") +
    ggtitle("Volcano plot")
```

### Heatmap of top DEGs

```r
top_genes <- head(order(res_shrink$padj), 40)
mat <- assay(vsd)[top_genes, ]
pheatmap(mat,
         scale = "row",
         main = "Top 40 DE genes")
```

::::::::::::::::::::::::::::::::::::::: callout

## Interpretation reminder

- Use padj < 0.05 and a meaningful LFC threshold.  
- Shrunken LFC values are more stable for visualization.  
- Always verify sample QC before interpreting DE results.

:::::::::::::::::::::::::::::::::::::::

## Summary

::::::::::::::::::::::::::::::::::::: keypoints

- Use featureCounts output as input for DESeq2.  
- Attach metadata before creating the DESeq2 dataset.  
- Perform QC using library sizes, PCA, and sample distances.  
- Run the DESeq2 model and extract contrasts.  
- Apply LFC shrinkage for stable effect size estimates.  
- Visualize results using MA plots, volcano plots, and heatmaps.

::::::::::::::::::::::::::::::::::::::::::::::::

