---
source: Rmd
title: "B. Differential expression using DESeq2 (Kallisto pathway)"
teaching: 40
exercises: 45
author:
  - Arun Seetharam
  - Michael Gribskov (contributed material)
---

:::::::::::::::::::::::::::::::::::::: questions

- How do we import transcript-level quantification from Kallisto into DESeq2?
- What exploratory analyses should we perform before differential expression testing?
- How do we perform differential expression analysis with DESeq2 using tximport data?
- How do we visualize and interpret DE results from transcript-based quantification?
- What are the key differences between genome-based and transcript-based DE workflows?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Load tximport data from Kallisto quantification into DESeq2.
- Perform quality control visualizations (boxplots, density plots, PCA).
- Apply variance stabilizing transformation for exploratory analysis.
- Run differential expression analysis with appropriate contrasts.
- Create volcano plots to visualize DE results.
- Export annotated DE results for downstream analysis.

::::::::::::::::::::::::::::::::::::::::::::::::


## Attribution

This section is adapted from materials developed by **Michael Gribskov**, Professor of Computational Genomics & Systems Biology at Purdue University.
 
Original and related materials are available via the **CGSB Wiki**:
<https://cgsb.miraheze.org/wiki/Main_Page>

## Introduction

In Episode 04b, we quantified transcript expression using **Kallisto** and summarized the results to gene-level counts using **tximport**. In this episode, we use those gene-level estimates for differential expression analysis with **DESeq2**.

The workflow follows the same general pattern as Episode 05 (genome-based), but with important differences:

1. **Input data**: We use the `txi` object from tximport rather than raw counts from featureCounts.
2. **Count handling**: Kallisto estimates are model-based (not integer counts), and DESeq2 handles this appropriately via `DESeqDataSetFromTximport()`.
3. **Bootstrap support**: Kallisto's uncertainty estimates (from bootstraps) can be used with sleuth for transcript-level DE.

::::::::::::::::::::::::::::::::::::::: callout

## When to use this pathway

Use this transcript-based workflow when:

- You quantified with Salmon, Kallisto, or similar tools.
- You want to leverage transcript-level bias corrections.
- You did not generate BAM files and cannot use featureCounts.

The genome-based workflow (Episode 05) is preferred when you have BAM files and need splice junction information or plan to visualize alignments.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: prereq

## What you need for this episode

- `txi.rds` file generated from tximport in Episode 04b
- A `samples.csv` file describing the experimental groups
- RStudio session via Open OnDemand

If you haven't created the sample metadata file, create `scripts/samples.csv` with:

```
sample,condition
WT_Bcell_mock_rep1,WT_mock
WT_Bcell_mock_rep2,WT_mock
WT_Bcell_mock_rep3,WT_mock
WT_Bcell_mock_rep4,WT_mock
WT_Bcell_IR_rep1,WT_IR
WT_Bcell_IR_rep2,WT_IR
WT_Bcell_IR_rep3,WT_IR
WT_Bcell_IR_rep4,WT_IR
```

Create a results directory:

```bash
mkdir -p results/deseq2_kallisto
```

:::::::::::::::::::::::::::::::::::::::

## Step 1: Load packages and data

Start your RStudio session via Open OnDemand as described in Episode 05, then load the required packages:

```r
library(DESeq2)
library(tximport)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(readr)
library(dplyr)

setwd(paste0(Sys.getenv("SCRATCH"), "/rnaseq-workshop"))
```

Load the tximport object created in Episode 04b:

```r
txi <- readRDS("results/kallisto_quant/txi.rds")
```

Examine the structure of the tximport object:

```r
names(txi)
```

```
[1] "abundance"           "counts"              "length"
[4] "countsFromAbundance"
```

```r
head(txi$counts)
```

```
                      WT_Bcell_mock_rep1 WT_Bcell_mock_rep2 WT_Bcell_mock_rep3 WT_Bcell_mock_rep4
ENSMUSG00000000001.5          473.4426         553.7849         826.9969         676.9101
ENSMUSG00000000003.16           0.0000           0.0000           0.0000           0.0000
ENSMUSG00000000028.16          40.5783          33.9665          55.3055          97.8451
ENSMUSG00000000031.20           0.0000           0.0000           0.0000           0.0000
```

::::::::::::::::::::::::::::::::::::::: callout

## Understanding the tximport object

The `txi` object contains several components:

- **counts**: Gene-level estimated counts (used by DESeq2).
- **abundance**: Gene-level TPM values (within-sample normalized).
- **length**: Average transcript length per gene (used for length bias correction).
- **countsFromAbundance**: Method used to generate counts (default: `"no"`).

When using `DESeqDataSetFromTximport()`, the `length` matrix is automatically used to correct for gene-length bias during normalization. This is why we use the default `countsFromAbundance = "no"` in Episode 04b—DESeq2 handles length correction internally, so pre-scaling counts would apply the correction twice.

:::::::::::::::::::::::::::::::::::::::

Load sample metadata:

```r
coldata <- read.csv(
    "scripts/samples.csv",
    row.names = 1,
    header = TRUE,
    stringsAsFactors = TRUE
)
coldata$condition <- as.factor(coldata$condition)

coldata
```

```
                   condition
WT_Bcell_mock_rep1   WT_mock
WT_Bcell_mock_rep2   WT_mock
WT_Bcell_mock_rep3   WT_mock
WT_Bcell_mock_rep4   WT_mock
WT_Bcell_IR_rep1     WT_IR
WT_Bcell_IR_rep2     WT_IR
WT_Bcell_IR_rep3     WT_IR
WT_Bcell_IR_rep4     WT_IR
```

Verify that sample names match between tximport and metadata:

```r
all(colnames(txi$counts) == rownames(coldata))
```

```
[1] TRUE
```

## Step 2: Create DESeq2 object from tximport

The key difference from the genome-based workflow is using `DESeqDataSetFromTximport()` instead of `DESeqDataSetFromMatrix()`:

```r
dds <- DESeqDataSetFromTximport(
    txi,
    colData = coldata,
    design = ~ condition
)
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use DESeqDataSetFromTximport?

This function:

- Automatically handles non-integer counts from Salmon/Kallisto.
- Preserves transcript length information for accurate normalization.
- Incorporates the `txi$length` matrix to account for gene-length bias.

Using `DESeqDataSetFromMatrix()` with tximport counts would lose this information.

:::::::::::::::::::::::::::::::::::::::

Filter lowly expressed genes using group-aware filtering:

```r
# Group-aware filtering: keep genes with >= 10 counts in at least 4 samples
# (the size of the smallest experimental group)
min_samples <- 4
min_counts <- 10
keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds <- dds[keep, ]
dim(dds)
```

```
[1] 18924     8
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use group-aware filtering?

A simple sum filter (`rowSums(counts(dds)) >= 10`) can be problematic:

- A gene with 10 total counts across 8 samples averages ~1.25 counts/sample—too low to be informative.
- It may remove genes expressed in only one condition (biologically interesting!).

Group-aware filtering (`rowSums(counts >= threshold) >= min_group_size`) ensures:

- Each kept gene has meaningful expression in at least one experimental group.
- Genes with condition-specific expression are retained.
- The threshold is interpretable (e.g., "at least 10 counts in at least 4 samples").

:::::::::::::::::::::::::::::::::::::::

Estimate size factors:

```r
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
```

```
WT_Bcell_mock_rep1 WT_Bcell_mock_rep2 WT_Bcell_mock_rep3 WT_Bcell_mock_rep4
         0.7232591          0.8443123          1.2021389          1.1142267
  WT_Bcell_IR_rep1   WT_Bcell_IR_rep2   WT_Bcell_IR_rep3   WT_Bcell_IR_rep4
         1.0516283          1.1245413          0.9832145          1.0567892
```

## Step 3: Exploratory data analysis

### Raw count distributions

Visualize the distribution of raw counts across samples:

```r
counts_melted <- melt(
    log10(counts(dds) + 1),
    varnames = c("gene", "sample"),
    value.name = "log10_count"
)

ggplot(counts_melted, aes(x = sample, y = log10_count, fill = sample)) +
    geom_boxplot(show.legend = FALSE) +
    theme_minimal() +
    labs(
        x = "Sample",
        y = "log10(count + 1)",
        title = "Raw count distributions (Kallisto)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Density plot of raw counts:

```r
ggplot(counts_melted, aes(x = log10_count, fill = sample)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(
        x = "log10(count + 1)",
        y = "Density",
        title = "Count density distributions"
    )
```

### Variance stabilizing transformation

Apply VST for exploratory analysis:

```r
vsd <- vst(dds, blind = TRUE)
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use VST?

Raw counts have a strong mean-variance relationship: highly expressed genes have higher variance. VST removes this dependency, making distance-based methods (PCA, clustering) more reliable.

Setting `blind = TRUE` ensures the transformation is not influenced by the experimental design, which is appropriate for quality control.

:::::::::::::::::::::::::::::::::::::::

### Sample-to-sample distances

```r
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    main = "Sample distance heatmap (Kallisto)"
)
```

The distance heatmap shows how similar samples are to each other. Samples from the same condition should cluster together.

### PCA plot

```r
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = name)) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw() +
    ggtitle("PCA of Kallisto-quantified samples")
```

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Interpret exploratory plots

Using the distance heatmap and PCA plot:

1. Do samples cluster by experimental condition?
2. Are there any outlier samples that don't group with their replicates?
3. How much variance is explained by PC1? What might this represent biologically?

::::::::::::::::::::::::::::::::::: solution

Example interpretation:

1. Samples should cluster by condition (WT_mock vs WT_IR) in both the heatmap and PCA.
2. If all replicates cluster together, there are no obvious outliers.
3. PC1 typically captures the largest source of variation. If it separates conditions, the experimental treatment is the dominant signal.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Step 4: Differential expression analysis

Run the full DESeq2 pipeline:

```r
dds <- DESeq(dds)
```

```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
```

Inspect dispersion estimates:

```r
plotDispEsts(dds)
```

::::::::::::::::::::::::::::::::::::::: callout

## Interpreting the dispersion plot

- **Black dots**: Gene-wise dispersion estimates.
- **Red line**: Fitted trend (shrinkage target).
- **Blue dots**: Final shrunken estimates.

A good fit shows the red line passing through the center of the black cloud, with blue dots closer to the line than the original black dots.

:::::::::::::::::::::::::::::::::::::::

Extract results for the contrast of interest:

```r
res <- results(
    dds,
    contrast = c("condition", "WT_IR", "WT_mock")
)

summary(res)
```

```
out of 21847 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3842, 18%
LFC < 0 (down)     : 3956, 18%
outliers [1]       : 0, 0%
low counts [2]     : 4182, 19%
(mean count < 5)
```

Order results by adjusted p-value:

```r
res_ordered <- res[order(res$padj), ]
head(res_ordered)
```

### Apply log fold change shrinkage

LFC shrinkage improves estimates for genes with low counts or high dispersion.

First, check the available coefficients:

```r
resultsNames(dds)
```

```
[1] "Intercept"                              "condition_WT_IR_vs_WT_mock"
```

Apply shrinkage using `apeglm` (recommended for standard contrasts):

```r
res_shrunk <- lfcShrink(
    dds,
    coef = "condition_WT_IR_vs_WT_mock",
    type = "apeglm"
)
```

::::::::::::::::::::::::::::::::::::::: callout

## Why shrink log fold changes?

Genes with low counts can have unreliably large fold changes. Shrinkage:

- Reduces noise in LFC estimates.
- Improves ranking for downstream analyses (e.g., GSEA).
- Does not affect p-values or significance calls.

## Choosing a shrinkage method

- **`apeglm`** (recommended): Fast, well-calibrated, uses coefficient name (`coef`).
- **`ashr`**: Required for complex contrasts that can't be specified with `coef` (e.g., interaction terms).
- **`normal`**: Legacy method, generally not recommended.

For simple two-group comparisons like ours, `apeglm` is preferred.

:::::::::::::::::::::::::::::::::::::::

## Step 5: Summarize and visualize results

Create a summary table:

```r
log2fc_cut <- log2(1.5)

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)

summary_table <- tibble(
    total_genes = nrow(res_df),
    sig = sum(res_df$padj < 0.05, na.rm = TRUE),
    up = sum(res_df$padj < 0.05 & res_df$log2FoldChange > log2fc_cut, na.rm = TRUE),
    down = sum(res_df$padj < 0.05 & res_df$log2FoldChange < -log2fc_cut, na.rm = TRUE)
)

print(summary_table)
```

### Volcano plot

```r
res_df <- res_df %>%
    mutate(
        sig = case_when(
            padj <= 0.05 & log2FoldChange >= log2fc_cut ~ "up",
            padj <= 0.05 & log2FoldChange <= -log2fc_cut ~ "down",
            TRUE ~ "ns"
        )
    )

ggplot(
    res_df,
    aes(
        x = log2FoldChange,
        y = -log10(padj),
        col = sig
    )
) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(
        "up" = "firebrick",
        "down" = "dodgerblue3",
        "ns" = "grey70"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", color = "grey40") +
    theme_classic() +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    ggtitle("Volcano plot: WT_IR vs WT_mock (Kallisto)")
```

::::::::::::::::::::::::::::::::::::::: callout

## Interpreting the volcano plot

- **X-axis**: Direction and magnitude of change (positive = higher in WT_IR, i.e., upregulated by radiation).
- **Y-axis**: Statistical significance (-log10 scale, higher = more significant).
- **Colored points**: Genes passing both fold change and significance thresholds.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Compare with genome-based results

If you also ran Episode 05 (genome-based workflow):

1. Are the number of DE genes similar between the two approaches?
2. Do the top DE genes overlap?
3. Are there systematic differences in fold change estimates?

::::::::::::::::::::::::::::::::::: solution

Expected observations:

1. The number of DE genes should be broadly similar, though not identical.
2. Most top DE genes should overlap between methods.
3. Kallisto may produce slightly different LFC estimates due to its different quantification algorithm.

Both approaches are valid; consistency between them increases confidence in the results.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Step 6: Save results

Save the full results table:

```r
write_tsv(
    res_df,
    "results/deseq2_kallisto/DESeq2_kallisto_results.tsv"
)
```

Save significant genes only:

```r
sig_res <- res_df %>%
    filter(
        padj <= 0.05,
        abs(log2FoldChange) >= log2fc_cut
    )

write_tsv(
    sig_res,
    "results/deseq2_kallisto/DESeq2_kallisto_sig.tsv"
)
```

Save the DESeq2 object for downstream analysis:

```r
saveRDS(dds, "results/deseq2_kallisto/dds_kallisto.rds")
```

::::::::::::::::::::::::::::::::::::::: discussion

## Genome-based vs. transcript-based: Which to choose?

| Aspect | Genome-based (Ep 05) | Transcript-based (Ep 05b) |
|--------|---------------------|--------------------------|
| Input | BAM files | FASTQ files |
| Speed | Slower (alignment + counting) | Faster (pseudo-alignment) |
| Storage | Large (BAM files) | Small (abundance.tsv files) |
| Bias correction | Limited | Sequence + GC bias |
| Novel transcripts | Can detect | Cannot detect |
| Visualization | IGV compatible | No BAM files |

For standard differential expression, both methods produce comparable results. Choose based on your specific needs and available resources.

:::::::::::::::::::::::::::::::::::::::

## Summary

::::::::::::::::::::::::::::::::::::: keypoints

- Salmon/Kallisto output is imported via `tximport` and loaded with `DESeqDataSetFromTximport()`.
- The tximport object preserves transcript length information for accurate normalization.
- Exploratory analysis (PCA, distance heatmaps) should precede differential expression testing.
- DESeq2 handles the statistical analysis identically to the genome-based workflow.
- LFC shrinkage improves fold change estimates for low-count genes.
- Results from transcript-based and genome-based workflows should be broadly concordant.

::::::::::::::::::::::::::::::::::::::::::::::::
