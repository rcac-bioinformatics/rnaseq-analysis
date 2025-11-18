---
source: Rmd
title: "Gene set enrichment analysis"
teaching: 35
exercises: 45
---

:::::::::::::::::::::::::::::::::::::: questions

- How do we identify functional pathways and biological themes from DE genes.
- How do we run over representation analysis for GO, KEGG, and MSigDB sets.
- How do we create simple visualizations of enriched terms.
- How do we interpret enrichment results in a biologically meaningful way.

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Prepare a gene list from DESeq2 results.
- Run GO enrichment using clusterProfiler.
- Run KEGG pathway enrichment with organism specific data.
- Perform MSigDB Hallmark enrichment.
- Produce barplots and dotplots for enriched terms.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

Differential expression gives a list of up regulated and down regulated genes.  
To understand the underlying biology, we examine which **functional gene sets** are enriched.  
This process is called **over representation analysis** (ORA).  
ORA identifies pathways or gene sets that contain more DE genes than expected by chance.

Here we will use **clusterProfiler**, **org.Mm.eg.db**, and **msigdbr** to run GO, KEGG, and Hallmark analyses on mouse data.

::::::::::::::::::::::::::::::::::::::: callout

## What ORA requires

- A background universe of all tested genes (usually all genes in the DESeq2 dataset).
- A subset of genes considered significant (for example padj < 0.05).
- A gene identifier type that matches the annotation database (Ensembl or Entrez).

:::::::::::::::::::::::::::::::::::::::

## Step 1: Load DE results and prepare gene lists

Start your R session:

```bash
sinteractive -A rcac -p cpu -N 1 -n 4 --time=02:00:00
module load biocontainers
module load r-rnaseq
R
```

Load data:

```r
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(readr)
```

Load DESeq2 results from Episode C:

```r
res <- read_csv("DESeq2_results_shrunken.csv", show_col_types = FALSE)
```

Define:

- universe = all genes tested  
- significant set = padj < 0.05  

```r
universe <- res$X1
sig <- res %>% filter(padj < 0.05) %>% pull(X1)
```

Convert Ensembl IDs to Entrez (required for GO and KEGG):

```r
id_map <- bitr(universe,
               fromType = "ENSEMBL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)
```

Merge IDs into results:

```r
res2 <- left_join(res, id_map, by = c("X1" = "ENSEMBL"))
sig_entrez <- res2 %>% filter(padj < 0.05) %>% pull(ENTREZID)
universe_entrez <- res2$ENTREZID
```

## Step 2: GO enrichment

Run GO Biological Process ORA:

```r
ego <- enrichGO(gene = sig_entrez,
                universe = universe_entrez,
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                readable = TRUE)
```

View top terms:

```r
head(ego)
```

Plot:

```r
dotplot(ego, showCategory = 20) +
    ggtitle("GO Biological Process enrichment")
```

::::::::::::::::::::::::::::::::::::::: discussion

## Why GO is useful

GO BP terms provide high level summaries of biological activity.  
Strongly biased terms (for example mitochondrial or translation) often correspond to global shifts in cell state or stress.

:::::::::::::::::::::::::::::::::::::::

## Step 3: KEGG enrichment

KEGG requires Entrez IDs and organism code `"mmu"` for mouse.

```r
ekegg <- enrichKEGG(gene = sig_entrez,
                    universe = universe_entrez,
                    organism = "mmu")
```

Plot KEGG enrichment:

```r
barplot(ekegg, showCategory = 15)
```

::::::::::::::::::::::::::::::::::::::: callout

## KEGG caveat

KEGG Entrez annotations may not perfectly match recent Ensembl builds.  
If a gene is missing, it will be silently dropped.

:::::::::::::::::::::::::::::::::::::::

## Step 4: MSigDB Hallmark enrichment

Load Hallmark gene sets for mouse:

```r
hallmark <- msigdbr(species = "Mus musculus", category = "H")
hallmark_list <- split(hallmark$entrez_gene, hallmark$gs_name)
```

Run ORA:

```r
ehall <- enricher(gene = sig_entrez,
                  universe = universe_entrez,
                  TERM2GENE = hallmark_list)
```

Plot:

```r
dotplot(ehall, showCategory = 20) +
    ggtitle("MSigDB Hallmark enrichment")
```

### Inspect enriched Hallmark sets

```r
head(ehall)
```

::::::::::::::::::::::::::::::::::::::: callout

## Why Hallmark sets are recommended

Hallmark pathways are curated to reduce redundancy.  
They are easier to interpret compared to long lists of GO terms.

:::::::::::::::::::::::::::::::::::::::

## Step 5: Saving results

```r
write_csv(as.data.frame(ego), "GO_enrichment.csv")
write_csv(as.data.frame(ekegg), "KEGG_enrichment.csv")
write_csv(as.data.frame(ehall), "Hallmark_enrichment.csv")
```

## Step 6: Interpretation guidelines

Meaningful enrichment patterns should:

- involve pathways relevant to your experimental design  
- contain multiple DE genes in the same direction  
- show consistent themes across GO, KEGG, and Hallmark sets  

Pathways enriched only by one DE gene or very small sets should be interpreted cautiously.

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: interpret enrichment results

Using your GO, KEGG, and Hallmark outputs:

1. Identify the pathways with the smallest adjusted p values.  
2. Are the enriched pathways consistent with the expected phenotype.  
3. Which pathway collections (GO, KEGG, Hallmark) give the clearest signal.  

::::::::::::::::::::::::::::::::::: solution

Example interpretation:

- Several GO terms related to RNA processing and cell cycle appear enriched.  
- Hallmark sets provide clearer summaries of biological themes compared to GO.  
- KEGG pathways show fewer hits but corroborate the GO findings.  

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Summary

::::::::::::::::::::::::::::::::::::: keypoints

- ORA requires a universe and a significant gene set.  
- Ensembl IDs must be converted to Entrez for GO and KEGG.  
- clusterProfiler provides GO BP and KEGG ORA.  
- Hallmark sets offer high quality curated pathways for interpretation.  
- Dotplots and barplots are effective for visualizing enrichment results.

::::::::::::::::::::::::::::::::::::::::::::::::
