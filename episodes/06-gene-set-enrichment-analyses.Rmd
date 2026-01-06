---
source: Rmd
title: "Gene set enrichment analysis"
teaching: 45
exercises: 50
---

:::::::::::::::::::::::::::::::::::::: questions

- What is over-representation analysis and how does it work statistically?
- How do we identify functional pathways and biological themes from DE genes?
- How do we run enrichment analysis for GO, KEGG, and MSigDB gene sets?
- How do we interpret and compare enrichment results across different databases?
- What are the limitations and best practices for enrichment analysis?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the statistical basis of over-representation analysis (ORA).
- Prepare gene lists and background universes from DESeq2 results.
- Run GO enrichment using clusterProfiler and interpret ontology categories.
- Run KEGG pathway enrichment with organism-specific data.
- Perform MSigDB Hallmark enrichment for curated biological signatures.
- Produce and interpret barplots, dotplots, and enrichment maps.
- Compare results across multiple gene set databases.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

Differential expression analysis identifies individual genes that change between conditions, but interpreting long gene lists is challenging. **Gene set enrichment analysis** provides biological context by asking: *Are genes with specific functions over-represented among the differentially expressed genes?*

This episode covers **over-representation analysis (ORA)**, the most common approach for interpreting DE results.

::::::::::::::::::::::::::::::::::::::: prereq

## What you need for this episode

- DE results from Episode 05 or 05b (`DESeq2_results_joined.tsv` or `DESeq2_salmon_results.tsv`)
- RStudio session via Open OnDemand
- Internet connection (for KEGG queries)

## Required packages

- clusterProfiler (Bioconductor)
- org.Mm.eg.db (mouse annotation)
- msigdbr (MSigDB gene sets)
- enrichplot (visualization)
- dplyr, ggplot2, readr

:::::::::::::::::::::::::::::::::::::::

## Background: Understanding over-representation analysis

### What is ORA?

ORA tests whether a predefined gene set contains more DE genes than expected by chance. It answers the question: *If I randomly selected the same number of genes from my background, how often would I see this many genes from pathway X?*

### The statistical test

ORA uses the **hypergeometric test** (equivalent to Fisher's exact test), which calculates the probability of observing *k* or more genes from a pathway given:

- **N**: Total genes in the background (universe)
- **K**: Genes in the pathway of interest
- **n**: Total DE genes in your list
- **k**: DE genes that are in the pathway

```
                    Pathway genes    Non-pathway genes
DE genes                 k              n - k
Non-DE genes           K - k          N - K - n + k
```

::::::::::::::::::::::::::::::::::::::: callout

## The hypergeometric distribution

The p-value represents the probability of drawing *k* or more pathway genes when randomly selecting *n* genes from a population of *N* genes containing *K* pathway members:

$$P(X \geq k) = \sum_{i=k}^{\min(K,n)} \frac{\binom{K}{i}\binom{N-K}{n-i}}{\binom{N}{n}}$$

This is a one-tailed test asking if the pathway is over-represented (enriched) in your gene list.

:::::::::::::::::::::::::::::::::::::::

### Multiple testing correction

With thousands of gene sets tested, many false positives would occur by chance. We apply **false discovery rate (FDR)** correction using the Benjamini-Hochberg method:

1. Rank all p-values from smallest to largest.
2. Calculate adjusted p-value: $p_{adj} = p \times \frac{n}{rank}$
3. Filter results by adjusted p-value threshold (typically 0.05).

::::::::::::::::::::::::::::::::::::::: callout

## Why the universe matters

The background universe should include **all genes that could have been detected as DE** in your experiment, typically all genes tested by DESeq2.

Common mistakes:
- Using the whole genome (inflates significance).
- Using only expressed genes (may miss relevant pathways).
- Using a different species' gene list.

The universe defines what "expected by chance" means. An inappropriate universe leads to misleading results.

:::::::::::::::::::::::::::::::::::::::

### Gene set databases

Different databases capture different aspects of biology:

| Database | Description | Best for |
|----------|-------------|----------|
| GO BP | Biological Process ontology | Functional processes |
| GO MF | Molecular Function ontology | Protein activities |
| GO CC | Cellular Component ontology | Subcellular localization |
| KEGG | Curated pathways | Metabolic/signaling pathways |
| Reactome | Detailed pathway reactions | Mechanistic understanding |
| MSigDB Hallmark | 50 curated signatures | High-level biological states |
| MSigDB C2 | Curated gene sets | Canonical pathways |

::::::::::::::::::::::::::::::::::::::: discussion

## Which database should I use?

There's no single "best" database. Consider:

- **GO**: Comprehensive but redundant; many overlapping terms.
- **KEGG**: Well-curated pathways; limited coverage; requires Entrez IDs.
- **Hallmark**: Reduced redundancy; easier to interpret; may miss specific pathways.

Best practice: Run multiple databases and look for **convergent signals** across them.

:::::::::::::::::::::::::::::::::::::::

## Step 1: Load DE results and prepare gene lists

Start your R session:

```bash
sinteractive -A rcac -p cpu -N 1 -n 4 --time=02:00:00
module load biocontainers
module load r-rnaseq
R
```

Load packages:

```r
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(readr)

setwd(paste0(Sys.getenv("SCRATCH"), "/rnaseq-workshop"))
```

Load DESeq2 results from Episode 05:

```r
res <- read_tsv("results/deseq2/DESeq2_results_joined.tsv", show_col_types = FALSE)
head(res)
```

Define the gene lists:

```r
# Universe: all genes tested
universe <- res$ensembl_gene_id_version

# Significant genes: padj < 0.05 and |log2FC| > log2(1.5)
log2fc_cut <- log2(1.5)
sig_genes <- res %>%
    filter(padj < 0.05, abs(log2FoldChange) > log2fc_cut) %>%
    pull(ensembl_gene_id_version)

length(universe)
length(sig_genes)
```

```
[1] 15839
[1] 4543
```

::::::::::::::::::::::::::::::::::::::: callout

## Choosing significance thresholds

The threshold for "significant" genes affects enrichment results:

- **Stringent (padj < 0.01, |LFC| > 1)**: Fewer genes, clearer signals, may miss subtle effects.
- **Lenient (padj < 0.1, no LFC cutoff)**: More genes, noisier results, higher sensitivity.

There's no universally correct threshold. Consider your experimental question and validate key findings.

:::::::::::::::::::::::::::::::::::::::

### Convert gene IDs

Most enrichment tools require **Entrez IDs**. Our data uses Ensembl IDs with version numbers, so we need to convert:

```r
# Remove version suffix for mapping
universe_clean <- gsub("\\.\\d+$", "", universe)
sig_clean <- gsub("\\.\\d+$", "", sig_genes)

# Map Ensembl to Entrez
id_map <- bitr(
    universe_clean,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
)

head(id_map)
```

```
           ENSEMBL ENTREZID
1 ENSMUSG00000051285   381290
2 ENSMUSG00000025900    14468
3 ENSMUSG00000033793    70675
4 ENSMUSG00000033774    18777
5 ENSMUSG00000025903    14459
6 ENSMUSG00000090031   628071
```

::::::::::::::::::::::::::::::::::::::: callout

## Gene ID conversion challenges

Not all genes map successfully:

- Some Ensembl IDs lack Entrez equivalents.
- ID mapping databases may be outdated.
- Multi-mapping (one Ensembl → multiple Entrez) can occur.

`bitr()` drops unmapped genes silently. Check how many genes were lost:

```r
message(sprintf("Mapped %d of %d genes (%.1f%%)",
    nrow(id_map), length(universe_clean),
    100 * nrow(id_map) / length(universe_clean)))
```

:::::::::::::::::::::::::::::::::::::::

Create final gene lists with Entrez IDs:

```r
# Universe with Entrez IDs
universe_entrez <- id_map$ENTREZID

# Significant genes with Entrez IDs
sig_entrez <- id_map %>%
    filter(ENSEMBL %in% sig_clean) %>%
    pull(ENTREZID)

length(sig_entrez)
```

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Prepare up-regulated and down-regulated gene lists

Create separate gene lists for:
1. Up-regulated genes (log2FC > 0 and padj < 0.05)
2. Down-regulated genes (log2FC < 0 and padj < 0.05)

Why might analyzing these separately be informative?

::::::::::::::::::::::::::::::::::: solution

```r
# Up-regulated genes
up_genes <- res %>%
    filter(padj < 0.05, log2FoldChange > log2fc_cut) %>%
    pull(ensembl_gene_id_version)
up_clean <- gsub("\\.\\d+$", "", up_genes)
up_entrez <- id_map %>%
    filter(ENSEMBL %in% up_clean) %>%
    pull(ENTREZID)

# Down-regulated genes
down_genes <- res %>%
    filter(padj < 0.05, log2FoldChange < -log2fc_cut) %>%
    pull(ensembl_gene_id_version)
down_clean <- gsub("\\.\\d+$", "", down_genes)
down_entrez <- id_map %>%
    filter(ENSEMBL %in% down_clean) %>%
    pull(ENTREZID)

length(up_entrez)
length(down_entrez)
```

Analyzing up and down separately can reveal:
- Different biological processes activated vs. suppressed.
- Clearer signals when combining masks opposing effects.
- Direction-specific pathway involvement.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Step 2: GO enrichment analysis

### Run GO Biological Process enrichment

```r
ego_bp <- enrichGO(
    gene = sig_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
)

head(ego_bp, 10)
```

::::::::::::::::::::::::::::::::::::::: callout

## Understanding enrichGO output

Key columns in the results:

- **ID**: GO term identifier (e.g., GO:0006915)
- **Description**: Human-readable term name
- **GeneRatio**: DE genes in term / total DE genes with annotation
- **BgRatio**: Term genes in universe / total annotated genes in universe
- **pvalue**: Raw hypergeometric p-value
- **p.adjust**: BH-corrected p-value
- **qvalue**: Alternative FDR estimate
- **geneID**: List of DE genes in this term
- **Count**: Number of DE genes in this term

:::::::::::::::::::::::::::::::::::::::

### Visualize GO results

Dotplot showing top enriched terms:

```r
dotplot(ego_bp, showCategory = 20) +
    ggtitle("GO Biological Process enrichment")
```

The dotplot shows:
- **X-axis**: Gene ratio (proportion of DE genes in the term)
- **Y-axis**: GO term description
- **Color**: Adjusted p-value
- **Size**: Number of genes

Barplot alternative:

```r
barplot(ego_bp, showCategory = 15) +
    ggtitle("GO BP enrichment - Top 15 terms")
```

### Reduce GO redundancy

GO terms are hierarchical and often redundant. Use `simplify()` to remove similar terms:

```r
ego_bp_simple <- simplify(
    ego_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
)

dotplot(ego_bp_simple, showCategory = 15) +
    ggtitle("GO BP enrichment (simplified)")
```

::::::::::::::::::::::::::::::::::::::: callout

## Why simplify GO results?

The GO hierarchy means related terms often appear together:
- "regulation of cell death" and "positive regulation of cell death"
- "response to stimulus" and "response to external stimulus"

`simplify()` removes redundant terms based on semantic similarity, making results easier to interpret. The `cutoff` parameter (0-1) controls how aggressively terms are merged.

:::::::::::::::::::::::::::::::::::::::

### Enrichment map visualization

For many significant terms, an enrichment map shows relationships:

```r
ego_bp_em <- pairwise_termsim(ego_bp_simple)
emapplot(ego_bp_em, showCategory = 30) +
    ggtitle("GO BP enrichment map")
```

The enrichment map clusters similar terms together, revealing functional themes.

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Compare GO ontologies

Run enrichment for GO Molecular Function (MF) and Cellular Component (CC):

```r
ego_mf <- enrichGO(gene = sig_entrez, universe = universe_entrez,
                   OrgDb = org.Mm.eg.db, ont = "MF", readable = TRUE)
ego_cc <- enrichGO(gene = sig_entrez, universe = universe_entrez,
                   OrgDb = org.Mm.eg.db, ont = "CC", readable = TRUE)
```

1. How do the top terms differ between BP, MF, and CC?
2. Are there consistent themes across ontologies?
3. Which ontology provides the most interpretable results for your experiment?

::::::::::::::::::::::::::::::::::: solution

Example interpretation:

- **BP** terms describe processes (e.g., "cell cycle", "RNA processing").
- **MF** terms describe molecular activities (e.g., "RNA binding", "kinase activity").
- **CC** terms describe locations (e.g., "nucleus", "mitochondrion").

Consistent themes across ontologies strengthen biological interpretation. For example, if BP shows "translation" enrichment, MF might show "ribosome binding", and CC might show "ribosome" or "cytoplasm".

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Step 3: KEGG pathway enrichment

KEGG provides curated pathway maps connecting genes to biological systems.

```r
ekegg <- enrichKEGG(
    gene = sig_entrez,
    universe = universe_entrez,
    organism = "mmu",  # Mouse
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)

head(ekegg)
```

::::::::::::::::::::::::::::::::::::::: callout

## KEGG organism codes

Common organism codes:
- `"hsa"`: Human
- `"mmu"`: Mouse
- `"rno"`: Rat
- `"dme"`: Drosophila
- `"sce"`: Yeast
- `"ath"`: Arabidopsis

Find your organism: `search_kegg_organism("species name")`

:::::::::::::::::::::::::::::::::::::::

Visualize KEGG results:

```r
barplot(ekegg, showCategory = 15) +
    ggtitle("KEGG pathway enrichment")
```

::::::::::::::::::::::::::::::::::::::: callout

## KEGG limitations

Be aware of KEGG-specific issues:

1. **ID coverage**: KEGG uses its own gene IDs; some Entrez IDs may not map.
2. **Update frequency**: KEGG pathways are updated periodically; results may vary.
3. **License**: KEGG has usage restrictions for commercial applications.
4. **Incomplete mapping**: Some genes lack pathway assignments.

Check dropped genes:

```r
message(sprintf("KEGG analysis used %d of %d input genes",
    length(ekegg@gene), length(sig_entrez)))
```

:::::::::::::::::::::::::::::::::::::::

## Step 4: MSigDB Hallmark enrichment

MSigDB Hallmark gene sets are 50 curated signatures representing well-defined biological states and processes.

### Load Hallmark gene sets

```r
hallmark <- msigdbr(species = "Mus musculus", category = "H")
head(hallmark)
```

Prepare gene set list:

```r
hallmark_list <- hallmark %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(term = gs_name, gene = entrez_gene)
```

### Run Hallmark enrichment

```r
ehall <- enricher(
    gene = sig_entrez,
    universe = universe_entrez,
    TERM2GENE = hallmark_list,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)

head(ehall)
```

Visualize:

```r
dotplot(ehall, showCategory = 20) +
    ggtitle("MSigDB Hallmark enrichment")
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use Hallmark gene sets?

Advantages of Hallmark over GO/KEGG:

1. **Reduced redundancy**: Each signature is distinct and non-overlapping.
2. **Expert curation**: Refined from multiple sources to represent coherent biology.
3. **Interpretability**: 50 sets is manageable for detailed examination.
4. **Reproducibility**: Well-documented and stable definitions.

Common Hallmark sets include:
- HALLMARK_INFLAMMATORY_RESPONSE
- HALLMARK_P53_PATHWAY
- HALLMARK_OXIDATIVE_PHOSPHORYLATION
- HALLMARK_MYC_TARGETS_V1

:::::::::::::::::::::::::::::::::::::::

## Step 5: Compare results across databases

Create a summary of top enriched terms from each database:

```r
# Extract top 10 from each
top_go <- head(ego_bp_simple@result, 10) %>%
    mutate(Database = "GO BP") %>%
    dplyr::select(Database, Description, p.adjust, Count)

top_kegg <- head(ekegg@result, 10) %>%
    mutate(Database = "KEGG") %>%
    dplyr::select(Database, Description, p.adjust, Count)

top_hall <- head(ehall@result, 10) %>%
    mutate(Database = "Hallmark") %>%
    dplyr::select(Database, ID, p.adjust, Count) %>%
    dplyr::rename(Description = ID)

# Combine
comparison <- bind_rows(top_go, top_kegg, top_hall)
print(comparison, n = 30)
```

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Identify convergent themes

Looking at your GO, KEGG, and Hallmark results:

1. Are there biological themes that appear across multiple databases?
2. Are there database-specific findings that don't replicate elsewhere?
3. How would you prioritize which findings to investigate further?

::::::::::::::::::::::::::::::::::: solution

Convergent signals to look for:

- If GO shows "cell cycle" terms, KEGG might show "Cell cycle" pathway, and Hallmark might show "HALLMARK_G2M_CHECKPOINT".
- If GO shows "immune response", KEGG might show "Cytokine-cytokine receptor interaction", and Hallmark might show "HALLMARK_INFLAMMATORY_RESPONSE".

Database-specific findings may indicate:
- Pathway coverage differences.
- Annotation biases.
- Potentially spurious signals (require validation).

Prioritization strategy:
1. Focus on themes appearing in multiple databases.
2. Consider effect size (how many genes, fold enrichment).
3. Assess biological plausibility given experimental design.
4. Validate key findings with independent methods.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Step 5b: Direction-aware enrichment analysis

ORA on all significant genes (up + down combined) can mask biological signals when up-regulated and down-regulated genes have different functions. Analyzing directions separately often reveals clearer patterns.

### Prepare direction-specific gene lists

```r
# Up-regulated genes
up_genes <- res %>%
    filter(padj < 0.05, log2FoldChange > log2fc_cut) %>%
    pull(ensembl_gene_id_version)
up_clean <- gsub("\\.\\d+$", "", up_genes)
up_entrez <- id_map %>%
    filter(ENSEMBL %in% up_clean) %>%
    pull(ENTREZID)

# Down-regulated genes
down_genes <- res %>%
    filter(padj < 0.05, log2FoldChange < -log2fc_cut) %>%
    pull(ensembl_gene_id_version)
down_clean <- gsub("\\.\\d+$", "", down_genes)
down_entrez <- id_map %>%
    filter(ENSEMBL %in% down_clean) %>%
    pull(ENTREZID)

message(sprintf("Up-regulated: %d genes, Down-regulated: %d genes",
    length(up_entrez), length(down_entrez)))
```

### Run enrichment for each direction

```r
# GO enrichment for up-regulated genes
ego_up <- enrichGO(
    gene = up_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    readable = TRUE
)

# GO enrichment for down-regulated genes
ego_down <- enrichGO(
    gene = down_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    readable = TRUE
)
```

### Compare directions

```r
# Combine top results for comparison
up_top <- head(ego_up@result, 10) %>%
    mutate(Direction = "Up-regulated") %>%
    dplyr::select(Direction, Description, p.adjust, Count)

down_top <- head(ego_down@result, 10) %>%
    mutate(Direction = "Down-regulated") %>%
    dplyr::select(Direction, Description, p.adjust, Count)

direction_comparison <- bind_rows(up_top, down_top)
print(direction_comparison, n = 20)
```

::::::::::::::::::::::::::::::::::::::: callout

## When to analyze directions separately

Separate analysis is particularly useful when:

- You expect different biological processes to be activated vs. suppressed.
- Combined analysis shows few significant results (opposing signals may cancel out).
- You want to understand mechanism (what's gained vs. lost).

Separate analysis is less useful when:

- You have very few DE genes in one direction.
- The biological question is about overall pathway disruption.

:::::::::::::::::::::::::::::::::::::::

## Step 5c: Gene Set Enrichment Analysis (GSEA)

ORA has a fundamental limitation: it requires an arbitrary cutoff to define "significant" genes. **GSEA** uses the full ranked gene list without cutoffs.

### How GSEA differs from ORA

| Aspect | ORA | GSEA |
|--------|-----|------|
| Input | List of significant genes | All genes, ranked by fold change |
| Cutoff | Required (padj, LFC) | None |
| Sensitivity | May miss coordinated small changes | Detects subtle but consistent shifts |
| Gene weights | All DE genes treated equally | Magnitude of change matters |

### Prepare ranked gene list for GSEA

GSEA requires a named vector of gene scores (typically log2 fold change or a signed significance score):

```r
# Create ranked list using signed -log10(pvalue) * sign(LFC)
# This weights by both significance and direction
res_gsea <- res %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
        ensembl_clean = gsub("\\.\\d+$", "", ensembl_gene_id_version),
        rank_score = -log10(pvalue) * sign(log2FoldChange)
    )

# Map to Entrez IDs
res_gsea <- res_gsea %>%
    inner_join(id_map, by = c("ensembl_clean" = "ENSEMBL"))

# Create named vector (required format for GSEA)
gene_list <- res_gsea$rank_score
names(gene_list) <- res_gsea$ENTREZID

# Sort by rank score (required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)

head(gene_list)
tail(gene_list)
```

### Run GSEA with clusterProfiler

```r
gsea_go <- gseGO(
    geneList = gene_list,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)

head(gsea_go)
```

### Visualize GSEA results

```r
# Dot plot of enriched gene sets
dotplot(gsea_go, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GSEA: GO Biological Process")
```

The `.sign` column indicates whether the gene set is enriched in up-regulated (`activated`) or down-regulated (`suppressed`) genes.

### GSEA enrichment plot for a specific term

```r
# Enrichment plot for top gene set
gseaplot2(gsea_go, geneSetID = 1, title = gsea_go$Description[1])
```

::::::::::::::::::::::::::::::::::::::: callout

## Interpreting GSEA results

Key metrics:

- **NES (Normalized Enrichment Score)**: Direction and strength of enrichment. Positive = enriched in up-regulated genes; Negative = enriched in down-regulated genes.
- **pvalue/p.adjust**: Statistical significance.
- **core_enrichment**: Genes contributing most to the enrichment signal ("leading edge").

GSEA is particularly valuable when:

- Changes are subtle but coordinated across many genes.
- You want to avoid arbitrary significance cutoffs.
- You're comparing with published gene signatures.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Compare ORA and GSEA results

1. Do the same pathways appear in both ORA and GSEA results?
2. Does GSEA identify any pathways that ORA missed?
3. For pathways appearing in both, is the direction (NES sign) consistent with the ORA direction analysis?

::::::::::::::::::::::::::::::::::::::: solution

Expected observations:

- Many top pathways should overlap between ORA and GSEA.
- GSEA may identify additional pathways where changes are subtle but coordinated.
- The NES sign (positive/negative) should match whether the pathway was enriched in up- or down-regulated genes in ORA.
- Discrepancies may indicate pathways where the signal comes from many small changes rather than a few large ones.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Step 6: Save results

Save all enrichment results:

```r
# Create output directory
dir.create("results/enrichment", showWarnings = FALSE)

# GO results
write_csv(as.data.frame(ego_bp), "results/enrichment/GO_BP_enrichment.csv")
write_csv(as.data.frame(ego_bp_simple), "results/enrichment/GO_BP_simplified.csv")

# KEGG results
write_csv(as.data.frame(ekegg), "results/enrichment/KEGG_enrichment.csv")

# Hallmark results
write_csv(as.data.frame(ehall), "results/enrichment/Hallmark_enrichment.csv")

# Direction-aware results
write_csv(as.data.frame(ego_up), "results/enrichment/GO_BP_up_regulated.csv")
write_csv(as.data.frame(ego_down), "results/enrichment/GO_BP_down_regulated.csv")

# GSEA results
write_csv(as.data.frame(gsea_go), "results/enrichment/GSEA_GO_BP.csv")
```

Save plots:

```r
# GO dotplot
pdf("results/enrichment/GO_BP_dotplot.pdf", width = 10, height = 8)
dotplot(ego_bp_simple, showCategory = 20) +
    ggtitle("GO Biological Process enrichment")
dev.off()

# KEGG barplot
pdf("results/enrichment/KEGG_barplot.pdf", width = 10, height = 6)
barplot(ekegg, showCategory = 15) +
    ggtitle("KEGG pathway enrichment")
dev.off()

# Hallmark dotplot
pdf("results/enrichment/Hallmark_dotplot.pdf", width = 10, height = 8)
dotplot(ehall, showCategory = 20) +
    ggtitle("MSigDB Hallmark enrichment")
dev.off()
```

## Interpretation guidelines

### Signs of meaningful enrichment

Enrichment results are more reliable when:

- **Multiple related terms** appear (not just one isolated hit).
- **Consistent themes** emerge across databases.
- **Biological plausibility**: Results align with known biology of your system.
- **Reasonable gene counts**: Terms with many genes (>5-10) are more robust.
- **Moderate p-values**: Extremely low p-values (e.g., 10^-50) may indicate annotation bias.

### Common pitfalls

1. **Over-interpretation**: A significant p-value doesn't prove causation.
2. **Ignoring effect size**: A term with 3/1000 genes is less meaningful than 30/100.
3. **Cherry-picking**: Report all significant results, not just expected ones.
4. **Universe mismatch**: Using wrong background inflates or deflates significance.
5. **Outdated annotations**: Gene function annotations improve over time.

### Reporting enrichment results

In publications, include:

- Method and software versions used.
- Background universe definition.
- Multiple testing correction method.
- Full results tables (supplementary).
- Number of input genes and those successfully mapped.

## Summary

::::::::::::::::::::::::::::::::::::: keypoints

- ORA tests if gene sets contain more DE genes than expected by chance using the hypergeometric test.
- The background universe must include all genes that could have been detected as DE.
- Multiple testing correction (FDR) is essential when testing thousands of gene sets.
- GO provides comprehensive but redundant functional annotation; use `simplify()` to reduce redundancy.
- KEGG provides curated pathway maps but has limited gene coverage and requires Entrez IDs.
- MSigDB Hallmark sets offer 50 high-quality, non-redundant biological signatures.
- **Direction-aware analysis** (up vs. down separately) often reveals clearer biological signals.
- **GSEA** avoids arbitrary cutoffs by using the full ranked gene list; it detects subtle but coordinated changes.
- Convergent findings across multiple databases and methods strengthen biological interpretation.
- Always report methods, thresholds, and full results for reproducibility.

::::::::::::::::::::::::::::::::::::::::::::::::
