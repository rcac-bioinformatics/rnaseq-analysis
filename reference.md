---
title: Reference
---


Here are selected method-specific references covering the major components of RNA-seq analysis, including alignment, quantification, quality control, differential expression, functional interpretation, and batch correction. These papers represent widely used tools and foundational methods in the field.


### General RNA-seq references

1. [A beginner’s guide to analysis of RNA sequencing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC6096346/) : K. M. Koch et al.
   A practical overview of the key steps in a typical RNA-seq analysis: library preparation, QC, alignment, counting, differential expression, and pitfalls.
2. [RNA sequencing (RNA-seq) methods & workflows](https://www.illumina.com/techniques/sequencing/rna-sequencing.html) : Illumina
   A vendor-agnostic guide describing different RNA-seq workflows, including bulk mRNA-seq, total RNA-seq, strand specificity, and exploratory considerations for experimental design.
3. [RNA sequencing data: Hitchhiker’s guide to expression, interpretation and beyond](https://www.annualreviews.org/content/journals/10.1146/annurev-biodatasci-072018-021255) : K. Van den Berge et al.
   A detailed review covering design, quantification, differential expression, and emerging types of RNA-seq (e.g., long-read, single-cell).


---

### Core differential expression methods

1. **DESeq2**
   Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 2014.
   [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

2. **edgeR**
   Robinson MD, McCarthy DJ, Smyth GK. *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.* Bioinformatics, 2010.
   [https://academic.oup.com/bioinformatics/article/26/1/139/182458](https://academic.oup.com/bioinformatics/article/26/1/139/182458)

3. **limma-voom**
   Law CW, Chen Y, Shi W, Smyth GK. *voom: precision weights unlock linear model analysis tools for RNA-seq read counts.* Genome Biology, 2014.
   [https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)

4. **Sleuth (for Salmon/Kallisto quantification)**
   Pimentel H, Bray NL, Puente S, Melsted P, Pachter L. *Differential analysis of RNA-seq incorporating quantification uncertainty.* Nature Methods, 2017.
   [https://www.nature.com/articles/nmeth.4324](https://www.nature.com/articles/nmeth.4324)

---


### Alignment methods

1. **STAR aligner**
   Dobin A et al. *STAR: ultrafast universal RNA-seq aligner.* Bioinformatics, 2013.
   [https://academic.oup.com/bioinformatics/article/29/1/15/272537](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

2. **HISAT2**
   Kim D, Langmead B, Salzberg SL. *HISAT: a fast spliced aligner with low memory requirements.* Nature Methods, 2015.
   [https://www.nature.com/articles/nmeth.3317](https://www.nature.com/articles/nmeth.3317)

3. **TopHat2 (historical, not recommended now but still cited)**
   Kim D et al. *TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions.* Genome Biology, 2013.
   [https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r36](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r36)


---

### Quantification methods

1. **Salmon (quasi-mapping and lightweight alignment)**
   Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. *Salmon provides fast and bias-aware quantification of transcript expression.* Nature Methods, 2017.
   [https://www.nature.com/articles/nmeth.4197](https://www.nature.com/articles/nmeth.4197)

2. **Kallisto (pseudo-alignment)**
   Bray NL, Pimentel H, Melsted P, Pachter L. *Near-optimal probabilistic RNA-seq quantification.* Nature Biotechnology, 2016.
   [https://www.nature.com/articles/nbt.3519](https://www.nature.com/articles/nbt.3519)

3. **RSEM (alignment-based quantification)**
   Li B, Dewey CN. *RSEM: accurate transcript quantification from RNA-seq data with or without a reference genome.* BMC Bioinformatics, 2011.
   [https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)

---


### Quality control and preprocessing

1. **FastQC**
   Andrews S. *FastQC: A quality control tool for high throughput sequence data.* Babraham Institute, 2010.
   [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

2. **Trim Galore** (cutadapt + FastQC wrapper)
   Krueger F. *Trim Galore.* Babraham Institute.
   [https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

3. **cutadapt**
   Martin M. *Cutadapt removes adapter sequences from high-throughput sequencing reads.* EMBnet Journal, 2011.
   [https://journal.embnet.org/index.php/embnetjournal/article/view/200](https://journal.embnet.org/index.php/embnetjournal/article/view/200)

4. **MultiQC**
   Ewels P, Magnusson M, Lundin S, Käller M. *MultiQC: summarize analysis results for multiple tools and samples in a single report.* Bioinformatics, 2016.
   [https://academic.oup.com/bioinformatics/article/32/19/3047/2196507](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)
---


### Functional analysis and interpretation

1. **GOseq (for GO analysis accounting for transcript length bias)**
   Young MD, Wakefield MJ, Smyth GK, Oshlack A. *Gene ontology analysis for RNA-seq: accounting for selection bias.* Genome Biology, 2010.
   [https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14)

2. **clusterProfiler (widely used for GO/KEGG)**
   Yu G, Wang LG, Han Y, He QY. *clusterProfiler: an R package for comparing biological themes among gene clusters.* OMICS, 2012.
   [https://www.liebertpub.com/doi/10.1089/omi.2011.0118](https://www.liebertpub.com/doi/10.1089/omi.2011.0118)

3. **fgsea (fast GSEA)**
   Korotkevich G, Sukhov V, Sergushichev A. *Fast gene set enrichment analysis.* bioRxiv, 2016.
   [https://www.biorxiv.org/content/10.1101/060012v3](https://www.biorxiv.org/content/10.1101/060012v3)

4. **KEGG**
   Kanehisa M, Goto S. *KEGG: Kyoto Encyclopedia of Genes and Genomes.* Nucleic Acids Research, 2000.
   [https://academic.oup.com/nar/article/28/1/27/2384398](https://academic.oup.com/nar/article/28/1/27/2384398)

---


### Splicing and isoform analysis (if you want to include optional advanced topics)

1. **DEXSeq (differential exon usage)**
   Anders S, Reyes A, Huber W. *Detecting differential usage of exons from RNA-seq data.* Genome Research, 2012.
   [https://genome.cshlp.org/content/22/10/2008](https://genome.cshlp.org/content/22/10/2008)

2. **SUPPA2 (isoform-level splicing changes)**
   Trincado JL et al. *SUPPA2 provides fast, accurate, and uncertainty-aware differential splicing analysis.* Genome Biology, 2018.
   [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1417-1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1417-1)

3. **StringTie2 (transcript assembly and quantification)**
   Kovaka S et al. *Transcriptome assembly from long-read RNA-seq alignments with StringTie2.* Genome Biology, 2019.
   [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1)

---


### Batch correction and confounder modeling

1. **ComBat**
   Johnson WE, Li C, Rabinovic A. *Adjusting batch effects in microarray expression data using empirical Bayes methods.* Biostatistics, 2007.
   [https://academic.oup.com/biostatistics/article/8/1/118/252073](https://academic.oup.com/biostatistics/article/8/1/118/252073)

   (Still widely used for RNA-seq after variance stabilizing transformation.)

2. **RUVSeq**
   Risso D, Ngai J, Speed TP, Dudoit S. *Normalization of RNA-seq data using factor analysis of control genes or samples.* Nature Biotechnology, 2014.
   [https://www.nature.com/articles/nbt.2931](https://www.nature.com/articles/nbt.2931)
---


## Glossary

### Experimental concepts

**adapter trimming**
Removal of sequencing adapters and low-quality bases from raw reads before alignment, typically using tools such as `fastp` or `Trim Galore`.

**alignment (read mapping)**
Assigning sequencing reads to a reference genome or transcriptome using aligners like HISAT2 or STAR.

**biological replicate**
Independent samples from the same biological condition, used to estimate natural biological variation.

**batch effects**
Confounding technical variation introduced by differences in library prep, sequencing runs, or operators; must be accounted for in experimental design and statistical models.

**bulk RNA-seq**
Conventional RNA-seq performed on pooled cells or tissues, as opposed to single-cell or spatial RNA-seq.

**cDNA**
Complementary DNA synthesized from RNA templates during library preparation.

**library preparation**
Experimental conversion of RNA into a sequencing-ready library, including fragmentation, adapter ligation, and amplification.

**paired-end sequencing**
Sequencing both ends of DNA/cDNA fragments to improve mapping accuracy and detection of splice junctions.

**read depth / sequencing depth**
The total number of reads generated per sample; higher depth increases power to detect lowly expressed or differentially expressed genes.

**stranded library / strand specificity**
Preservation of the transcriptional strand of origin during library prep, allowing direction-specific quantification.

---


### File formats and data organization

**FASTQ**
A text format containing raw sequencing reads with per-base quality scores.

**GTF / GFF**
Annotation file formats describing gene and transcript coordinates, essential for read counting and feature quantification.

**counts table**
A matrix of raw read counts (rows = genes or transcripts, columns = samples) used for downstream differential expression analysis.

---


### Quantification and normalization

**featureCounts**
A program that summarizes aligned reads to annotated genomic features, producing a counts table for DE analysis.

**gene-level vs. transcript-level quantification**
Gene-level quantification collapses isoforms into a single total count per gene; transcript-level quantification measures individual isoform expression, critical for splicing or isoform-specific analyses.

**normalization**
Adjusting for library size and composition so samples are comparable. Common methods include TMM (edgeR), median-ratio (DESeq2), and TPM for within-sample comparisons. Use **raw counts** with DESeq2/edgeR for statistical modeling—FPKM and TPM are not appropriate for DE testing.

**FPKM / TPM**
FPKM = Fragments Per Kilobase per Million; TPM = Transcripts Per Million. TPM normalizes read counts across transcript lengths and library size and is preferred for *within-sample* comparisons, but **not** for differential expression across samples.

**pseudo-alignment**
Lightweight transcript quantification without full alignment, used by tools such as Salmon or Kallisto.

**reference-based vs. de novo RNA-seq**
Reference-based workflows align reads to a known genome or transcriptome. De novo workflows assemble transcripts directly from reads when no reference exists.

**transcriptome assembly**
Reconstruction of transcripts from RNA-seq reads, either de novo or guided by a reference.

---


### Differential expression and statistics

**differential expression (DE)**
Statistical identification of genes or transcripts with significant changes in expression between experimental conditions.

**dispersion**
A parameter in negative-binomial count models representing variability beyond Poisson noise; estimated in DESeq2 and edgeR.

**log₂ fold change (LFC)**
A measure of change in expression between two conditions on a log₂ scale; values > 1 or < –1 typically indicate meaningful differences.

**FDR / adjusted p-value**
The false discovery rate controls for multiple testing (Benjamini–Hochberg correction). Genes are typically considered significant at FDR < 0.05.

**normalization factor**
Scaling coefficient applied during model fitting to account for library size or composition bias (internally computed by DESeq2/edgeR).


---

### Interpretation and downstream analysis

**functional enrichment analysis**
Identification of Gene Ontology (GO) terms, pathways, or functional categories overrepresented among differentially expressed genes.

**gene ontology (GO)**
A controlled vocabulary describing gene functions in three domains: biological process, molecular function, and cellular component.

**KEGG**
The Kyoto Encyclopedia of Genes and Genomes; maps genes to metabolic and signaling pathways.

**visualization**
Graphical summaries such as PCA plots, volcano plots, or heatmaps used to assess overall patterns and highlight DE genes.

**workflow reproducibility**
Ensuring analyses can be replicated by others using automated, version-controlled workflows (e.g., Snakemake, Nextflow, Apptainer).



