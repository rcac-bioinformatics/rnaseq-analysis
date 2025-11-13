---
title: Reference
---

## Additional Reading

For those interested in further exploring RNA-seq analysis and workflows, here are some recommended resources:

1. **[A beginner’s guide to analysis of RNA sequencing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC6096346/)** – K. M. Koch et al.
   A practical overview of the key steps in a typical RNA-seq analysis: library preparation, QC, alignment, counting, differential expression, and pitfalls. ([PMC][1])
2. **[RNA sequencing (RNA-seq) methods & workflows](https://www.illumina.com/techniques/sequencing/rna-sequencing.html)** – Illumina
   A vendor-agnostic guide describing different RNA-seq workflows, including bulk mRNA-seq, total RNA-seq, strand specificity, and exploratory considerations for experimental design. ([Illumina][2])
3. **[RNA sequencing data: Hitchhiker’s guide to expression, interpretation and beyond](https://www.annualreviews.org/content/journals/10.1146/annurev-biodatasci-072018-021255)** – K. Van den Berge et al.
   A detailed review covering design, quantification, differential expression, and emerging types of RNA-seq (e.g., long-read, single-cell). ([Annual Reviews][3])

## Glossary

**adapter trimming**
The process of removing sequencing-adapter sequences (and optionally low-quality bases) from raw reads prior to alignment.

**alignment (read mapping)**
Placing sequencing reads against a reference genome or transcriptome to determine their origin and count them.

**biological replicate**
Independent experimental samples from the same condition, used to estimate biological variability.

**bulk RNA-Seq**
Standard RNA-seq in which RNA is extracted from a population of cells or tissue, as opposed to single-cell or spatial methods.

**cDNA**
Complementary DNA synthesised from RNA (often via reverse transcription) for sequencing.

**counts table**
A matrix of read counts (rows = genes or transcripts, columns = samples) used for downstream differential expression analysis.

**differential expression (DE)**
Identifying genes or transcripts that change in expression levels between conditions, typically with statistical testing.

**FASTQ**
A file format containing sequencing reads and their associated quality scores.

**featureCounts**
A software tool for summarising reads mapped to genomic features (genes, exons) to produce count data.

**FPKM / TPM**
Normalization units for RNA-seq: FPKM (Fragments Per Kilobase of transcript per Million mapped reads), TPM (Transcripts Per Million). Note: TPM is preferred for comparisons *within* samples, but care is needed when comparing *between* samples. ([galaxyproject.org][4])

**library preparation**
The experimental steps that convert RNA into a ready-to-sequence format (e.g., reverse transcription, fragmentation, adapter ligation, amplification).

**paired-end sequencing**
Sequencing both ends of DNA/cDNA fragments, allowing improved mapping and detection of splice junctions or structural events.

**pseudo-alignment**
A lightweight read-mapping approach (e.g., with kallisto or Salmon) that assigns reads to transcripts without full alignment.

**quality control (QC)**
Diagnostic checks on raw or mapped sequencing data (e.g., base quality, GC content, duplication, adapter contamination) to identify potential problems early.

**reference-based RNA-Seq**
An RNA-seq workflow that aligns reads against a known reference genome or transcriptome. ([galaxyproject.org][4])

**replicate variability**
Biological or technical variability between replicates; important to model for robust differential expression analysis.

**single-cell RNA-Seq (scRNA-Seq)**
High-throughput sequencing of individual cells to assess transcriptomic heterogeneity at single-cell resolution.

**stranded library / strand specificity**
Library preparation preserving the information about which DNA strand the RNA originated from, enabling better transcript annotation and expression quantification. ([galaxyproject.org][4])

**transcriptome**
The complete set of RNA transcripts (coding and non-coding) present in a cell or tissue at a given time.

**transcriptome assembly**
Reconstructing RNA transcripts from reads, typically when no good reference genome exists (de novo), or optionally guided by a genome.

**upstream splice junction**
The boundary upstream of an intron/exon junction; identifying these is critical for correct transcript model assignment.

**workflow reproducibility**
Use of documented, automated pipelines (e.g., with Snakemake, Nextflow, containers) to ensure RNA-seq analyses can be repeated by others.


[1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6096346/?utm_source=chatgpt.com "A Beginner's Guide to Analysis of RNA Sequencing Data"
[2]: https://www.illumina.com/techniques/sequencing/rna-sequencing.html?utm_source=chatgpt.com "RNA Sequencing | RNA-Seq methods & workflows"
[3]: https://www.annualreviews.org/content/journals/10.1146/annurev-biodatasci-072018-021255?utm_source=chatgpt.com "RNA Sequencing Data: Hitchhiker's Guide to Expression ..."
[4]: https://galaxyproject.org/tutorials/rb-rnaseq/?utm_source=chatgpt.com "RNAseq - an introduction"
