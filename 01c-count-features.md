---
source: Rmd
title: Quantifying Gene Expression
teaching: 20
exercises: 15
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do we quantify gene expression from aligned reads?
- What tools can we use to count reads mapped to genes or other genomic features?
- What options do we need to consider for accurate read counting?
- How do we interpret read count outputs for downstream analyses?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the process of read counting and its role in RNA-seq analysis.
- Learn how to use popular tools like `featureCounts` to count reads for genes or other genomic features.
- Recognize common challenges, such as handling overlapping features or isoforms, multimapped reads, and learn how to address them.
- Interpret read count outputs and prepare data for downstream analyses like differential expression.


::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction


In this episode, we will focus on counting the number of reads mapped to each genomic feature, such as genes or exons, a critical step in RNA-seq analysis. Accurate read counting provides the basis for downstream analyses, including differential gene expression. We will explore tools and methods used to assign reads to features and prepare the data for meaningful biological insights.

There are quite a few tools available for counting reads in RNA-seq analysis, each with its strengths and limitations.

:::::::::::::::::::::::::::::::::::::::  prereq

## What programs should I use for counting?

Some popular tools include:

- `featureCounts`: A fast and efficient tool for counting reads aligned to genomic features such as genes, exons, or transcripts.
- `HTSeq-count`: Widely used for counting reads mapped to genes or other genomic features in RNA-seq experiments.
- `RSEM`: For quantifying transcript-level expression using RNA-seq data, particularly for isoform-level counts.
- `Salmon` & `Kallisto`: Pseudo-alignment tools that provide transcript-level quantification directly from raw reads.
- `STAR`: Includes options like `--quantMode GeneCounts` for gene-level read counting during the alignment process.

:::::::::::::::::::::::::::::::::::::::

In this episode, we will focus on using `featureCounts`, a popular tool for counting reads mapped to genomic features. We will explore its features, parameters, and output formats to understand how read counting influences downstream analyses like differential expression.


## Organize your data

We will need `bam` files for read counting. These were generated in the previous episode, by mapping RNAseq reads to the genome. We will softlink them from the `mapping` directory to the `counting` directory before proceeding.

```bash
cd rcac_rnaseq
ln -s ../mapping/*.bam counts/
```


## Counting reads with `featureCounts`

`featureCounts` is a widely used tool for counting reads mapped to genomic features, such as genes, exons, or transcripts. It is part of the [Subread](https://subread.sourceforge.net/) package and provides a fast and efficient way to quantify gene expression from RNA-seq data. Genomic features are defined in a GTF (Gene Transfer Format) file, which contains the genomic coordinates of genes, exons, and other features. `bam` files provide the read alignments to the reference genome. Using the location information from the GTF file, `featureCounts` assigns read counts to features and generates count matrices for downstream analyses.

:::::::::::::::::::::::::::::::::::::::  prereq

## What type of features should I count reads for (genes, exons, transcripts)?

When counting RNA-seq reads, the most common feature to count is genes, as it provides gene-level expression data for downstream analysis. However, if your study focuses on alternative splicing or transcript isoforms, counting exons or transcripts might be more appropriate. The choice depends on the biological question: gene-level counts for overall expression, exon-level for splicing events, or transcript-level for isoform-specific expression.

:::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::  challenge

## What options do you need to specify for counting reads?

When using `featureCounts`, you need to specify several options to ensure accurate read counting. Take a look at the [featureCounts manual](https://subread.sourceforge.net/featureCounts.html) and identify the key parameters you need to provide for counting reads. There are also example commands specified in the manual that you can use as a reference. 


::::::::::::::::::::::::::::::::::: solution

The key parameters you need to specify for counting reads with `featureCounts` include:

- `-a`: Annotation file in GTF format containing genomic features.
- `-o`: Output file to store read counts.
- `-t`: Feature type to count (e.g., `gene`, `exon`, `transcript`).
- `-g`: Attribute type to group features (e.g., `gene_id`, `transcript_id`).
- `-T`: Number of threads to use for parallel processing.
- `-p`: Count read pairs instead of individual reads (for paired-end data).
- `--countReadPairs`: Count read pairs instead of individual reads (when using paired-end data).
- `--primary`: Count only primary alignments.


::::::::::::::::::::::::::::::::::: 

:::::::::::::::::::::::::::::::::::::::



Putting it all together, the `featureCounts` command to count reads for genes in a BAM file would look like this:

```bash
module load biocontainers
module load subread
featureCounts \
   -p \ # input is paired-end data
   --countReadPairs \ # count read pairs instead of individual reads
   -t exon \ # feature type to count
   -g gene_id \ # attribute type to group features
   -a data/annotation.gtf \ # annotation file in GTF format
   -o counts/counts.txt \ # output file for read counts
   counts/*.bam # input BAM files
```

You can do this as a job file and submit it to the cluster for processing. The output will be a tab-delimited file with read counts for each gene in each sample.


```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --job-name=count-features
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.out
ml biocontainers
ml subread
gtf="/home/aseethar/rcac_rnaseq/data/gencode.vM35.primary_assembly.basic.annotation.gtf"
featureCounts \
   -p \  
   --countReadPairs \
   -t exon \
   -T 16 \
   -g gene_id \
   -a ${gtf} \
   -o counts/counts.txt \
   counts/*.bam
```

You can submit this job file to the cluster using the `sbatch` command:

```bash
sbatch count_features.sh
```

This will run the `featureCounts` command on the cluster, counting reads for each gene in the BAM files provided.


## Understanding the output

The output of `featureCounts` is a tab-delimited file containing read counts for each feature in each sample. The first few lines of the output file look like this:

```bash
Geneid  Chr Start   End Strand  Length  Sample1 Sample2 Sample3
ENSMUSG00000000001  1   3073253 3074322 +   1069    0   0   0
ENSMUSG00000000003  1   3102016 3102125 -   110    0   0   0
ENSMUSG00000000028  1   3491927 3513553 -   21627   0   0   0
ENSMUSG00000000031  1   3461304 3513553 -   52250   0   0   0
ENSMUSG00000000037  1   1442067 1467620 +   25554   0   0   0
```

The columns in the output file include:

- `Geneid`: Gene identifier from the GTF file.
- `Chr`: Chromosome where the gene is located.
- `Start`: Start position of the gene.
- `End`: End position of the gene.
- `Strand`: Strand where the gene is located.
- `Length`: Length of the gene.
- `Sample1`, `Sample2`, `Sample3`, etc.: Read counts for each sample.


::::::::::::::::::::::::::::::::::::: keypoints 

- Quantifying gene expression involves counting the number of reads mapped to genomic features like genes or exons, which provides data for downstream analyses, such as differential gene expression.

- Tools like `featureCounts`, `RSEM`, `HTSeq-count`, `Salmon`, and `Kallisto` are commonly used for counting reads, each offering different advantages based on feature type and data structure.

- Depending on the research question, reads can be counted at various feature levels, such as genes, exons, or transcripts, to capture overall gene expression or more granular details like splicing events.

- The output of a read counting tool typically includes read counts for each feature across all samples, serving as input for further analysis, such as identifying differentially expressed genes.

::::::::::::::::::::::::::::::::::::::::::::::::

