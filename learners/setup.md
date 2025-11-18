---
title: Setup
---

## Instructors

1. **Arun Seetharam, Ph.D.**: Arun is a lead bioinformatics scientist at Purdue University’s Rosen Center for Advanced Computing. With extensive expertise in comparative genomics, genome assembly, annotation, single-cell genomics,  NGS data analysis, metagenomics, proteomics, and metabolomics. Arun supports a diverse range of bioinformatics projects across various organisms, including human model systems.

2. **Michael Carlson, Ph.D.**: Michael has a background in computational physics, specifically hypersonic materials. He also leads many introductory workshops in the High-Performance Computing domain.


## Schedule

| **Time**  | **Session**  |
|:---|-------------|
| **8:30 AM** | Arrival & Setup |
| **9:00 AM** | **Introduction to RNA-seq Analysis:** Overview of experimental design, biological replicates, sequencing depth, and workflow overview (QC → alignment → quantification → differential expression → visualization) |
| **9:45 AM** | **Data Preparation & Quality Control:** Exploring raw FASTQ files, running FastQC and MultiQC, trimming low-quality reads with fastp |
| **10:30 AM** | **Break** |
| **10:45 AM** | **Read Alignment & Quantification:**  Mapping reads with STAR, building genome indices, generating gene-level counts with featureCounts and Salmon |
| **12:00 PM** | **Lunch Break** |
| **1:00 PM** | **Differential Expression Analysis (DESeq2):** Importing count data into R, data normalization, exploring variance, identifying significantly expressed genes |
| **2:15 PM** | **Break** |
| **2:30 PM** | **Functional Interpretation & Visualization:** Gene ontology enrichment, pathway analysis, and visualizing results (volcano plots, heatmaps, PCA) |
| **3:30 PM** | **Integration & Reproducibility:** Creating reproducible pipelines with Snakemake/Nextflow templates, managing metadata and documentation |
| **4:00 PM** | **Wrap-Up & Discussion:** Troubleshooting, best practices, and Q&A |


## What is not covered

1. Raw data generation, library preparation, or experimental design optimization
2. De novo transcriptome assembly (e.g., Trinity) or genome-guided transcript reconstruction
3. Single-cell RNA-seq or spatial transcriptomics analysis
4. Alternative splicing, isoform quantification, or long-read transcript analysis
5. Advanced visualization dashboards or interactive analysis tools (e.g., Shiny, iDEP)

---

## Pre-requisites

1. Basic understanding of genomics concepts (genes, transcripts, and genome structure)
2. Familiarity with the command line interface (Linux/Unix shell)
3. Prior exposure to basic bioinformatics tools and file formats (FASTA, GFF, FASTQ)

---

## Data sets

To copy only the training data:

```bash
rsync -avP /depot/workshop/data/rnaseq_workshop ${RCAC_SCRATCH}/
```

A completed version of the workshop data is available at:

```
/depot/workshop/data/rnaseq_workshop-results
```

You can copy it to your scratch space using:

```bash
rsync -avP /depot/workshop/data/rnaseq_workshop-results ${RCAC_SCRATCH}/
```

Use this folder **only if you are unable to complete the exercises during the workshop**.
You will only need one directory on the Gilbreth cluster. See below for details.


---

## Software setup

::::::::::::::::::::::::::::::::::::::: discussion

### Details

SSH key setup for different systems is provided in the expandable sections below.
Follow the instructions for your operating system to configure passwordless access.

:::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::: solution

### Windows

Open **PowerShell** or **Git Bash** and run:

```bash
ssh-keygen -b 4096 -t rsa
type .ssh\id_rsa.pub | ssh trainXX@scholar.rcac.purdue.edu "mkdir -p ~/.ssh; cat >> ~/.ssh/authorized_keys"
```

:::::::::::::::::::::::::

:::::::::::::::: solution

### macOS

Open **Terminal** and run:

```bash
ssh-keygen -b 4096 -t rsa
cat .ssh/id_rsa.pub | ssh trainXX@scholar.rcac.purdue.edu "mkdir -p ~/.ssh; cat >> ~/.ssh/authorized_keys"
```

:::::::::::::::::::::::::

:::::::::::::::: solution

### Linux

Open a terminal and run:

```bash
ssh-keygen -b 4096 -t rsa
cat .ssh/id_rsa.pub | ssh trainXX@scholar.rcac.purdue.edu "mkdir -p ~/.ssh; cat >> ~/.ssh/authorized_keys"
```

:::::::::::::::::::::::::
