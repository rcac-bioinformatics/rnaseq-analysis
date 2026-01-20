---
title: Setup
---

## Instructors

1. **Arun Seetharam, Ph.D.**: Arun is a lead bioinformatics scientist at Purdue University’s Rosen Center for Advanced Computing. With extensive expertise in comparative genomics, genome assembly, annotation, single-cell genomics,  NGS data analysis, metagenomics, proteomics, and metabolomics. Arun supports a diverse range of bioinformatics projects across various organisms, including human model systems.

2. **Tomas Ratkus, Ph.D.**: is a Research Solutions Engineer with Rosen Center for Advanced Computing (RCAC), and a Senior Laboratory Operations Specialist in the Department of Earth, Atmospheric, and Planetary Sciences (EAPS). He primarily supports earth system and climate modeling researchers through high performance computing, software environments, and workflow support.

3. **Rose Wilfong, Ph.D.**: works jointly with the Bindley Bioscience Center on campus. Prior to joining, she did graduate research in bioinformatics, cheminformatics, and cryo-ET.


## Schedule (01/22/2026)


| **Time**     | **Session**                                                                                                                                                                                          |
|:---|-------------|
| **8:30 AM**  | Arrival & Setup                                                                                                                                                                                      |
| **9:00 AM**  | **Introduction to RNA-seq Analysis:** Experimental design, biological replicates, sequencing depth, and overview of the analysis workflow (QC → alignment → quantification → DE)                     |
| **9:45 AM**  | **Data Preparation & Quality Control:** Inspecting raw FASTQ files, running FastQC and MultiQC, trimming with fastp                                                                                  |
| **10:30 AM** | **Break**                                                                                                                                                                                            |
| **10:45 AM** | **Read Alignment & Quantification:** Mapping with STAR, building indices, generating gene-level counts with featureCounts (Kallisto covered conceptually only)                                         |
| **12:00 PM** | **Lunch Break**                                                                                                                                                                                      |
| **1:00 PM**  | **Differential Expression Analysis (DESeq2):** Importing counts into R, normalization, exploratory plots (VST, distance heatmaps, PCA), and identifying significantly differentially expressed genes |
| **2:15 PM**  | **Break**                                                                                                                                                                                            |
| **2:30 PM**  | **Visualization & Interpretation:** Volcano plots, heatmaps, PCA review, summary tables. Introduction to gene set enrichment methods (ORA/GSEA) with pointers to explore independently.              |
| **3:30 PM**  | **Wrap-Up & Discussion:** Review of workflow, troubleshooting common issues, recommended next steps                                                                                                  |
| **4:00 PM**  | End of Workshop                                                                                                                                                                                      |


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
rsync -avP /scratch/negishi/aseethar/rnaseq-workshop ${RCAC_SCRATCH}/
```

A completed version of the workshop data is available at:

```
/depot/workshop/data/rnaseq-workshop_results
```

You can copy it to your scratch space using:

```bash
rsync -avP /scratch/negishi/aseethar/rnaseq-workshop/rnaseq-workshop_results ${RCAC_SCRATCH}/
```

Use this folder **only if you are unable to complete the exercises during the workshop**.



---

## Software setup

::::::::::::::::::::::::::::::::::::::: discussion

## Details


SSH key setup for different systems is provided in the expandable sections below.
Follow the instructions for your operating system to configure passwordless access.

:::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::: solution

### Windows

Open **PowerShell** or **Git Bash** and run:

```bash
ssh-keygen -b 4096 -t rsa
type .ssh\id_rsa.pub | ssh boiler@scholar.rcac.purdue.edu "mkdir -p ~/.ssh; cat >> ~/.ssh/authorized_keys"
```

:::::::::::::::::::::::::

:::::::::::::::: solution

### macOS

Open **Terminal** and run:

```bash
ssh-keygen -b 4096 -t rsa
cat .ssh/id_rsa.pub | ssh boiler@scholar.rcac.purdue.edu "mkdir -p ~/.ssh; cat >> ~/.ssh/authorized_keys"
```

:::::::::::::::::::::::::

:::::::::::::::: solution

### Linux

Open a terminal and run:

```bash
ssh-keygen -b 4096 -t rsa
cat .ssh/id_rsa.pub | ssh boiler@scholar.rcac.purdue.edu "mkdir -p ~/.ssh; cat >> ~/.ssh/authorized_keys"
```

:::::::::::::::::::::::::
