---
title: Setup
---

## Instructors

1. **Arun Seetharam, Ph.D.**: Arun is a lead bioinformatics scientist at Purdue University’s Rosen Center for Advanced Computing. With extensive expertise in comparative genomics, genome assembly, annotation, single-cell genomics,  NGS data analysis, metagenomics, proteomics, and metabolomics. Arun supports a diverse range of bioinformatics projects across various organisms, including human model systems.

2. **Michael Carlson**, Ph.D.: Michael is a Senior Computational Scientist at Purdue University's Rosen Center for Advanced Computing (RCAC). Michael has a background in computational physics, specifically hypersonic materials. He also leads many introductory workshops in the High-Performance Computing domain.


## Schedule (10/06/2026)


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

## Prerequisites

This workshop assumes:

- **Basic Linux/command-line skills**: navigating directories, running commands, editing files
- **Basic R skills**: installing packages, reading/writing data, creating plots
- **A Purdue HPC account**: access to the Negishi cluster or Scholar (provided)
- **An SSH client**: terminal (macOS/Linux) or PuTTY/MobaXterm (Windows)
- **Genomics knowledge**: understanding of genes, transcripts, and genome structure

No prior experience with single-cell RNA-seq is required.

## SSH Setup

You need SSH access to the Negishi cluster for Episode 2 (raw data processing) and to copy the workshop data. Follow the instructions for your operating system below.

::::::::::::::::::::::::::::::::::::::: discussion

## Connecting to the Cluster

You will need SSH access to the Negishi cluster at Purdue. Choose the instructions for your operating system below.

:::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::: solution

### Windows

1. Download and install [MobaXterm](https://mobaxterm.mobatek.net/) (recommended) or [PuTTY](https://www.putty.org/)
2. Open MobaXterm and click **Session > SSH**
3. Set **Remote host** to `negishi.rcac.purdue.edu`
4. Check **Specify username** and enter your Purdue career account username
5. Click **OK** and enter your password when prompted
6. Complete Microsoft two-factor authentication

:::::::::::::::::::::::::

:::::::::::::::: solution

### macOS

1. Open **Terminal** (Applications > Utilities > Terminal)
2. Connect to the cluster:

```bash
ssh your_username@negishi.rcac.purdue.edu
```

1. Enter your password (characters may not appear, but your password is being entered)
2. Complete Microsoft two-factor authentication

:::::::::::::::::::::::::

:::::::::::::::: solution

### Linux

1. Open your terminal emulator
2. Connect to the cluster:

```bash
ssh your_username@negishi.rcac.purdue.edu
```

1. Enter your password (characters may not appear, but your password is being entered)
2. Complete Microsoft two-factor authentication

:::::::::::::::::::::::::


## Data Setup

To copy only the training data:

```bash
rsync -avP /scratch/negishi/aseethar/rnaseq-workshop ${RCAC_SCRATCH}/
```

A completed version of the workshop data is available at:

```
/scratch/negishi/aseethar/rnaseq-workshop_results
```

You can copy it to your scratch space using:

```bash
rsync -avP /scratch/negishi/aseethar/rnaseq-workshop_results ${RCAC_SCRATCH}/
```

Use this folder **only if you are unable to complete the exercises during the workshop**.



