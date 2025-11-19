---
source: Rmd
title: Downloading and organizing files
teaching: 20
exercises: 15
---

:::::::::::::::::::::::::::::::::::::: questions

- What files are required to process raw RNA seq reads?
- Where do we obtain raw reads, reference genomes, annotations, and transcript sequences?
- How should project directories be organized to support a smooth workflow?
- What practical considerations matter when downloading and preparing RNA seq data?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Identify the essential inputs for RNA seq analysis: FASTQ files, reference genome, annotation, and transcript sequences.
- Learn where and how to download RNA seq datasets and reference materials.
- Create a clean project directory structure suitable for quality control, alignment, and quantification.
- Understand practical issues such as compressed FASTQ files and adapter trimming.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

Before we can perform quality control, mapping, or quantification, we must first gather the files required for RNA seq analysis and organize them in a reproducible way. This episode introduces the essential inputs of an RNA seq workflow and demonstrates how to obtain and structure them on a computing system.

We begin by discussing raw FASTQ reads, reference genomes, annotation files, and transcript sequences. We then download the dataset used throughout the workshop and set up the directory structure that supports downstream steps.

Most FASTQ files arrive compressed as `.fastq.gz` files. Modern RNA seq tools accept compressed files directly, so uncompressing them is usually unnecessary.

::::::::::::::::::::::::::::::::::::::: prereq

## Should I uncompress FASTQ files?

In almost all RNA seq workflows, you should keep FASTQ files compressed. Tools such as FastQC, HISAT2, STAR, salmon, and fastp can read `.gz` files directly. Keeping files compressed saves space and reduces input and output overhead.

:::::::::::::::::::::::::::::::::::::::

Many library preparation protocols introduce adapter sequences at the ends of reads. For most alignment based analyses, explicit trimming is not required because aligners soft clip adapter sequences automatically.

::::::::::::::::::::::::::::::::::::::: prereq

## Should I trim adapters for RNA seq analysis?

In most cases, no. Aligners such as HISAT2 and STAR soft clip adapter sequences. Overly aggressive trimming can reduce mapping rates or distort read length distributions. Trimming is needed only for specific applications such as transcript assembly where uniform read lengths are important.  
More information: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4766705/>

:::::::::::::::::::::::::::::::::::::::

## Downloading the data and project organization

For this workshop, we use a publicly available dataset that investigates the effects of cardiomyocyte specific Nat10 knockout on heart development in mice (BioProject PRJNA1254208, GEO GSE295332). The study includes three biological replicates per group, sequenced on an Illumina NovaSeq X Plus in paired end mode.

Before downloading files, we create a reproducible directory layout for raw data, scripts, mapping outputs, and count files.


::::::::::::::::::::::::::::::::::::::: challenge

## Create the directories needed for this episode

Create a working directory for the workshop (e.g., `rnaseq-workshop`). Inside it, create four subdirectories:

- `data` : raw FASTQ files, reference genome, annotation  
- `scripts` : custom scripts, SLURM job files  
- `results` : alignment, counts and various other outputs  


::::::::::::::::::::::::::::::::::: solution


```bash
cd $RCAC_SCRATCH
mkdir -p rnaseq-workshop/{data,scripts,results}
```

The resulting directory structure should look like this:

```bash
aseethar@scholar-fe03:[aseethar] $ pwd
/scratch/scholar/aseethar
aseethar@scholar-fe03:[aseethar] $ tree rnaseq-workshop/
rnaseq-workshop/
├── results
├── data
└── scripts
```
:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## What files do I need for RNA-seq analysis?

::::::::::::::::::::::::::::::::::::::: prereq

You will need:


1. Raw reads (FASTQ)
   - Usually `.fastq.gz` files containing nucleotide sequences and per-base quality scores.

**_For alignment based workflows_**

2. Reference genome (FASTA)
   - Contains chromosome or contig sequences.

3. Annotation file (GTF/GFF)
   - Describes gene models: exons, transcripts, coordinates.

**_For transcript-level quantification workflows_**

4. Transcript sequences (FASTA)  
   - Required for transcript based quantification with tools such as salmon, kallisto, and RSEM.

Reference files must be consistent with each other (same genome version).
To ensure a smooth analysis workflow, keep these files logically organized. Good directory structure minimizes mistakes, simplifies scripting, and supports reproducibility.

:::::::::::::::::::::::::::::::::::::::


## Downloading the data

We now download:

- the GRCm39 reference genome  
- the corresponding GENCODE annotation
- transcript sequences (if using transcript based quantification)
- raw FASTQ files from SRA

### Reference genome and annotation

We will use **GENCODE GRCm39 (M38)** as the reference.

::::::::::::::::::::::::::::::::::::::: challenge

## Annotation

Visit the GENCODE mouse page: <https://www.gencodegenes.org/mouse/>.  
Which annotation file (GTF/GFF3 files section) should you download, and why?

::::::::::::::::::::::::::::::::::: solution

Use the **basic gene annotation in GTF format** (e.g., `gencode.vM38.primary_assembly.basic.annotation.gtf.gz`).  
It contains essential gene and transcript models without pseudogenes or other biotypes not typically used for RNA-seq quantification.

<div class="figure" style="text-align: center">
<img src="fig/01_fq2counts/gencode-gtf.png" alt="Basic gene annotation in GTF format"  />
<p class="caption">Basic gene annotation in GTF format</p>
</div>

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: challenge

## Reference genome

Select the corresponding **Genome sequence (GRCm39) FASTA** (e.g., `GRCm39.primary_assembly.genome.fa.gz`).

::::::::::::::::::::::::::::::::::: solution

<div class="figure" style="text-align: center">
<img src="fig/01_fq2counts/gencode-genome.png" alt="Reference genome file in FASTA format"  />
<p class="caption">Reference genome file in FASTA format</p>
</div>

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

### Transcript sequences

Transcript level quantification requires a transcript FASTA file consistent with the annotation.
GENCODE provides transcript sequences for all transcripts in GRCm39 (`gencode.vM38.transcripts.fa.gz`).

::::::::::::::::::::::::::::::::::::::: challenge

## Transcript sequence file

Where on the GENCODE page can you find the transcript FASTA file?

::::::::::::::::::::::::::::::::::: solution

The transcript FASTA file is available in the same release directory as the genome and annotation files (under FASTA files). It is named `gencode.vM38.transcripts.fa.gz`.

<div class="figure" style="text-align: center">
<img src="fig/01_fq2counts/gencode-transcripts.png" alt="Transcript sequences file in FASTA format"  />
<p class="caption">Transcript sequences file in FASTA format</p>
</div>

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::


### Downloading all reference files

```bash
cd ${SCRATCH}/rnaseq-workshop

GTFlink="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.basic.annotation.gtf.gz"
GENOMElink="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz"
CDSLink="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz"

# Download to 'data'
wget -P data ${GENOMElink}
wget -P data ${GTFlink}
wget -P data ${CDSLink}

# Uncompress
cd data
gunzip GRCm39.primary_assembly.genome.fa.gz
gunzip gencode.vM38.primary_assembly.basic.annotation.gtf.gz
gunzip gencode.vM38.transcripts.fa.gz
```
The transcript file header has multiple pieces of information which often interfere with downstream analysis. We create a simplified version containing only the transcript IDs.

```bash
grep ">" gencode.vM38.transcripts.fa |head -1
cut -f 1 -d "|" gencode.vM38.transcripts.fa > gencode.vM38.transcripts-clean.fa
grep ">" gencode.vM38.transcripts-clean.fa |head -1
```


## Raw reads

The experiment we use investigates how **NAT10** influences **mouse heart development**.

> **Study:** The effects of NAT10 on heart development in mice  
> **Organism:** Mus musculus  
> **Design:** Nat10^flox/flox vs Nat10-CKO (n=3 per group), P10 hearts  
> **Sequencing:** Illumina NovaSeq X Plus, paired-end  
> **BioProject:** PRJNA1254208  
> **GEO:** GSE295332  

<div class="figure" style="text-align: center">
<img src="fig/01_fq2counts/geo-db.png" alt="FASTQ files from GEO"  />
<p class="caption">FASTQ files from GEO</p>
</div>

::::::::::::::::::::::::::::::::::::::: challenge

## Where do you find FASTQ files?

GEO provides metadata and experimental details, but FASTQ files are stored in the **SRA**.
How do you obtain the complete list of sequencing runs?

::::::::::::::::::::::::::::::::::: solution

Navigate to the SRA page for the BioProject and use **Run Selector** to list runs.  
On the Run Selector page, click **Accession List** to download all SRR accession IDs.

<div class="figure" style="text-align: center">
<img src="fig/01_fq2counts/srr-runs.png" alt="FASTQ files from SRA Run Selector"  />
<p class="caption">FASTQ files from SRA Run Selector</p>
</div>

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

### Download FASTQ files with SRA Toolkit

```bash
sinteractive -A rcac -p cpu -N 1 -n 4 --time=1:00:00
cd ${SCRATCH}/rnaseq-workshop/data
module load biocontainers
module load sra-tools
while read SRR; do 
   fasterq-dump --threads 4 --progress $SRR;
   gzip ${SRR}_1.fastq
   gzip ${SRR}_2.fastq
done<SRR_Acc_List.txt
```

Since downloading may take a long time, we provide a pre-downloaded copy of all FASTQ files. The directory, after downloading, should look like this:

```bash
data
├── annot.tsv
├── gencode.vM38.primary_assembly.basic.annotation.gtf
├── gencode.vM38.transcripts.fa
├── GRCm39.primary_assembly.genome.fa
├── mart.tsv
├── SRR33253285_1.fastq.gz
├── SRR33253285_2.fastq.gz
├── SRR33253286_1.fastq.gz
├── SRR33253286_2.fastq.gz
├── SRR33253287_1.fastq.gz
├── SRR33253287_2.fastq.gz
├── SRR33253288_1.fastq.gz
├── SRR33253288_2.fastq.gz
├── SRR33253289_1.fastq.gz
├── SRR33253289_2.fastq.gz
├── SRR33253290_1.fastq.gz
├── SRR33253290_2.fastq.gz
└── SRR_Acc_List.txt
```

::::::::::::::::::::::::::::::::::::: keypoints

- RNA-seq requires three main inputs: FASTQ reads, a reference genome (FASTA), and gene annotation (GTF/GFF).
- _In lieu_ of a reference genome and annotation, transcript sequences (FASTA) can be used for transcript-level quantification.
- Keeping files compressed and well-organized supports reproducible analysis.
- Reference files must match in genome build and version.
- Public data repositories such as GEO and SRA provide raw reads and metadata.
- Tools like `wget` and `fasterq-dump` enable programmatic, reproducible data retrieval.

::::::::::::::::::::::::::::::::::::::::::::::::
