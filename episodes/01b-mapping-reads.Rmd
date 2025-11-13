---
source: Rmd
title: Mapping reads to a reference genome
teaching: 15
exercises: 30
---

:::::::::::::::::::::::::::::::::::::: questions 

- How can we map raw reads to a reference genome?
- What programs are available for mapping reads to a reference genome?
- What opitons and parameters need to be considered when mapping?
- How do to access the mapping results?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the process of mapping raw reads to a reference genome.
- Fine tune the mapping process by adjusting parameters.
- Access the mapping results.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

In this episode, we will learn how to map raw reads to a reference genome. Mapping is the process of aligning short reads to a reference genome. The goal of mapping is to determine the location of the reads in the reference genome. We will use the RNA-seq data that we downloaded (soft-linked) in the previous episode and map them to the reference genome using the STAR aligner. STAR is a fast and accurate RNA-seq aligner that can align reads to the reference genome with high sensitivity and specificity.



:::::::::::::::::::::::::::::::::::::::  prereq

## What mapping programs should I use?

One major requirement for a suitable RNAseq mapping program is that it should be splice-aware. This means the program must be able to handle reads that span exon-exon junctions, which is crucial for RNAseq data since the reads are derived from spliced transcripts. `STAR` is a widely used splice-aware aligner for RNAseq, known for its speed and accuracy when mapping reads to a reference genome. Other popular RNAseq aligners include `HISAT2` and, for certain applications, `Bowtie2`. It’s important not to use DNA aligners like `BWA`, which are not designed to handle spliced reads.

:::::::::::::::::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::::::::  prereq

## How does read quality affect alignment performance?

The quality of the reads can have a significant impact on the alignment performance and thus your downstream analysis. 

- **Uneven read quality**: poor quality of the reads can result in partial alignments or truncated read mapping, leaving you with large portions of the read unmapped.
- **Overrepresented sequences**: such as ribosomal RNA, or microbial contamination, can lead to too many multi-mapped reads.
- **Poor library preparation**: can result in RNAseq reads that consist primarily of adapter sequences without any actual biological content - this can result in too many unmapped reads.
- **Over-processing** Short reads from over-trimming can shorten reads to the point where they are too short to align uniquely, leading to either multimapping or unmapped reads.

It is crucial to check the quality of the reads before mapping and also alignment statistics after mapping.

:::::::::::::::::::::::::::::::::::::::


## Preparing for read mapping

Before we start mapping the reads, we need to prepare the reference genome. Mapping programs require the reference genome to be indexed. Indexing is the process of creating a data structure that allows the mapping program to quickly search the reference genome for the best alignment of the reads. `STAR` requires the reference genome to be indexed before mapping. Mapping programs can also use `GTF` or `GFF` files to guide the alignment process. These files contain information about the gene models in the reference genome, such as the location of exons and introns, that can speed up the alignment process. 

We will use the reference genome and annotations that we downloaded in the previous episode and index it using `STAR`.

## Indexing

STAR has `--runMode` option where you can alter the function of the program. When `genomeGenerate` is specified with this option, it can be used for indexing. 

:::::::::::::::::::::::::::::::::::::::  challenge

## Generating genome indexes

There are many options for the `STAR` program to index the reference genome. You will need the reference genome fasta file that needs to indexed, but there are many other options as well. Can you [check the manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) and find options that are important for indexing the reference genome? 



::::::::::::::::::::::::::::::::::: solution

In page 5, Section 2, "Generating genome indexes" discuss the options for indexing the reference genome. You will need to use `--runMode genomeGenerate` for indexing. Important options for indexing are:

- `runThreadN` number of threads to use
- `genomeDir` output directory for the genome indexes
- `genomeFastaFiles` reference genome fasta file
- `sjdbGTFfile` the GTF file
- `sjdbOverhang` default value of 100 works well, for most cases this can be omitted

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

We will also need SLURM directives to submit the job to the cluster. The number of threads to use can be specified using the `--cpus-per-task` directive. The script to index the reference genome using `STAR` would look like this:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --job-name=index_genome
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.out
ml biocontainers
ml star
genome="/home/aseethar/rcac_rnaseq/data/GRCm39.primary_assembly.genome.fa"
genomeDir="/home/aseethar/rcac_rnaseq/data/star_index"
gtf="/home/aseethar/rcac_rnaseq/data/gencode.vM35.primary_assembly.basic.annotation.gtf"
STAR \
   --runThreadN 4\
   --runMode genomeGenerate \
   --genomeDir ${genomeDir} \
   --genomeFastaFiles ${genome} \
   --sjdbGTFfile ${gtf}
```

The script can be saved in a file, say `index_genome.sh`, within the `scripts` directory. It should be submitted as a slurm job after making it executable.

```bash
chmod +x index_genome.sh
sbatch index_genome.sh
```

It takes about 30-45 minutes to complete the indexing process.

:::::::::::::::::::::::::::::::::::::::  prereq

## Is GTF file necessary for indexing?

The GTF file is not necessary for indexing the reference genome, but will be useful for alignment process. The GTF file contains information about the gene models in the reference genome, such as the location of exons and introns, that can speed up the alignment process. If you have a GTF file, it is recommended to use it when indexing the reference genome.

:::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  prereq

## Can I use GFF3 file instead?

Yes, you can use GFF3 file instead of GTF file. STAR can handle both GTF and GFF3 files. In addition to the aforementioned options, for GFF3 formatted annotations you need to use `--sjdbGTFtagExonParentTranscript Parent`. But if you have issues with GFF3, you can convert it to GTF using `gffread` from the `Cufflinks` package.

:::::::::::::::::::::::::::::::::::::::




## Mapping reads to the reference genome

After indexing, our genome is ready for mapping reads. The `STAR` program has many options that can be used to fine-tune the mapping process. You can [check the manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) and find options that are important for mapping reads to the reference genome. 


For mapping, you will need to use `--runMode alignReads` important options for aligning reads are:

- `runThreadN` number of threads to use
- `readFilesIn` fastq files to be aligned, can be single (1 file) or paired-end (2 files)
- `genomeDir` reference genome index directory

For output options, you can use:

- `outFileNamePrefix` prefix for all output files
- `outSAMtype`: output filetype (SAM default), but can be BAM, CRAM, etc.
- `outSAMunmapped`: handling of unmapped reads

:::::::::::::::::::::::::::::::::::::::  prereq

## What other mapping options should I consider?

STAR's default parameters are optimized for mammalian genomes. If you are using for other genomes like plant, fungi or any other species, you may need to fine-tune your mapping options:

- When using compressed reads, option `readFilesCommand` with _UncompressionCommand_ can be used. If using gzipped files (`*.gz`) use `--readFilesCommand zcat` OR `--readFilesCommand gunzip -c`. For `bzip2` compressed files, use `--readFilesCommand bunzip2 -c`.
- For controling the multi-mapped reads, options like `outFilterMultimapNmax` (default 10) and `winAnchorMultimapNmax` can be used. Useful when you have polyploid genomes or repetitive regions.
- For controlling mapping sensitivity, options like `outFilterMismatchNmax` (default 10) and `outFilterMismatchNoverLmax` can be used.
- Shorter reads may require lower `outFilterMismatchNmax` and `outFilterMismatchNoverLmax` values.

:::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Script for mapping reads

Write a script to map the reads to the reference genome using `STAR`. You will need the reference genome index directory, fastq files to be aligned, and output directory. Can you write a mapping command for `STAR`?


::::::::::::::::::::::::::::::::::: solution

```bash
STAR --runMode alignReads \
     --runThreadN 12 \
     --genomeDir /home/aseethar/rcac_rnaseq/data/star_index \
     --readFilesIn SRR1234567_1.fastq.gz SRR1234567_2.fastq \
     --outFileNamePrefix SRR1234567 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within
```
:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

Since we have 45 sets of fastq files, we can write an array job for mapping all the reads. The script to map the reads using `STAR` would look like this:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=8:00:00
#SBATCH --job-name=read_mapping
#SBATCH --array=1-45    
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.out
ml biocontainers
ml star
# Define paths
FASTQ_DIR="/home/aseethar/rcac_rnaseq/data/fastq"
GENOME_INDEX="/home/aseethar/rcac_rnaseq/data/star_index"
OUTPUT_DIR="/home/aseethar/rcac_rnaseq/data/mapping"
# Create an array of FASTQ file pairs
FASTQ_R1=($FASTQ_DIR/*_R1_*.fastq.gz)   # List of all R1 files
FASTQ_R2=($FASTQ_DIR/*_R2_*.fastq.gz)   # List of all R2 files
# Get the FASTQ file for this array job
R1=${FASTQ_R1[$SLURM_ARRAY_TASK_ID-1]}
R2=${FASTQ_R2[$SLURM_ARRAY_TASK_ID-1]}
# Extract sample name from FASTQ file
SAMPLE=$(basename $R1 _R1_001.fastq.gz)
# Run STAR alignment for the current pair of FASTQ files
STAR --runThreadN 4 \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $R1 $R2 \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/$SAMPLE \
     --outSAMunmapped Within \
     --outSAMtype BAM SortedByCoordinate
```

The job can be submitted as an array job using the following command:

```bash
chmod +x map_reads.sh
sbatch --array=1-45 map_reads.sh
```

Since we have already finished mapping the reads, we can skip this and move on to the next step.

:::::::::::::::::::::::::::::::::::::::  prereq

## Can I use STAR for counting reads?

Yes, STAR can be used for counting reads. STAR can output the number of reads that map to each gene in the reference genome. This can be useful for differential gene expression analysis. You can use the `--quantMode` option to enable read counting. The `--quantMode` option can be set to `GeneCounts` to output the number of reads that map to each gene.

:::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::: keypoints 


- Mapping Raw Reads: The process of mapping involves aligning short reads from sequencing data to a reference genome to determine their locations. STAR, a splice-aware RNA-seq aligner, is commonly used due to its speed and accuracy in handling exon-exon junctions.
   
- Impact of Read Quality: Poor read quality, such as over-represented sequences or low-quality bases, can lead to unmapped or multi-mapped reads. It is important to check the quality of reads before mapping and review alignment statistics after mapping to ensure data quality.

- Reference Genome Indexing: Before mapping, the reference genome needs to be indexed using tools like STAR. Indexing creates a data structure to facilitate rapid alignment, and additional files like GTF/GFF annotations can improve alignment accuracy.

- Optimizing Mapping: Fine-tuning mapping parameters (e.g., handling multi-mapped reads, controlling sensitivity) based on the dataset is essential. The STAR aligner offers numerous options that allow adjustment for specific genome types and read qualities.


::::::::::::::::::::::::::::::::::::::::::::::::

