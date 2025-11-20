---
source: Rmd
title: "A. Genome-based quantification (STAR + featureCounts)"
teaching: 30
exercises: 40
---

:::::::::::::::::::::::::::::::::::::: questions

- How do we map RNA-seq reads to a reference genome?
- How do we determine library strandness before alignment?
- How do we quantify reads per gene using featureCounts?
- How do we submit mapping jobs on RCAC clusters with SLURM?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Check RNA-seq library strandness using Salmon.
- Index a reference genome for STAR.
- Map reads to the genome using STAR.
- Quantify gene level counts with featureCounts.
- Interpret basic alignment and counting statistics.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

In this episode we move from raw FASTQ files to genome based gene counts.  
This workflow has three major steps:

1. Determine library strandness: STAR does not infer strandness automatically. `featureCounts` requires the correct setting.
2. Map reads to the genome using STAR: STAR is a splice aware aligner designed for RNA-seq data.
3. Count mapped reads per gene using `featureCounts`

This represents one of the two standard quantification strategies in RNA-seq.  
The alternative method (B. transcript based quantification with Salmon or Kallisto) will be covered in Episode 4B.

## Step 1: Checking RNA-seq library strandness

::::::::::::::::::::::::::::::::::::::: callout

## Why strandness matters

RNA-seq libraries can be:

- unstranded  
- stranded forward  
- stranded reverse  

Stranded protocols are common in modern kits, but old GEO and SRA datasets are often unstranded.  
featureCounts must be given the correct strand setting. If the wrong setting is used, most reads will be counted against incorrect genes.

STAR does not infer strandness automatically, so we must determine it **before** alignment.

:::::::::::::::::::::::::::::::::::::::

### Detecting strandness with Salmon

Salmon can infer strandness directly from raw FASTQ files without any alignment.  
We run Salmon with `-l A`, which tells it to explore all library types.

```bash
sinteractive -A rcac -p cpu -N 1 -n 4 --time=1:00:00

cd $SCRATCH/rnaseq-workshop
mkdir -p results/strand_check

module load biocontainers
module load salmon

salmon index --transcripts data/gencode.vM38.transcripts-clean.fa \
             --index data/salmon_index \
             --threads 4

salmon quant --index data/salmon_index \
             --libType A \
             --mates1 data/SRR33253286_1.fastq.gz \
             --mates2 data/SRR33253286_2.fastq.gz \
             --output results/strand_check \
             --threads 4
```
The important output is in `results/strand_check/lib_format_counts.json`:

::::::::::::::::::::::::::::::::::::: spoiler

## What does this file look like?

```json
{
    "read_files": "[ data/SRR33253286_1.fastq.gz, data/SRR33253286_2.fastq.gz]",
    "expected_format": "IU",
    "compatible_fragment_ratio": 1.0,
    "num_compatible_fragments": 21088696,
    "num_assigned_fragments": 21088696,
    "num_frags_with_concordant_consistent_mappings": 19116702,
    "num_frags_with_inconsistent_or_orphan_mappings": 2510542,
    "strand_mapping_bias": 0.5117980078362889,
    "MSF": 0,
    "OSF": 0,
    "ISF": 9783890,
    "MSR": 0,
    "OSR": 0,
    "ISR": 9332812,
    "SF": 1234245,
    "SR": 1276297,
    "MU": 0,
    "OU": 0,
    "IU": 0,
    "U": 0
}

```


> ### Interpreting `lib_format_counts.json`
>
> This JSON file reports Salmon's automatic library type detection.
>
> - **`"expected_format": "IU"`**: This is Salmon's conclusion. "I" means **Inward** (correct for paired-end reads) and "U" means **Unstranded**.
> - **`"strand_mapping_bias": 0.512`**: This is the key evidence. A value near 0.5 (50%) indicates that reads mapped equally to both the sense and antisense strands, the definitive sign of an **unstranded** library.
> - **`"ISF"`** and **`"ISR"`**: These are the counts for Inward-Sense-Forward (9.78M) and Inward-Sense-Reverse (9.33M) fragments. Because these values are almost equal, they confirm the ~50/50 split seen in the strand bias.

Common results:

| Code | Meaning                       |
| ---- | ----------------------------- |
| IU   | unstranded                    |
| ISR  | stranded reverse (paired end) |
| ISF  | stranded forward (paired end) |
| SR   | stranded reverse (single end) |
| SF   | stranded forward (single end) |

**Conclusion:** The data is **unstranded**. We will record this and supply it to STAR/featureCounts.

:::::::::::::::::::::::::::::::::::::

## Step 2: Preparing the reference genome for alignment

STAR requires a genome index. Indexing builds a searchable data structure that speeds up alignment.

Required input:

- genome FASTA file
- annotation file (GTF or GFF3)
- output directory for the index

Optionally, annotation improves splice junction detection and yields better alignments.


::::::::::::::::::::::::::::::::::::::: callout

## Why STAR?

- Splice-aware  
- Extremely fast  
- Produces high-quality alignments  
- Well supported and widely used  
- Compatible with downstream tools such as featureCounts, StringTie, and RSEM

Optional alternatives include HISAT2 and Bowtie2, but STAR remains the standard for high-depth Illumina RNA-seq.

:::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Read the STAR manual

Check the **"Generating genome indexes"** section of the  
[STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

Which options are required for indexing?

::::::::::::::::::::::::::::::::::: solution

Key options:

- `--runMode genomeGenerate`  
- `--genomeDir` directory to store the index  
- `--genomeFastaFiles` the FASTA file  
- `--sjdbGTFfile` annotation (recommended)  
- `--runThreadN` number of threads  
- `--sjdbOverhang` read length minus one (100 is a safe default)

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::


### Indexing the Genome with STAR

Below is a minimal SLURM job script for building the STAR genome index. Create a file named `index_genome.sh` in `$SCRATCH/rcac_rnaseq/scripts` with this content:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=workshop
#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=star_index
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load star

data_dir="${SCRATCH}/rnaseq-workshop/data"
genome="${data_dir}/GRCm39.primary_assembly.genome.fa"
gtf="${data_dir}/gencode.vM38.primary_assembly.basic.annotation.gtf"
indexdir="${data_dir}/star_index"

STAR \
    --runMode genomeGenerate \
    --runThreadN ${SLURM_CPUS_ON_NODE} \
    --genomeDir ${indexdir} \
    --genomeFastaFiles ${genome} \
    --sjdbGTFfile ${gtf}

```

Submit the job:

```bash
sbatch index_genome.sh
```

Indexing requires about 30 to 45 minutes.

::::::::::::::::::::::::::::::::::::::: spoiler

## What will the output look like?

The index directory is now ready for use in mapping and should have contents like this:
Location: `$SCRATCH/rcac_rnaseq/data`
Contents of `star_index/`:

```bash
star_index/
├── chrLength.txt
├── chrNameLength.txt
├── chrName.txt
├── chrStart.txt
├── exonGeTrInfo.tab
├── exonInfo.tab
├── geneInfo.tab
├── Genome
├── genomeParameters.txt
├── Log.out
├── SA
├── SAindex
├── sjdbInfo.txt
├── sjdbList.fromGTF.out.tab
├── sjdbList.out.tab
└── transcriptInfo.tab
```
:::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::::: callout

## Do I need the GTF for indexing?

No, but it is *strongly recommended*.  
Including it improves splice junction detection and increases mapping accuracy.

:::::::::::::::::::::::::::::::::::::::


## Step 3: Mapping reads with STAR

Once the genome index is ready, we can align reads.

For aligning reads, important STAR options are:

- `--runThreadN`
- `--genomeDir`
- `--readFilesIn`
- `--readFilesCommand zcat` for gzipped FASTQ files
- `--outSAMtype BAM SortedByCoordinate`
- `--outSAMunmapped Within`

STAR's defaults are tuned for mammalian genomes. If working with plants, fungi, or repetitive genomes, additional tuning may be required.


::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Write a STAR alignment command

Write a basic STAR alignment command using your genome index and paired FASTQ files.

::::::::::::::::::::::::::::::::::: solution


```bash
STAR \
    --runThreadN 12 \
    --genomeDir data/star_index \
    --readFilesIn data/SRR1234567_1.fastq.gz data/SRR1234567_2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix results/mapping/SRR1234567 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within
```

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::


### Mapping many FASTQ files with a SLURM array

When having multiple samples, array jobs simplify submission.

First, we will create `samples.txt` file with just the sample names

```bash
cd $SCRATCH/rnaseq-workshop/data
ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > ${SCRATCH}/rnaseq-workshop/scripts/samples.txt
```

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=workshop
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --job-name=read_mapping
#SBATCH --array=1-6
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load star

FASTQ_DIR="$SCRATCH/rnaseq-workshop/data"
GENOME_INDEX="$SCRATCH/rnaseq-workshop/data/star_index"
OUTDIR="$SCRATCH/rnaseq-workshop/results/mapping"

mkdir -p $OUTDIR

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

R1=${FASTQ_DIR}/${SAMPLE}_1.fastq.gz
R2=${FASTQ_DIR}/${SAMPLE}_2.fastq.gz

STAR \
    --runThreadN ${SLURM_CPUS_ON_NODE} \
    --genomeDir $GENOME_INDEX \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $OUTDIR/$SAMPLE \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate
```

Submit:

```bash
sbatch map_reads.sh
```


::::::::::::::::::::::::::::::::::::::: discussion

## Why use an array job?

Think about this:

- We have *six* samples, each with two FASTQ files  
- Having six separate slurm jobs means:  
  - Copy-pasting the same command six times  
  - Higher chance of typos or mistakes  
- if mapping with single slurm file, sequential execution = slower
- Samples can run independently  

Why is a job array a good choice?


Arrays allow us to run identical commands across many inputs:

- All samples are processed the same way  
- Runs in parallel  
- Simplifies submission  
- Reduces copy-paste errors  
- Minimizes runtime cost on the cluster


:::::::::::::::::::::::::::::::::::::::


## Step 4: Assessing alignment quality

STAR produces several useful files, particularly:

- `Log.final.out`  
  Contains mapping statistics, including % uniquely mapped.

- `Aligned.sortedByCoord.out.bam`  
  Ready for downstream quantification.

To summarize alignment statistics across samples, we use MultiQC:

```bash
cd $SCRATCH/rnaseq-workshop
module load biocontainers
module load multiqc

multiqc results/mapping -o results/qc_alignment
```

### Key metrics to review

- **Uniquely mapped reads** — typically 60–90% is good  
- **Multi-mapped reads** — depends on genome, but >20% is concerning  
- **Unmapped reads** — check reasons (mismatches? too short?)  
- **Mismatch rate** — high values indicate quality or contamination issues  

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Examine your alignment metrics

Using your MultiQC alignment report:

1. Which sample has the highest uniquely mapped percentage?  
2. Are any samples clear outliers?  
3. Does mapped/unmapped reads appear consistent across samples?

::::::::::::::::::::::::::::::::::: solution

Example interpretation for this dataset:

- All samples have ~80% uniquely mapped reads.
- No sample is an outlier.
- mapped/unmapped reads consistent across samples

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::



## Step 5: Quantifying gene counts with featureCounts

Once the reads have been aligned and sorted, we can quantify how many read pairs map to each gene. This gives us gene-level counts that we will later use for differential expression. For alignment-based RNA-seq workflows, `featureCounts` is the recommended tool because it is fast, robust, and compatible with most gene annotation formats.

`featureCounts` uses two inputs:

1. A set of aligned BAM files (sorted by coordinate)
2. A GTF annotation file describing gene models

It assigns each read (or read pair) to a gene based on the exon boundaries.

::::::::::::::::::::::::::::::::::::::: callout

## Important: stranded or unstranded?

`featureCounts` needs to know whether your dataset is stranded.

- `-s 0` unstranded
- `-s 1` stranded forward
- `-s 2` stranded reverse

From our Salmon strandness check, this dataset behaves as **unstranded**, so we will use `-s 0`.

:::::::::::::::::::::::::::::::::::::::

### Running featureCounts

We will create a new directory for count output:

```bash
cd $SCRATCH/rnaseq-workshop/scripts
mkdir -p $SCRATCH/rnaseq-workshop/results/counts
```

Below is an example SLURM script to count all BAM files at once.

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=workshop
#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=featurecounts
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load subread

GTF="$SCRATCH/rnaseq-workshop/data/gencode.vM38.annotation.gtf"
BAMS="$SCRATCH/rnaseq-workshop/results/mapping/*.bam"
output_dir="$SCRATCH/rnaseq-workshop/results/counts"

featureCounts \
  -T ${SLURM_CPUS_PER_TASK} \
  -p \
  -s 0 \
  -t exon \
  -g gene_id \
  -a ${GTF} \
  -o ${output_dir}/gene_counts.txt \
  ${BAMS}
```

Submit the job:

```bash
sbatch count_features.sh
```

`gene_counts.txt` will contain:

- one row per gene
- one column per sample

::::::::::::::::::::::::::::::::::::::: callout
## Generate a clean count matrix

```bash
cd $SCRATCH/rnaseq-workshop/results/counts
cut -f 1,7- gene_counts.txt |\
  sed 's+/scratch/negishi/'$USER'/rnaseq-workshop/results/mapping/++g' |\
  grep -v "^#" |\
  sed 's/Aligned.sortedByCoord.out.bam//g' \
  > gene_counts_clean.txt
```
:::::::::::::::::::::::::::::::::::::::

## Step 6: QC on counts

Before proceeding to differential expression, it is good practice to perform some QC on the count data. We can again use MultiQC to summarize featureCounts statistics:

```bash
cd $SCRATCH/rnaseq-workshop
mkdir -p results/qc_counts
module load biocontainers
module load multiqc
multiqc results/counts -o results/qc_counts
```

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: Examine your counting metrics

Using your MultiQC featureCounts report:

1. Which sample has the lowest percent assigned, and what are common reasons for lower assignment?
2. Which unassigned category is the largest, and what does it indicate?
3. Does the sample with the most assigned reads also have the highest percent assigned?
4. What does "Unassigned: No Features" mean, and why might it occur?
5. Are multi-mapped reads a concern for RNA-seq differential expression?

::::::::::::::::::::::::::::::::::: solution

Example interpretation for this dataset:

1. SRR33253289 has the lowest assignment rate. Causes include low complexity, incomplete annotation, or more intronic reads.
2. Unmapped or No Features dominate. These reflect reads that do not align or do not overlap annotated exons.
3. No. Total assigned reads depend on depth, while percent assigned reflects library quality.
4. Reads mapped outside annotated exons, often from intronic regions, unannotated transcripts, or incomplete annotation.
5. Multi-mapping is expected for paralogs and repeats. It is usually fine as long as the rate is consistent across samples.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::



## Summary

::::::::::::::::::::::::::::::::::::: keypoints

- STAR is a fast splice-aware aligner widely used for RNA-seq.
- Genome indexing requires FASTA and optionally GTF for improved splice detection.
- Mapping is efficiently performed using SLURM array jobs.
- Alignment statistics (from `Log.final.out` + MultiQC) must be reviewed.
- Strandness should be inferred using aligned BAM files.
- Final BAM files are ready for downstream counting.
- `featureCounts` is used to obtain gene-level counts for differential expression.
- Exons are counted and summed per gene.
- The output is a simple matrix for downstream statistical analysis.
::::::::::::::::::::::::::::::::::::::::::::::::
