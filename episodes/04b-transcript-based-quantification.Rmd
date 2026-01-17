---
source: Rmd
title: "B. Transcript-based quantification (Kallisto)"
teaching: 30
exercises: 40
---

:::::::::::::::::::::::::::::::::::::: questions

- How do we quantify expression without genome alignment?
- What inputs does Kallisto require?
- How do we run Kallisto for paired-end RNA-seq data?
- How do we interpret transcript-level outputs (TPM, est_counts)?
- How do we summarize transcripts to gene-level counts using tximport?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Build a transcriptome index for Kallisto.
- Quantify transcript abundance directly from FASTQ files.
- Understand key Kallisto output files.
- Summarize transcript-level estimates to gene-level counts.
- Prepare counts for downstream differential expression.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

In the previous episode, we mapped reads to the **genome** using STAR and generated **gene-level counts** with featureCounts.

In this episode, we use an alternative quantification strategy:

1. Build a *transcriptome* index
2. Quantify expression directly from FASTQ files using **Kallisto** (*without alignment*)
3. Summarize transcript estimates back to gene-level counts

This workflow is faster, uses less storage, and models transcript-level uncertainty.

::::::::::::::::::::::::::::::::::::::: callout

## Why Kallisto?

Kallisto uses pseudo-alignment to rapidly quantify transcript abundances without traditional read mapping. Key advantages include:

- **Speed**: Kallisto is extremely fast, typically completing in 2-3 minutes per sample.
- **Accuracy**: Comparable accuracy to alignment-based methods for quantification.
- **Bootstrap support**: Native support for uncertainty estimation via bootstraps (useful for sleuth).
- **Low resource usage**: No large BAM files generated, minimal storage requirements.

Kallisto is ideal when you do not need alignment files, want fast quantification, or plan to use sleuth for differential transcript analysis.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: callout

## When should I use transcript-based quantification?

Kallisto is ideal when:

- You do not need alignment files.
- You want transcript-level TPMs.
- You want fast quantification.
- Storage is limited (no BAM files generated).
- You plan to use sleuth for differential transcript expression.

It is *not* ideal if you need splice junctions, variant calling, or visualization in IGV (all those depend on alignment files).

:::::::::::::::::::::::::::::::::::::::


## Step 1: Preparing the transcriptome reference

Kallisto requires a **transcriptome FASTA file** (all annotated transcripts). We previously downloaded the transcripts file from GENCODE for this purpose (`gencode.vM38.transcripts-clean.fa`).

::::::::::::::::::::::::::::::::::::::: callout

## Why transcriptome choice matters

Transcript-level quantification **inherits all assumptions of the annotation**.
Missing or incorrect transcripts → biased TPM estimates.

:::::::::::::::::::::::::::::::::::::::

### Building the Kallisto index

The Kallisto index is lightweight and quick to build:

```bash
cd $SCRATCH/rnaseq-workshop
mkdir -p data/kallisto_index

module load biocontainers
module load kallisto

kallisto index \
    -i data/kallisto_index/transcripts.idx \
    data/gencode.vM38.transcripts-clean.fa
```

This creates an index file that Kallisto uses for pseudo-alignment. Index building typically takes 2-3 minutes for a mammalian transcriptome.

::::::::::::::::::::::::::::::::::::::: callout

## Kallisto vs Salmon index

Note that Kallisto and Salmon indices are **not interchangeable**. Each tool has its own index format:

- Kallisto: Single `.idx` file
- Salmon: Directory with multiple files

You must build a separate index for each tool.

:::::::::::::::::::::::::::::::::::::::

## Step 2: Quantifying transcript abundances

The core Kallisto command is `kallisto quant`. We use the transcript index and paired-end FASTQ files.

```bash
mkdir -p $SCRATCH/rnaseq-workshop/results/kallisto_quant
```

Example command:

```bash
kallisto quant \
    -i data/kallisto_index/transcripts.idx \
    -o results/kallisto_quant/WT_Bcell_mock_rep1 \
    -b 100 \
    -t 16 \
    data/WT_Bcell_mock_rep1_R1.fastq.gz \
    data/WT_Bcell_mock_rep1_R2.fastq.gz
```

Important flags:

* `-i` → path to the Kallisto index
* `-o` → output directory for this sample
* `-b 100` → number of bootstrap samples (useful for uncertainty estimation and sleuth)
* `-t 16` → number of threads to use

::::::::::::::::::::::::::::::::::::::: callout

## Strand-specific libraries

For **strand-specific** libraries, add the appropriate flag:

- `--rf-stranded` for reverse-stranded libraries (e.g., Illumina TruSeq stranded)
- `--fr-stranded` for forward-stranded libraries

Our dataset is **unstranded** (detected as `IU` in Step 0), so we omit these flags.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: callout

## Bootstrap samples

The `-b 100` flag generates 100 bootstrap samples. These are useful for:

- Estimating uncertainty in transcript abundance
- Required if you plan to use **sleuth** for differential expression
- Can be omitted if you only use DESeq2/edgeR (for faster runtime)

For DESeq2 analysis only, you can use `-b 0` or omit the flag entirely.

:::::::::::::::::::::::::::::::::::::::

## Step 3: Running Kallisto for many samples

We use a SLURM array job for efficiency and consistency.

First, create a sample list (or reuse the one from Episode 4A):

```bash
cd $SCRATCH/rnaseq-workshop/data
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > $SCRATCH/rnaseq-workshop/scripts/samples.txt
```

Array job:

```bash
!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=workshop
#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=kallisto_quant
#SBATCH --array=1-8
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load kallisto

DATA="$SCRATCH/rnaseq-workshop/data"
INDEX="$SCRATCH/rnaseq-workshop/data/kallisto_index/transcripts.idx"
OUT="$SCRATCH/rnaseq-workshop/results/kallisto_quant"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

mkdir -p ${OUT}/${SAMPLE}

R1=${DATA}/${SAMPLE}_R1.fastq.gz
R2=${DATA}/${SAMPLE}_R2.fastq.gz

kallisto quant \
    -i $INDEX \
    -o $OUT/$SAMPLE \
    -b 100 \
    -t ${SLURM_CPUS_ON_NODE} \
    $R1 $R2 &> $OUT/$SAMPLE/${SAMPLE}.log
```

Submit:

```bash
sbatch quant_kallisto.sh
```

The typical run time is ~2-3 minutes per sample (faster than Salmon).

::::::::::::::::::::::::::::::::::::::: discussion

## Why use an array job here?

Kallisto runs very quickly and produces small output. The advantage of array jobs is *consistency*: all samples are quantified with identical parameters and metadata. This ensures that transcript-level TPM and count estimates are directly comparable across samples.

:::::::::::::::::::::::::::::::::::::::

### Inspect results

After Kallisto finishes, each sample directory contains a small set of output files:

```
WT_Bcell_mock_rep1/
├── abundance.tsv          # IMPORTANT: transcript-level TPM, counts, effective lengths
├── abundance.h5           # IMPORTANT: HDF5 format with bootstrap data (for sleuth)
└── run_info.json          # IMPORTANT: run parameters and mapping statistics
```

### What is inside *abundance.tsv*?

This file contains the transcript-level quantification results:

* **target_id**: transcript ID
* **length**: transcript length
* **eff_length**: effective length adjusted for fragment distribution
* **est_counts**: estimated number of reads assigned to that transcript
* **tpm**: within-sample normalized abundance (Transcripts Per Million)

Example:

```text
target_id                  length  eff_length  est_counts  tpm
ENSMUST00000193812.2       1070    775.000     0.000       0.000000
ENSMUST00000082908.3       110     110.000     0.000       0.000000
ENSMUST00000162897.2       4153    3483.630    1.000       0.013994
ENSMUST00000159265.2       2989    2694.000    0.000       0.000000
ENSMUST00000070533.5       3634    3339.000    0.000       0.000000
```

### Key QC metrics to review

Open `run_info.json` and check:

- **n_processed**: total number of reads processed
- **n_pseudoaligned**: number of reads that pseudoaligned to transcripts
- **p_pseudoaligned**: proportion of reads pseudoaligned (mapping rate)

These metrics provide a first-pass sanity check before summarizing the results to gene-level counts.

::::::::::::::::::::::::::::::::::::::: callout

## What does Kallisto count?

`est_counts` is a model-based estimate, not raw counts of reads. TPM values are *within-sample* normalized and should not be directly compared across samples for statistical testing.

:::::::::::::::::::::::::::::::::::::::

## Step 4: Summarizing to gene-level counts (`tximport`)

Most differential expression tools (DESeq2, edgeR) require **gene-level** counts. We convert transcript estimates → gene-level matrix.

### Prepare tx2gene mapping

We need a mapping file that links transcript IDs to gene IDs. The GTF annotation file contains this information.

```bash
cd $SCRATCH/rnaseq-workshop/data
awk -F'\t' '$3=="transcript" {
    match($9, /transcript_id "([^"]+)"/, tx);
    match($9, /gene_id "([^"]+)"/, gene);
    print tx[1] "\t" gene[1]
}' gencode.vM38.primary_assembly.basic.annotation.gtf > tx2gene.tsv
```

Next, we process this mapping in R. Load the modules and start the R session:

```bash
module load biocontainers
module load r-rnaseq
R
```

In the R session, read in the `tx2gene` mapping and run tximport:

```r
setwd(paste0(Sys.getenv("SCRATCH"), "/rnaseq-workshop"))
library(readr)
library(tximport)

tx2gene <- read_tsv("data/tx2gene.tsv", col_types = "cc", col_names = FALSE)
samples <- read_tsv("scripts/samples.txt", col_names = FALSE)
files <- file.path("results/kallisto_quant", samples$X1, "abundance.h5")
names(files) <- samples$X1

txi <- tximport(files,
                type = "kallisto",
                tx2gene = tx2gene)
saveRDS(txi, file = "results/kallisto_quant/txi.rds")
```

`txi$counts` now contains gene-level counts suitable for DESeq2.

```r
> head(txi$counts)
                      WT_Bcell_mock_rep1 WT_Bcell_mock_rep2 WT_Bcell_mock_rep3 WT_Bcell_mock_rep4
ENSMUSG00000000001.5          473.4426         553.7849         826.9969         676.9101
ENSMUSG00000000003.16           0.0000           0.0000           0.0000           0.0000
ENSMUSG00000000028.16          40.5783          33.9665          55.3055          97.8451
ENSMUSG00000000031.20           0.0000           0.0000           0.0000           0.0000
ENSMUSG00000000037.18          19.2220           9.5997          16.9515          20.2265
ENSMUSG00000000049.12           0.0000           5.0221           1.9846           1.0164
                      WT_Bcell_IR_rep1 WT_Bcell_IR_rep2 WT_Bcell_IR_rep3 WT_Bcell_IR_rep4
ENSMUSG00000000001.5        724.1176        813.8484        689.2341        752.1567
ENSMUSG00000000003.16         0.0000          0.0000          0.0000          0.0000
ENSMUSG00000000028.16        41.6915         63.8123         48.3421         55.7892
ENSMUSG00000000031.20         0.0000          0.0000          0.0000          0.0000
ENSMUSG00000000037.18        19.6807         18.9451         22.3456         17.8923
ENSMUSG00000000049.12         0.9937          9.1869          2.3456          4.5678
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use the default countsFromAbundance setting?

We use the default `countsFromAbundance = "no"` because `DESeqDataSetFromTximport()` (used in Episode 05b) automatically incorporates the `txi$length` matrix to correct for transcript length bias.

Using `"lengthScaledTPM"` would apply length correction twice—once in tximport and again in DESeq2—potentially introducing bias. The default preserves Kallisto's original estimated counts while letting DESeq2 handle length normalization correctly.

**Note:** If you plan to use edgeR instead of DESeq2, you may want `countsFromAbundance = "lengthScaledTPM"` since edgeR's `DGEList` doesn't automatically use the length matrix.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: callout

## Using abundance.h5 vs abundance.tsv

tximport can read either format:

- **abundance.h5**: HDF5 format, includes bootstrap information, faster to read
- **abundance.tsv**: Plain text, easier to inspect manually

We use the `.h5` files here because they load faster and preserve bootstrap data for potential sleuth analysis. If you didn't generate bootstraps, you can use `.tsv` files instead.

:::::::::::::::::::::::::::::::::::::::

## Step 5: QC of quantification metrics

Use MultiQC to aggregate Kallisto results:

```bash
cd $SCRATCH/rnaseq-workshop
module load biocontainers
module load multiqc

multiqc results/kallisto_quant -o results/qc_kallisto
```

Inspect:

* pseudoalignment rates
* fragment length distributions
* consistency across samples

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: examine your Kallisto outputs

Using your MultiQC summary and Kallisto outputs:

1. What is the percent pseudoaligned for each sample, and are any samples outliers?
2. Are the pseudoalignment rates consistent across samples?
3. How similar are the fragment length distributions across samples?
4. Check `run_info.json` for one sample - what parameters were used?

::::::::::::::::::::::::::::::::::::::: solution

Example interpretation for this dataset:

1. All samples show about 85 to 90 percent pseudoaligned. No sample is an outlier and the range is narrow.
2. Pseudoalignment rates are consistent across all samples.
3. Fragment length distributions should be nearly identical across samples.
4. The `run_info.json` should show the index path, number of bootstraps, and thread count used.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::

## Summary

::::::::::::::::::::::::::::::::::::: keypoints

* Kallisto performs fast, alignment-free transcript quantification.
* Salmon is used for strand detection before running Kallisto.
* A transcriptome FASTA is required for building the Kallisto index.
* `kallisto quant` uses FASTQ files directly for pseudo-alignment.
* Transcript-level outputs include TPM and estimated counts.
* `tximport` converts transcript estimates to gene-level counts with `type = "kallisto"`.
* Gene-level counts from Kallisto are suitable for DESeq2.

::::::::::::::::::::::::::::::::::::::::::::::::
