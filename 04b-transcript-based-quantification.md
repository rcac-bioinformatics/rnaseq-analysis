---
source: Rmd
title: "4B. Transcript-based quantification (Kallisto)"
teaching: 30
exercises: 30
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

- **Speed**: Kallisto is fast; on this dataset one sample takes 1 to 2 minutes on 16 cores (without bootstraps).
- **Accuracy**: Comparable accuracy to alignment-based methods for quantification.
- **Bootstrap support**: Optional uncertainty estimation via bootstraps (needed for sleuth, not for DESeq2).
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

Building the index for the mouse transcriptome needs more memory than a small job provides, so we run it as a batch job, like the STAR index in Episode 4A. This version of Kallisto builds the index with a single thread, but on Negishi memory comes with the cores you request, so the script asks for 48 cores. Save it as `$SCRATCH/rnaseq-workshop/scripts/index_kallisto.sh`:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --job-name=kallisto_index
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load kallisto

cd $SCRATCH/rnaseq-workshop
mkdir -p data/kallisto_index

kallisto index \
    -i data/kallisto_index/transcripts.idx \
    data/gencode.vM38.transcripts-clean.fa
```

Submit it from the `scripts` directory:

```bash
cd $SCRATCH/rnaseq-workshop/scripts
sbatch index_kallisto.sh
```

This creates an index file that Kallisto uses for pseudo-alignment. Building it for the GENCODE mouse transcriptome takes about 5 minutes and about 7.5 GB of memory. In our tests a 4-core job, which gets about the same amount of memory on Negishi, was killed for running out of it. Check progress with `squeue -u $USER`; the log is written to `cluster-kallisto_index.<jobid>.out` in the `scripts` directory. Wait for the job to finish before quantifying.

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

For one sample, the command looks like this. You do not need to run it yourself: the array job in Step 3 runs it for every sample.

```bash
kallisto quant \
    -i data/kallisto_index/transcripts.idx \
    -o results/kallisto_quant/WT_Bcell_mock_rep1 \
    -b 0 \
    -t 16 \
    data/WT_Bcell_mock_rep1_R1.fastq.gz \
    data/WT_Bcell_mock_rep1_R2.fastq.gz
```

Important flags:

* `-i` → path to the Kallisto index
* `-o` → output directory for this sample
* `-b 0` → number of bootstrap samples (none here; see the callout below)
* `-t 16` → number of threads to use (the array job uses the cores it requests, `${SLURM_CPUS_ON_NODE}`)

::::::::::::::::::::::::::::::::::::::: callout

## Strand-specific libraries

For **strand-specific** libraries, add the appropriate flag:

- `--rf-stranded` for reverse-stranded libraries (e.g., Illumina TruSeq stranded)
- `--fr-stranded` for forward-stranded libraries

Our dataset is **unstranded** (Salmon reported `IU` in Episode 4A), so we omit these flags.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: callout

## Do you need bootstraps?

`-b N` makes Kallisto repeat the quantification on N resampled versions of the reads. The spread of the N estimates measures how uncertain each transcript's abundance is, which matters mostly when reads are shared between similar isoforms.

- **Needed** for sleuth, which uses the bootstrap variance in its model of transcript-level differential expression, and for other tools that model quantification uncertainty (for example fishpond/swish).
- **Not needed** for tximport plus DESeq2, as in this lesson. tximport sums transcripts to genes, where isoform ambiguity largely cancels out, and DESeq2 estimates variability from biological replicates and ignores bootstraps.

Each bootstrap round repeats the expectation-maximization step, so `-b 100` makes every sample several times slower and can push a job past its time limit on `standby`. We use `-b 0`. Use `-b 30` to `-b 100` if you plan a sleuth analysis.

:::::::::::::::::::::::::::::::::::::::

## Step 3: Running Kallisto for many samples

We use a SLURM array job for efficiency and consistency.

First, create a sample list (or reuse the one from Episode 4A):

```bash
cd $SCRATCH/rnaseq-workshop/data
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > $SCRATCH/rnaseq-workshop/scripts/samples.txt
```

Save the array job below as `$SCRATCH/rnaseq-workshop/scripts/quant_kallisto.sh`. Like `map_reads.sh` in Episode 4A, it reads `samples.txt` from the directory you submit from.

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=rcac-rnaseq
#SBATCH --qos=standby
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
    -b 0 \
    -t ${SLURM_CPUS_ON_NODE} \
    $R1 $R2 &> $OUT/$SAMPLE/${SAMPLE}.log
```

Submit from the `scripts` directory:

```bash
cd $SCRATCH/rnaseq-workshop/scripts
sbatch quant_kallisto.sh
```

Each sample takes 1 to 2 minutes on 16 cores.

::::::::::::::::::::::::::::::::::::::: discussion

## Why use an array job here?

Kallisto runs very quickly and produces small output. The advantage of array jobs is *consistency*: all samples are quantified with identical parameters and metadata. This ensures that transcript-level TPM and count estimates are directly comparable across samples.

:::::::::::::::::::::::::::::::::::::::

### Inspect results

After Kallisto finishes, each sample directory contains a small set of output files (plus the `.log` file that our array script writes):

```text
abundance.h5
abundance.tsv
run_info.json
WT_Bcell_mock_rep1.log
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
target_id	length	eff_length	est_counts	tpm
ENSMUST00000193812.2	1070	902.888	0	0
ENSMUST00000082908.3	110	10.4106	0	0
ENSMUST00000162897.2	4153	3985.89	0	0
ENSMUST00000159265.2	2989	2821.89	1	0.0101594
ENSMUST00000070533.5	3634	3466.89	0	0
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

We need a mapping file that links transcript IDs to gene IDs. It must cover every transcript in the Kallisto index, so we build it from the headers of the same transcript FASTA: field 1 is the transcript ID and field 2 the gene ID. The basic GTF lists only 127,934 of these 278,396 transcripts (46 percent) (the transcript FASTA is the comprehensive set, with extra isoforms); with a GTF-based map, tximport would silently drop the rest, together with the reads Kallisto assigned to them.

```bash
cd $SCRATCH/rnaseq-workshop/data
grep ">" gencode.vM38.transcripts.fa | cut -c2- | cut -d "|" -f 1,2 | tr "|" "\t" > tx2gene.tsv
head -3 tx2gene.tsv
```

```text
ENSMUST00000193812.2	ENSMUSG00000102693.2
ENSMUST00000082908.3	ENSMUSG00000064842.3
ENSMUST00000162897.2	ENSMUSG00000051951.6
```

Next, we process this mapping in R. Reading eight Kallisto results takes a few GB of memory, more than a login node allows, so start a short interactive session, then load the modules and start R:

```bash
sinteractive -A rcac-rnaseq -q standby -p cpu -N 1 -n 4 --time=1:00:00
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
head(txi$counts)
```

```text
                      WT_Bcell_IR_rep1 WT_Bcell_IR_rep2 WT_Bcell_IR_rep3
ENSMUSG00000000001.5         387.39815        329.43564        737.97318
ENSMUSG00000000003.16          0.00000          1.00000          0.00000
ENSMUSG00000000028.16         37.73844         23.20617         39.45979
ENSMUSG00000000031.20          2.00000          1.00000          1.00000
ENSMUSG00000000037.18          2.00000          4.00000          3.00000
ENSMUSG00000000049.12          0.00000          0.00000          1.00000
                      WT_Bcell_IR_rep4 WT_Bcell_mock_rep1 WT_Bcell_mock_rep2
ENSMUSG00000000001.5         654.17550           687.1189          620.35065
ENSMUSG00000000003.16          0.00000             0.0000            0.00000
ENSMUSG00000000028.16         29.71743            55.0000           39.20456
ENSMUSG00000000031.20          0.00000             0.0000            0.00000
ENSMUSG00000000037.18          8.00000             1.0000            1.00000
ENSMUSG00000000049.12          0.00000             0.0000            2.00000
                      WT_Bcell_mock_rep3 WT_Bcell_mock_rep4
ENSMUSG00000000001.5            661.3541          704.88506
ENSMUSG00000000003.16             0.0000            0.00000
ENSMUSG00000000028.16            40.0984           41.44387
ENSMUSG00000000031.20             0.0000            0.00000
ENSMUSG00000000037.18             2.0000            1.00000
ENSMUSG00000000049.12             0.0000            1.00000
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use the default countsFromAbundance setting?

We use the default `countsFromAbundance = "no"` because `DESeqDataSetFromTximport()` (used in Episode 5B) automatically incorporates the `txi$length` matrix to correct for transcript length bias.

Using `"lengthScaledTPM"` would apply length correction twice—once in tximport and again in DESeq2—potentially introducing bias. The default preserves Kallisto's original estimated counts while letting DESeq2 handle length normalization correctly.

**Note:** If you plan to use edgeR instead of DESeq2, you may want `countsFromAbundance = "lengthScaledTPM"` since edgeR's `DGEList` doesn't automatically use the length matrix.

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: callout

## Using abundance.h5 vs abundance.tsv

tximport can read either format:

- **abundance.h5**: HDF5 format, includes bootstrap information, faster to read
- **abundance.tsv**: Plain text, easier to inspect manually

We use the `.h5` files here because they load faster. Reading them requires the `rhdf5` package, which the `r-rnaseq` module provides. With `-b 0` the `.h5` and `.tsv` files hold the same estimates, so either works.

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

Interpretation for this dataset:

1. Pseudoalignment rates range from about 51 to 64 percent; none is an outlier. That is lower than the 80 to 90 percent often quoted for well-annotated transcriptomes: Kallisto only counts reads compatible with annotated mature transcripts, so reads from introns (pre-mRNA), intergenic regions, and rRNA not in the transcript FASTA are not pseudoaligned. Compare with the STAR unique mapping rates from Episode 4A, which include intronic reads.
2. Reasonably consistent: they span about 14 percentage points (standard deviation about 5), from mock_rep2 (51 percent) to IR_rep3 and IR_rep4 (64 percent). Every sample had 20 million read pairs processed.
3. The estimated mean fragment length is 159 to 168 bp in the mock samples and 181 to 201 bp in the IR samples. Within a group the samples agree, but the groups differ. This may reflect library preparation done in separate batches. Note such differences: when they line up with the experimental groups, their effect cannot be separated from the treatment.
4. `run_info.json` records the Kallisto and index versions, `n_bootstraps` (0 here), `n_processed`, `n_pseudoaligned`, `p_pseudoaligned`, and the exact command line in `call`.

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
