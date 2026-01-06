---
source: Rmd
title: "B. Transcript-based quantification (Salmon)"
teaching: 30
exercises: 40
---

:::::::::::::::::::::::::::::::::::::: questions

- How do we quantify expression without genome alignment?
- What inputs does Salmon require?
- How do we run Salmon for paired-end RNA-seq data?
- How do we interpret transcript-level outputs (TPM, NumReads)?
- How do we summarize transcripts to gene-level counts using tximport?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Build a transcriptome index for Salmon.
- Quantify transcript abundance directly from FASTQ files.
- Understand key Salmon output files.
- Summarize transcript-level estimates to gene-level counts.
- Prepare counts for downstream differential expression.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

In the previous episode, we mapped reads to the **genome** using STAR and generated **gene-level counts** with featureCounts.

In this episode, we use an alternative quantification strategy:

1. Build a *transcriptome* index  
2. Quantify expression directly from FASTQ files using **Salmon** (*without alignment*)  
3. Summarize transcript estimates back to gene-level counts

This workflow is faster, uses less storage, and models transcript-level uncertainty.

::::::::::::::::::::::::::::::::::::::: callout

## When should I use transcript-based quantification?

Salmon is ideal when:

- You do not need alignment files.
- You want transcript-level TPMs.
- You want fast quantification with bias correction.
- Storage is limited (no BAM files generated).

It is *not* ideal if you need splice junctions, variant calling, or visualization in IGV (all those depends on alignment files).

:::::::::::::::::::::::::::::::::::::::

## Step 1: Preparing the transcriptome reference

Salmon requires a **transcriptome FASTA file** (all annotated transcripts). We previously downloaded transcripts file from GENCODE for this purpose (`gencode.vM38.transcripts.fa`).



::::::::::::::::::::::::::::::::::::::: callout

## Why transcriptome choice matters

Transcript-level quantification **inherits all assumptions of the annotation**.  
Missing or incorrect transcripts → biased TPM estimates.

:::::::::::::::::::::::::::::::::::::::

### Building the Salmon index

A Salmon index is lightweight and quick to build. _You may have already built this in Episode 4A._

```bash
sinteractive -A workshop -p cpu -N 1 -n 4 --time=1:00:00

cd $SCRATCH/rnaseq-workshop
mkdir -p data/salmon_index

module load biocontainers
module load salmon

salmon index \
    --transcripts data/gencode.vM38.transcripts-clean.fa \
    --index data/salmon_index \
    --threads 4
```

This creates a directory containing Salmon's k-mer index of all transcripts.

::::::::::::::::::::::::::::::::::::::: spoiler

## What will the index contain?

```
salmon_index/
├── complete_ref_lens.bin        # IMPORTANT: full transcript lengths used in TPM/effLength calculations
├── ctable.bin                   # internal k-mer lookup table (not user-relevant)
├── ctg_offsets.bin              # internal offsets for reference sequences
├── duplicate_clusters.tsv       # internal duplicate k-mer cluster info
├── info.json                    # IMPORTANT: index metadata (k-mer size, transcripts, build parameters)
├── mphf.bin                     # minimal perfect hash function for fast lookup
├── pos.bin                      # internal k-mer position table
├── pre_indexing.log             # IMPORTANT: log of preprocessing and FASTA parsing steps
├── rank.bin                     # internal rank/select bitvector structure
├── refAccumLengths.bin          # cumulative reference lengths (internal)
├── ref_indexing.log             # IMPORTANT: main indexing log; good for debugging build problems
├── reflengths.bin               # IMPORTANT: transcript lengths used in quantification math
├── refseq.bin                   # IMPORTANT: transcript sequences stored in binary format (core index content)
├── seq.bin                      # internal sequence encoding structures
└── versionInfo.json             # IMPORTANT: Salmon version + environment info for reproducibility
```

:::::::::::::::::::::::::::::::::::::::

## Step 2: Quantifying transcript abundances

The core Salmon command is `salmon quant`.
We use the transcript index and paired-end FASTQ files.

```bash
mkdir -p $SCRATCH/rnaseq-workshop/results/salmon_quant
```

Example command:

```bash
salmon quant \
    --index salmon_index \
    --libType IU \
    --mates1 data/SRR33253286_1.fastq.gz \
    --mates2 data/SRR33253286_2.fastq.gz \
    --output salmon_quant/SRR33253286 \
    --threads 16 \
    --validateMappings \
    --seqBias \
    --gcBias
```

Important flags:

* `--libType IU` → **unstranded**, inward-oriented paired-end
  (Replace with `ISR`, `ISF`, etc. based on your strandness.)
* `--validateMappings` → improves mapping accuracy
* `--seqBias`, `--gcBias` → corrects sequence- and GC-related biases

::::::::::::::::::::::::::::::::::::::: callout

## Correct library type is essential

A wrong `--libType` silently distorts transcript abundances.
Use Salmon's *strandness detection* from Episode 4A if unsure.

:::::::::::::::::::::::::::::::::::::::

## Step 3: Running Salmon for many samples

Again, we use a SLURM array.

First, create a sample list (or reuse the one from Episode 4A):

```bash
cd $SCRATCH/rnaseq-workshop/data
ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > $SCRATCH/rnaseq-workshop/scripts/samples.txt
```

Array job:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=workshop
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --job-name=salmon_quant
#SBATCH --array=1-6
#SBATCH --output=cluster-%x.%j.out
#SBATCH --error=cluster-%x.%j.err

module load biocontainers
module load salmon

DATA="$SCRATCH/rnaseq-workshop/data"
INDEX="$SCRATCH/rnaseq-workshop/data/salmon_index"
OUT="$SCRATCH/rnaseq-workshop/results/salmon_quant"

mkdir -p $OUT

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

R1=${DATA}/${SAMPLE}_1.fastq.gz
R2=${DATA}/${SAMPLE}_2.fastq.gz

salmon quant \
  --index $INDEX \
  --libType IU \
  --mates1 $R1 \
  --mates2 $R2 \
  --output $OUT/$SAMPLE \
  --threads ${SLURM_CPUS_ON_NODE} \
  --validateMappings \
  --seqBias \
  --gcBias
```

Submit:

```bash
sbatch quant_salmon.sh
```
The typical run time is ~5 minutes per sample.

::::::::::::::::::::::::::::::::::::::: discussion

## Why use an array job here?

Unlike the STAR workflow in 4A, Salmon runs quickly and produces small output.
The advantage of array jobs here is not speed, but *consistency*, all samples are quantified with
identical parameters and metadata. This ensures that transcript-level TPM and count estimates
are directly comparable across samples.

:::::::::::::::::::::::::::::::::::::::



### Inspect results

After Salmon finishes, each sample directory contains a small set of output files.  
Only a few of these are important for routine analysis:

```
SRR33253285/
├── quant.sf                       # IMPORTANT: transcript-level TPM, counts, effective lengths
├── lib_format_counts.json         # IMPORTANT: inferred library type (IU, ISR, ISF)
├── aux_info
│   ├── meta_info.json             # IMPORTANT: mapping rate, fragment counts, bias flags
│   ├── fld.gz                     # fragment length distribution (binary format)
│   ├── expected_bias.gz           # bias model used for correction
│   ├── observed_bias.gz           # observed sequence and GC bias
│   └── ambig_info.tsv             # ambiguous equivalence class summaries (advanced QC only)
├── libParams
│   └── flenDist.txt               # fragment length distribution in a human-readable format
├── cmd_info.json                  # command and parameters used to run Salmon
└── logs
└── salmon_quant.log               # main log file for troubleshooting
```

### What is inside *quant.sf*?

This file contains the transcript-level quantification results:

* **Name**: transcript ID  
* **Length**: transcript length  
* **EffectiveLength**: length adjusted for fragment distribution and bias  
* **TPM**: within-sample normalized abundance  
* **NumReads**: estimated number of reads assigned to that transcript  


Example:

```text
Name                       Length  EffectiveLength TPM             NumReads
ENSMUST00000193812.2       1070    775.000         0.000000        0.000
ENSMUST00000082908.3       110     110.000         0.000000        0.000
ENSMUST00000162897.2       4153    3483.630        0.013994        1.000
ENSMUST00000159265.2       2989    2694.000        0.000000        0.000
ENSMUST00000070533.5       3634    3339.000        0.000000        0.000
ENSMUST00000192857.2       480     185.000         0.000000        0.000
ENSMUST00000195335.2       2819    2524.000        0.000000        0.000
ENSMUST00000192336.2       2233    1938.000        0.000000        0.000
ENSMUST00000194099.2       2309    2014.000        0.000000        0.000
```

### Key QC metrics to review

Open `aux_info/meta_info.json` and check:

- **percent_mapped**: proportion of fragments assigned to transcripts
- **library_types**: Salmon inferred library orientation, which should match your `--libType`
- **fragment length statistics**: mean and standard deviation of inferred fragment lengths

These metrics provide a first-pass sanity check before summarizing the results to gene-level counts.


::::::::::::::::::::::::::::::::::::::: callout

## What does Salmon count?

`NumReads` is a model-based estimate, not raw counts of reads. 
TPM values are *within-sample* normalized and should not be directly compared across samples.

:::::::::::::::::::::::::::::::::::::::

## Step 4: Summarizing to gene-level counts (`tximport`)

Most differential expression tools (DESeq2, edgeR) require **gene-level** counts.
We convert transcript estimates → gene-level matrix.

### Prepare tx2gene mapping

We will need a mapping file that links the transcript IDs to gene IDs. The `gencode.vM38.primary_assembly.basic.annotation.gtf` file contains this information.

```bash
cd $SCRATCH/rnaseq-workshop/data
awk -F'\t' '$3=="transcript" {
    match($9, /transcript_id "([^"]+)"/, tx);
    match($9, /gene_id "([^"]+)"/, gene);
    print tx[1] "\t" gene[1]
}' gencode.vM38.primary_assembly.basic.annotation.gtf > tx2gene.tsv
```

Next, we need to process this mapping in R. Load the modules and start the R session, first:

```bash
module load biocontainers
module load r-rnaseq
R
```

In the R session, read in the `tx2gene` mapping:

```r
setwd(paste0(Sys.getenv("SCRATCH"),"/rnaseq-workshop"))
library(readr)
library(tximport)

tx2gene <- read_tsv("data/tx2gene.tsv", col_types = "cc", col_names = FALSE)
samples <- read_tsv("scripts/samples.txt", col_names = FALSE)
files <- file.path("results/salmon_quant", samples$X1, "quant.sf")
names(files) <- samples$X1

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene)
saveRDS(txi, file = "results/salmon_quant/txi.rds")
```

`txi$counts` now contains gene-level counts suitable for DESeq2.

```r
> head(txi$counts)
                      SRR33253285 SRR33253286 SRR33253287 SRR33253288
ENSMUSG00000000001.5     473.4426  553.784929  826.996853  676.910126
ENSMUSG00000000003.16      0.0000    0.000000    0.000000    0.000000
ENSMUSG00000000028.16     40.5783   33.966518   55.305501   97.845071
ENSMUSG00000000031.20      0.0000    0.000000    0.000000    0.000000
ENSMUSG00000000037.18     19.2220    9.599673   16.951493   20.226549
ENSMUSG00000000049.12      0.0000    5.022073    1.984551    1.016411
                      SRR33253289 SRR33253290
ENSMUSG00000000001.5   724.117586  813.848375
ENSMUSG00000000003.16    0.000000    0.000000
ENSMUSG00000000028.16   41.691461   63.812279
ENSMUSG00000000031.20    0.000000    0.000000
ENSMUSG00000000037.18   19.680740   18.945058
ENSMUSG00000000049.12    0.993724    9.186882
```

::::::::::::::::::::::::::::::::::::::: callout

## Why use the default countsFromAbundance setting?

We use the default `countsFromAbundance = "no"` because `DESeqDataSetFromTximport()` (used in Episode 05b) automatically incorporates the `txi$length` matrix to correct for transcript length bias.

Using `"lengthScaledTPM"` would apply length correction twice—once in tximport and again in DESeq2—potentially introducing bias. The default preserves Salmon's original estimated counts while letting DESeq2 handle length normalization correctly.

**Note:** If you plan to use edgeR instead of DESeq2, you may want `countsFromAbundance = "lengthScaledTPM"` since edgeR's `DGEList` doesn't automatically use the length matrix.

:::::::::::::::::::::::::::::::::::::::

## Step 5: QC of quantification metrics

Use MultiQC to verify:

```bash
cd $SCRATCH/rnaseq-workshop
module load biocontainers
module load multiqc

multiqc results/salmon_quant -o results/qc_salmon
```

Inspect:

* mapping rates
* fragment length distributions
* library type consistency

::::::::::::::::::::::::::::::::::::::: challenge

## Exercise: examine your Salmon outputs

Using your MultiQC summary and Salmon fragment length plots:

1. What is the percent aligned for each sample, and are any samples outliers?
2. Is the inferred library type consistent across samples?
3. How similar are the fragment length distributions across samples?
4. Does the fragment length curve show unusual peaks or multimodality?
5. Do CFR and M bias values suggest any strand or positional bias?

::::::::::::::::::::::::::::::::::: solution

Example interpretation for this dataset:

1. All samples show about 89 to 91 percent aligned. No sample is an outlier and the range is narrow.
2. The inferred library type is consistent across samples and matches an unstranded protocol.
3. Fragment length distributions are nearly identical. All samples peak around the same length and have similar spread.
4. The fragment length curves are smooth and unimodal. There are no secondary peaks or irregular shoulders.
5. CFR is 100 percent for all samples and M bias is about 0.5, indicating balanced, strand-neutral coverage without detectable end bias.

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::


## Summary

::::::::::::::::::::::::::::::::::::: keypoints

* Salmon performs fast, alignment-free transcript quantification.
* A transcriptome FASTA is required for building the index.
* `salmon quant` uses FASTQ files directly with bias correction.
* Transcript-level outputs include TPM and estimated counts.
* `tximport` converts transcript estimates to gene-level counts.
* Gene-level counts from Salmon are suitable for DESeq2.

::::::::::::::::::::::::::::::::::::::::::::::::

