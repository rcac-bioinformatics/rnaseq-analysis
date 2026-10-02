# Adaptations: where the kit does not run episode code exactly as a learner would

The rule is learner fidelity: episode commands run verbatim, with the same module loads, paths, flags, and SLURM headers, in the real Negishi environment. Every block the kit runs is extracted from the episodes by `build_kit.py`; nothing is hand-copied. This file lists every place where running unattended forces a difference, and why. If an episode command fails as written, the run records the failure; the kit never works around it.

## Where the learner directory lives

| What | Learner | Kit | Why |
|---|---|---|---|
| `$SCRATCH`, `$RCAC_SCRATCH` | `/scratch/negishi/<user>` | `$KIT_SCRATCH_BASE/runs/<RUN_ID>/scratch` | For `aseethar`, `$RCAC_SCRATCH/rnaseq-workshop` *is* the staged source data. Pointing both variables at a fresh test root makes every `cd $SCRATCH/rnaseq-workshop` land in a fresh copy and guarantees the staged data is only read. `lib.sh` refuses to run if the learner directory resolves inside `$STAGED` or `$STAGED_RESULTS`. |
| `$USER` in literal paths | user name | `aseethar/rnaseq-kit/runs/<RUN_ID>/scratch` | Two places build the path from `/scratch/negishi/$USER` instead of `$SCRATCH`: the clean-count-matrix `sed` in 04a and `work_dir` in 05, 05b, and 06. Setting `USER` to the test root's path below `/scratch/negishi/` makes the unchanged code resolve to the test copy. In 04a it is set only around that one block. For a real learner both forms give the same path (preflight checks `$SCRATCH = /scratch/negishi/$USER`). |

## Interactive sessions become batch jobs

| Episode block | Learner | Kit |
|---|---|---|
| 03 FastQC, 04a Salmon, 04b tximport | `sinteractive -A ... -n 4 --time=...`, then commands | A batch job with the same account, QoS, partition, nodes, tasks, and time, parsed from the episode's `sinteractive` line. Lines before `sinteractive` (login node for a learner) run at the start of the same job. |
| 04b `R` (interactive, `r-rnaseq` module) | types the R blocks into R | `R --no-save --no-restore < rsessions/04b-tximport.R`: the episode's R blocks fed on standard input, as if typed. R stops at the first error, as `Rscript` does. |
| Login-node commands: MultiQC (03, 04a x2, 04b), clean count matrix (04a), `tx2gene` (04b) | login node | Short batch jobs (`04a-mapqc`, `04a-post`, `04b-qc`, start of `04b-tximport`), because they depend on earlier jobs finishing. |
| Login-node commands with no job dependency: `ls *_R1.fastq.gz ... > samples.txt` (04a, 04b), `mkdir -p .../results/counts` (04a), `mkdir -p .../results/kallisto_quant` (04b) | login node | Run by `submit_all.sh` on the login node just before the dependent `sbatch`. |

Learner-like shell: the generated sessions run in a child `bash` without `set -e`, `set -u`, or `pipefail` (as in a terminal), with an `ERR` trap that logs every failing command as `### KIT-ERR`. The step fails if any command failed. Lmod's `module` function is used as exported by `sbatch`; if it is not exported, the session sources Lmod's init file, which gives the same function a login shell has.

## The learner's own SLURM scripts

`index_genome.sh`, `map_reads.sh`, `count_features.sh`, `index_kallisto.sh`, and `quant_kallisto.sh` are copied verbatim from the episodes into `scripts/` (the "save it as ..." step) and submitted from `scripts/`, as the episodes say. `submit_all.sh` adds only submission options, never script content: `--parsable` (to read the job ID), `--dependency=afterok:<id>` (to chain the graph), `--kill-on-invalid-dep=yes` (so a failed step cancels its dependents instead of leaving them pending), and `--qos=<standby|normal>` (standby by default, which is what the scripts say).

## Blocks not run, or run differently

| Episode block | Kit | Why |
|---|---|---|
| 02 `fasterq-dump` download | not run | Marked "DO NOT RUN" (hours). Preflight tests compute-node access to NCBI SRA instead. |
| 02 reference download + header cleanup | run verbatim in a separate fresh directory (`$RUN/ep02_fresh`, by pointing `SCRATCH` there for that block) | The learner copy already has the uncompressed files (the episode's new callout says so). The run checks the links still work and compares the downloads with the staged copies (md5). |
| 02 subsampling spoiler | run on 10,000-read copies of the staged FASTQs, with `20000000` replaced by `5000` | Tests the fixed archive-and-rename logic with the real seqtk module in seconds. Workshop-only reference code. |
| 03 fastp spoiler | placeholder names `SRRXXXXXXX_*` replaced by `WT_Bcell_IR_rep1` and kit output paths | The block cannot run as written; the run records the fastp JSON the collect step asks for. |
| 04a STAR exercise solution | not run | Placeholder `SRR1234567`; the mapping array runs the same command for real. |
| 04b single-sample `kallisto quant` example | not run | The episode says not to run it; the array job runs the same command for every sample. |
| 04a prebuilt-index callout (`ln -s ...rnaseq-workshop_results/data/star_index`) | run only with `submit_all.sh --prebuilt-index`, instead of building | Default is to build, so the run measures index time. |
| 05 `design = ~ batch + condition` | not run | A fragment inside a callout, not runnable code. |
| 05 biomaRt spoiler | run in its own R session (`05-biomart`) after the setup blocks | Optional reading, and Ensembl is often unavailable. Kept out of the main 05 session so a BioMart outage cannot change the main results; its success, failure, and output (compared with `data/mart.tsv` in `metrics`) are reported separately. A spoiler failure is a WARN, not a FAIL. |

## R episodes in the Open OnDemand image

Each R step is a fresh `apptainer exec --cleanenv` of the OOD R 4.4.0 image with the host library bound at `/opt/R/host-site-library`, a new empty `HOME`, and a new empty `R_LIBS_USER` (a new learner). If preflight finds that R does not put the host library on `.libPaths()` by itself, the jobs also set `R_LIBS_SITE` to it and preflight reports a WARN so the OOD app's behavior can be confirmed. Resources: 4 cores, standby, as in the setup.md form.

`R/run_blocks.R` runs the blocks in order with `source(echo = TRUE, print.eval = TRUE)`, which prints commands and values as the console does. Differences from pasting into RStudio: an error stops the rest of that block (pasting would continue line by line), and the run then continues with the next block so later failures are also visible; warnings are printed where they occur; each block's plots go to PNG files (1400 x 1000 px, 150 dpi) instead of the Plots pane.

06 on the transcript track swaps the two `read_tsv()` lines exactly as the episode tells Kallisto-track learners, and runs in its own small directory (`scratch_t/rnaseq-workshop`, with `results/deseq2_kallisto` linked from the main copy) so its `results/enrichment` does not overwrite the genome-track one.

## Kit-only additions (never change learner files)

Marked `# KIT EXTRA` in the generated sessions: directory listings for output placeholders, a MultiQC `--export` of the FastQC plots into the kit's records (for regenerating `fig/02_qc/`), the md5 comparison of references, and the subsample check. The `metrics` step (`R/metrics.R`) reads saved outputs to report anchors, PCA variance, library sizes, and the genome-vs-transcript comparison.
