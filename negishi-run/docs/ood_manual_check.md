# Open OnDemand manual check

The kit runs the R episodes in the OOD image from batch jobs. This checklist covers what only a real browser session shows: the form, the launch time, and the RStudio experience. Do it once, after `submit_all.sh` has finished at least the `04a-post` step of a run (you need its count matrix). Fill in the blanks and bring this file back with the tarball.

Date and time: ________  Run ID used: ________  Browser and network: ________

## 1. Launch with the setup.md form values

- [ ] Go to <https://gateway.negishi.rcac.purdue.edu/>, log in.
- [ ] **Interactive Apps**: the app is listed as **RStudio (bioconductor)** under **Bioinformatics Apps** (not RStudio Server under GUIs). Matches `learners/fig/ood_rstudio_dropdown.png`? yes / no
- [ ] Fill the form exactly as in `learners/setup.md`: Partition `cpu`, Account `rcac-rnaseq`, QoS `standby`, Wall Time `4`, Cores `4`, R version `4.4.0-bioconductor`, "Load extra RCAC site library (advanced)" checked. Every field present with these names? yes / no. Differences: ________
- [ ] Screenshot the filled form (to compare with `learners/fig/ood_rstudio_resources.png`).
- [ ] Click **Launch** at time ________ ; status **Running** at ________ ; RStudio usable after **Connect to RStudio Server** at ________ . Launch-to-ready: ______ min.

## 2. Smoke test

In the RStudio console (adjust the path to where you copied the kit):

```r
source("~/rnaseq-analysis/negishi-run/tools/workshop_smoke_test.R")
```

- [ ] Last line says `All packages load`. Missing packages, if any: ________
- [ ] R and Bioconductor versions printed: ________ (expected R 4.4.0, Bioconductor 3.20)
- [ ] `.libPaths()` includes `/opt/R/host-site-library`? yes / no
- [ ] KEGG REST, Ensembl BioMart, MSigDB lines show `200`: ________
- [ ] Number of "package 'X' was built under R version ..." warnings when packages load: ______ (learners will see these)

## 3. Episode 5A interactively, up to the first plot

The episode builds its working directory as `/scratch/negishi/$USER/rnaseq-workshop`. For `aseethar` that is the **staged source data**, so point it at the kit run's learner copy first. Run this line before any episode code (replace `<RUN_ID>`):

```r
Sys.setenv(USER = "aseethar/rnaseq-kit/runs/<RUN_ID>/scratch")
file.path("/scratch/negishi", Sys.getenv("USER"), "rnaseq-workshop")   # must print the kit run path
```

Then open Episode 5A on the published site (or `episodes/05-deseq2-expression-analyses.Rmd`) and paste its R blocks into the console in order, starting at "Setup: load packages and data", skipping the biomaRt spoiler, until the gene biotype bar plot appears.

- [ ] Started pasting at ________ ; first plot (gene biotypes) shown at ________ : ______ min.
- [ ] Any error before the plot? Text: ________
- [ ] Warnings shown (copy the distinct ones): ________
- [ ] The plot matches the kit's `figures/episodes/fig/05_deseq/gene-biotypes.png`? yes / no
- [ ] Memory in use after the plot (Environment pane, memory usage widget): ______ GB

When done: `Sys.setenv(USER = "aseethar")`, then quit the session from the OOD **My Interactive Sessions** page so the job ends.
