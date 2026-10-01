# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Carpentries Workbench (sandpaper) lesson: "RNA-seq in practice: Hands-on workshop on RCAC systems". It is course material, not software. Learners run the shell steps on RCAC clusters (Negishi/Scholar, SLURM + `module load biocontainers <tool>`) and the R steps in RStudio via Open OnDemand.

## Build and preview

Run from the repo root in R (sandpaper version pinned in `.github/workflows/sandpaper-version.txt`, currently 0.17.1):

```r
sandpaper::build_lesson()   # render into site/ (gitignored except site/README.md)
sandpaper::serve()          # live preview, rebuilds on save
sandpaper::validate_lesson() # check fenced-div structure, links, image alt text
```

There are no tests. "Done" for content edits means `sandpaper::validate_lesson()` is clean and the page renders under `serve()`. CI (`.github/workflows/sandpaper-main.yaml`) builds and deploys on push to `main`; the other workflows are stock Carpentries files and should not be hand-edited.

## Structure that is not obvious from the tree

- `config.yaml` controls episode order and navigation. A new episode file is not shown until it is added to the `episodes:` list. Custom theme: `carpentry: 'rcac'`, `varnish: aseetharam/varnish`.
- Two parallel analysis tracks share episodes 01-03 and 06:
  - Genome track: `04a` (STAR + featureCounts) then `05` (DESeq2 on featureCounts).
  - Transcript track: `04b` (Kallisto + tximport) then `05b` (DESeq2 on tximport).
  Changes to sample names, directory layout, or file paths in one track usually need the matching change in the other, and in `06`, which consumes DE results.
- Episodes are `.Rmd` but R code is not executed at build time. Analysis code is in plain ```` ```r ```` / ```` ```bash ```` fences (display only). The only knitr chunks are `echo=FALSE` chunks calling `knitr::include_graphics("fig/...")`. Figures under `episodes/fig/<episode-topic>/` are pre-rendered screenshots/plots; if example code changes its output, the PNG must be regenerated manually. Do not convert display fences to `{r}` chunks, since the build machine has no data or Bioconductor stack.
- Lesson paths assume the workshop data layout under `/scratch/negishi/aseethar/rnaseq-workshop` (copied to `${RCAC_SCRATCH}` by learners, see `learners/setup.md`), and SLURM examples use `--account=rcac-rnaseq --qos=standby --partition=cpu`. Keep these consistent across episodes.
- Episode markup uses Workbench fenced divs (`questions`, `objectives`, `callout`, `challenge` with nested `solution`, `discussion`, `spoiler`, `keypoints`). Every episode needs front matter with `title`, `teaching`, `exercises`.
- `instructors/instructor-notes.md` is the published instructor page. `instructor-notes/` (top level) holds per-episode narration scripts and a timing table; it is gitignored, local only. If episode content or `teaching`/`exercises` times change, the matching narration file and the timing table in `instructor-notes/README.md` may need updating too.
