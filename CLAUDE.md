# CLAUDE.md

This file guides Claude Code (claude.ai/code) when working in this repository.

## What this is

A Carpentries Workbench (sandpaper) lesson: "RNA-seq in practice: Hands-on workshop on RCAC systems". It is course material, not software. Learners run the shell steps on Purdue RCAC Negishi (SLURM, `module load biocontainers` then the tool module) and the R steps in RStudio via Open OnDemand (`r-rnaseq` module). Canonical repo: `rcac-bioinformatics/rnaseq-analysis`, branch `main`, site built from `gh-pages`. Known issue: `config.yaml` `source:` still points to the retired `aseetharam.github.io/rcac_rnaseq_workshop` URL (tracked in updates.md).

## Dataset and validation anchors

- Data: GEO GSE71176, p53-mediated ionizing radiation response in mouse B cells (Tonelli et al. 2015). 8 paired-end samples (16 FASTQ files), subsampled to 20 M read pairs each with `seqtk sample -s 42`, so reruns are deterministic and comparable to the reference results.
- Reference: GENCODE vM38 on GRCm39 (primary assembly genome FASTA, basic annotation GTF, cleaned transcripts FASTA), fetched from ftp.ebi.ac.uk in episode 02.
- A correct end-to-end run shows IR vs mock upregulation of Cdkn1a, Mdm2, Bax, Bbc3, Pmaip1, Gadd45a, and enrichment of p53 signaling, apoptosis, cell cycle, DNA damage response. Use these as pass/fail anchors when testing; their absence means an upstream break, not "different but fine" results.

## Build and preview

Run from the repo root in R (sandpaper pinned in `.github/workflows/sandpaper-version.txt`, currently 0.17.1):

```r
sandpaper::build_lesson()    # render into site/ (gitignored except site/README.md)
sandpaper::serve()           # live preview, rebuilds on save
sandpaper::validate_lesson() # fenced-div structure, links, image alt text
```

There are no unit tests. "Done" for a content edit means `validate_lesson()` is clean and the page renders under `serve()`. CI (`.github/workflows/sandpaper-main.yaml`) builds and deploys on push to `main`; the other workflows are stock Carpentries files, do not hand-edit them.

## Structure that is not obvious from the tree

- `config.yaml` controls episode order and navigation. A new episode file is invisible until added to the `episodes:` list. Custom theme: `carpentry: 'rcac'`, `varnish: aseetharam/varnish`.
- Two parallel analysis tracks share episodes 01-03 and 06:
  - Genome track: `04a` (STAR + featureCounts) then `05` (DESeq2 on featureCounts).
  - Transcript track: `04b` (Kallisto + tximport) then `05b` (DESeq2 on tximport).
  A change to sample names, directory layout, or file paths in one track usually needs the matching change in the other and in `06`, which must run from either track's DE output.
- Episodes are `.Rmd` but analysis code is NOT executed at build time. It sits in plain ```` ```r ```` and ```` ```bash ```` display fences. The only knitr chunks are `echo=FALSE` calls to `knitr::include_graphics("fig/...")`. Consequences: (a) `knitr::purl()` extracts nothing, so to run or test episode code you must parse the fences out of the Rmd and classify them by language; (b) figures under `episodes/fig/<topic>/` are pre-rendered screenshots and plots, regenerated manually whenever code output changes; (c) never convert display fences to `{r}` chunks, the build machine has no data or Bioconductor stack.
- Lesson paths assume workshop data staged at `/scratch/negishi/aseethar/rnaseq-workshop` (learners rsync it to `${RCAC_SCRATCH}`; a completed copy lives at `/scratch/negishi/aseethar/rnaseq-workshop_results`, see `learners/setup.md`). SLURM examples use `--account=rcac-rnaseq --qos=standby --partition=cpu`; keep these identical across episodes. Negishi scratch purges inactive files, so staged data must be re-verified, and usually re-staged from Fortress or Depot, before every delivery.
- Episode markup uses Workbench fenced divs (`questions`, `objectives`, `callout`, `challenge` with nested `solution`, `discussion`, `spoiler`, `keypoints`). Every episode needs front matter with `title`, `teaching`, `exercises`.

## Delivery format

- Front matter totals ~575 min (~9.5 h) of teaching plus exercises. A delivery day gives ~5.5-6 h of instruction (9:00-16:00 with breaks and lunch), so live delivery is selective while the full lesson stays published for self-paced completion.
- Pattern from the 2026-01-22 delivery, planned again for 2026-10-06: 01 condensed, 02-03 condensed, 04a hands-on with 04b covered conceptually only, 05 hands-on, 06 as introduction with pointers. Fully self-paced: 04b, 05b, the rest of 06. The current split is recorded in updates.md.
- All episodes, both tracks, must keep working and be tested before a delivery, including the self-paced ones.
- `learners/setup.md` carries the day schedule (with its date) and the instructor list; update both for every delivery.

## Instructor materials (local only, gitignored)

- `instructor-notes/`: per-episode narration scripts plus a timing table in its README.md. Near-verbatim delivery scripts.
- `coaching/`: TTS-friendly prep text (plain `.txt`, flowing prose, no markup) that the instructor converts to audio and listens to before a delivery. Coaching is NOT narration; it is a deep conceptual briefing (background, anticipated questions, stumbling points, timing) so the instructor can teach without notes. Generation conventions live in `prompt-02-generate-coaching.md`. Add `coaching/` to `.gitignore` when created.
- If episode content or `teaching`/`exercises` times change, update the matching narration file, the timing table in `instructor-notes/README.md`, and regenerate the affected coaching file.

## Maintenance log

`updates.md` at the repo root is the running log of testing findings and revision opportunities, tagged `[blocker]`, `[should-fix]`, `[nice]`, or `[check]`. Log every content issue there before or alongside fixing it; do not silently fix and move on.
