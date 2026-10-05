# CLAUDE.md

This file guides Claude Code (claude.ai/code) when working in this repository.

## What this is

A Carpentries Workbench (sandpaper) lesson: "RNA-seq in practice: Hands-on workshop on RCAC systems". It is course material, not software. Learners run the shell steps on Purdue RCAC Negishi (SLURM, `module load biocontainers` then the tool module; 04b's tximport uses the `r-rnaseq` module) and the R steps (05, 05b, 06) in the Open OnDemand app RStudio (bioconductor), R 4.4.0 with Bioconductor 3.20. That app runs `/depot/itap/aseethar/images/rstudio_bioc_ood_rocky8_r4.4.0_s2025.05.0-496_tex.sif` with the host library `/apps/biocontainers/extras/r-ood/r4.4.0_s2025.05.0-496` bound at `/opt/R/host-site-library`. Canonical repo: `rcac-bioinformatics/rnaseq-analysis`, branch `main`, site built from `gh-pages`.

## Dataset and validation anchors

- Data: GEO GSE71176, p53-mediated ionizing radiation response in mouse B cells (Tonelli et al. 2015). 8 paired-end samples (16 FASTQ files), subsampled to 20 M read pairs each with `seqtk sample -s 42`, so reruns are deterministic and comparable to the reference results.
- Sample selection is a design decision, not an error: GSE71176 (SRA SRP061386) has 24 samples (wild-type and Trp53-/- mice; splenic B and non-B cells; mock and 7 Gy IR at 4 h; 4 replicates per wild-type group, 2 per knockout group). The workshop uses only the 8 wild-type B-cell runs: SRR2121778-SRR2121781 (WT_Bcell_mock_rep1-4) and SRR2121786-SRR2121789 (WT_Bcell_IR_rep1-4), for a simple two-group design with 4 replicates and clear positive controls. Excluded: all non-B-cell samples (8 WT) and all Trp53-/- samples (8). Never expand `SRR_Acc_List.txt` or the episodes back to 24; Episode 02 explains the full series in a callout.
- Sequencing batches (from the read names): IR_rep1, IR_rep2, mock_rep1, mock_rep2, mock_rep3 on flowcell C1PF5ACXX (run 160, lanes 6-8); IR_rep3 and mock_rep4 on C1U04ACXX (run 174); IR_rep4 on C1TY0ACXX (run 175). The C1PF5ACXX samples are the lower-quality FastQC group in Episode 03 and form the small within-group structure in the 05/05b PCA. Batch is not confounded with condition. The analysis deliberately uses `~ condition`; see the revision opportunities in the maintenance log before changing that.
- Reference: GENCODE vM38 on GRCm39 (primary assembly genome FASTA, basic annotation GTF, cleaned transcripts FASTA), fetched from ftp.ebi.ac.uk in episode 02.
- A correct end-to-end run shows significant IR vs mock upregulation (padj < 0.05, positive log2 fold change, on both tracks) of Cdkn1a, Mdm2, Bax, Bbc3, and Pmaip1, plus, among upregulated genes on both tracks: KEGG p53 signaling and HALLMARK_P53_PATHWAY (ORA), GO "DNA damage response, signal transduction by p53 class mediator" (ORA of upregulated genes and GSEA), and in GO BP GSEA with positive NES "intrinsic apoptotic signaling pathway in response to DNA damage" (apoptosis) and DNA damage / mitotic cell cycle checkpoint signaling (cell cycle). HALLMARK_APOPTOSIS, G2M, and E2F are not significant at 20 M read pairs (2026-10-02); do not use them as anchors. Use these as pass/fail anchors when testing; their absence or a negative sign means an upstream break (the inverted-contrast bug looked exactly like that), not "different but fine" results.
- Gadd45a is a direction-only anchor: positive log2 fold change but not significant at 20 M read pairs (2026-10-02: +0.45, padj 0.068 genome track; +0.49, padj 0.057 transcript track). Do not treat its non-significance as a failure.
- Tool versions (Negishi biocontainers defaults, 2026-10-02): fastqc 0.12.1, multiqc 1.23, fastp 0.23.2, star 2.7.11b, subread 2.0.1 (`-p` counts read pairs; no `--countReadPairs`), kallisto 0.48.0, salmon 1.10.1, sra-tools 2.11.0-pl5262, seqtk 1.4. OOD R packages: DESeq2 1.46.0, tximport 1.34.0, clusterProfiler 4.14.6, enrichplot 1.26.6, biomaRt 2.62.1, msigdbr 26.1.1.

## Build and preview

Run from the repo root in R (sandpaper pinned in `.github/workflows/sandpaper-version.txt`, currently 0.17.1):

```r
sandpaper::build_lesson()    # render into site/ (gitignored except site/README.md)
sandpaper::serve()           # live preview, rebuilds on save
sandpaper::validate_lesson() # fenced-div structure, links, image alt text
```

There are no unit tests. "Done" for a content edit means `validate_lesson()` is clean and the page renders under `serve()`. Before every delivery, and after any change to episode code, run the Negishi run kit (`negishi-run/README.md`): it regenerates every job from the episodes (`python3 negishi-run/build_kit.py`), runs both tracks as a new learner would on standby, and returns a tarball whose SUMMARY.md has per-step status, the anchors, runtimes, and the values for output blocks. Output blocks, numbers, and figures in the episodes come from that run, never from an off-cluster test. CI (`.github/workflows/sandpaper-main.yaml`) builds and deploys on push to `main`; the other workflows are stock Carpentries files, do not hand-edit them.

## Structure that is not obvious from the tree

- `config.yaml` controls episode order and navigation. A new episode file is invisible until added to the `episodes:` list. Custom theme: `carpentry: 'rcac'`, `varnish: aseetharam/varnish`.
- Two parallel analysis tracks share episodes 01-03 and 06:
  - Genome track: `04a` (STAR + featureCounts) then `05` (DESeq2 on featureCounts).
  - Transcript track: `04b` (Kallisto + tximport) then `05b` (DESeq2 on tximport).
  A change to sample names, directory layout, or file paths in one track usually needs the matching change in the other and in `06`, which must run from either track's DE output.
- Episodes are `.Rmd` but analysis code is NOT executed at build time. It sits in plain ```` ```r ```` and ```` ```bash ```` display fences. The only knitr chunks are `echo=FALSE` calls to `knitr::include_graphics("fig/...")`. Consequences: (a) `knitr::purl()` extracts nothing, so to run or test episode code you must parse the fences out of the Rmd and classify them by language; (b) figures under `episodes/fig/<topic>/` are pre-rendered screenshots and plots, regenerated manually whenever code output changes; (c) never convert display fences to `{r}` chunks, the build machine has no data or Bioconductor stack.
- Lesson paths assume workshop data staged on Depot at `/depot/workshop/data/rnaseq-workshop` (group-readable; restaged with `negishi-run/restage.sh` and checked with `negishi-run/check_staged.sh`; learners rsync it to `${RCAC_SCRATCH}`; a completed copy lives at `/depot/workshop/data/rnaseq-workshop_results`, see `learners/setup.md`). SLURM examples use `--account=rcac-rnaseq --qos=standby --partition=cpu`; keep these identical across episodes. Never set `--mem`: Negishi allocates about 1.9 GB per core, so memory-heavy steps request more cores (STAR and kallisto indexes use 48). The staged copy must contain exactly what `negishi-run/staged_manifest.tsv` lists (FASTQ, references, `SRR_Acc_List.txt`, `tx2gene.tsv`, `mart.tsv`, `annot.tsv`, `scripts/`); learners create indexes, `results/`, and job logs themselves. Depot is not purged, but run `check_staged.sh` before every delivery. The completed-results copy is on Depot too (group-readable, without the FASTQ and reference files, which are in the staged copy); learners copy single files from it, not the whole directory. Regenerate it with `negishi-run/regen_results.sh` whenever lesson code changes its outputs.
- Figures under `episodes/fig/` come in two kinds. Plots of this workshop's data (most of `02_qc/`, `05_deseq/`, `06-enrich/`) are regenerated from the Negishi run kit. Contrast examples are intentional and must not be replaced or "corrected" to this dataset: the FastQC cartoons by Zandra Selina (CC BY 4.0) and plots from other datasets in `fig/02_qc/examples/`, plus `fig/02_qc/fastqc_per_base_sequence_content.png` (cartoon) and `fig/02_qc/fastqc-qual.png` (150 bp dataset). Their captions say they are schematic or from a different dataset.
- Episode markup uses Workbench fenced divs (`questions`, `objectives`, `callout`, `challenge` with nested `solution`, `discussion`, `spoiler`, `keypoints`). Every episode needs front matter with `title`, `teaching`, `exercises`.

## Delivery format

- Front matter totals ~575 min (~9.5 h) of teaching plus exercises. A delivery day gives ~5.5-6 h of instruction (9:00-16:00 with breaks and lunch), so live delivery is selective while the full lesson stays published for self-paced completion.
- Pattern from the 2026-01-22 delivery, planned again for 2026-10-06: 01 condensed, 02-03 condensed, 04a hands-on with 04b covered conceptually only, 05 hands-on, 06 as introduction with pointers. Fully self-paced: 04b, 05b, the rest of 06. The current split is recorded in `maintenance/log/updates.md`.
- All episodes, both tracks, must keep working and be tested before a delivery, including the self-paced ones.
- `learners/setup.md` carries the day schedule (with its date) and the instructor list; update both for every delivery.

## Instructor materials (local only, gitignored)

- `instructor-notes/`: per-episode narration scripts plus a timing table in its README.md. Near-verbatim delivery scripts.
- `coaching/`: TTS-friendly prep text (plain `.txt`, flowing prose, no markup) that the instructor converts to audio and listens to before a delivery. Coaching is NOT narration; it is a deep conceptual briefing (background, anticipated questions, stumbling points, timing) so the instructor can teach without notes. Generation conventions live in `maintenance/prompts/prompt-02-generate-coaching.md`. Add `coaching/` to `.gitignore` when created.
- If episode content or `teaching`/`exercises` times change, update the matching narration file, the timing table in `instructor-notes/README.md`, and regenerate the affected coaching file.

## Maintenance log and what the site publishes

sandpaper renders every `.md`/`.Rmd` file at the repo root and one folder below it (except files named README or CONTRIBUTING), and CI publishes whatever is tracked. Keep non-lesson Markdown out of that reach: this file lives in `.claude/` (hidden folders are skipped), the maintenance log and prompts two levels down in `maintenance/log/` and `maintenance/prompts/`, and the kit docs in `negishi-run/docs/`.


`maintenance/log/updates.md` is the running log of testing findings and revision opportunities, tagged `[blocker]`, `[should-fix]`, `[nice]`, or `[check]`. Log every content issue there before or alongside fixing it; do not silently fix and move on.
