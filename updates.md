# updates.md

Running log of testing findings and revision opportunities for this lesson. Severity tags: `[blocker]` breaks a learner or the live delivery, `[should-fix]` wrong or stale but survivable, `[nice]` polish, `[check]` needs verification before it becomes a finding. Append dated entries; move resolved items to Done with the fixing commit.

## Delivery constraint and selective teaching plan (2026-10-06)

Front matter totals ~575 min of teaching plus exercises; the day allows ~5.5-6 h of instruction. Everything stays published and everything gets tested; only part is taught live.

| Episode | Live on 2026-10-06 | Self-paced |
|---|---|---|
| 01 Intro to RNA-seq | Yes, condensed to ~45 min | Full episode |
| 02 Downloading and organizing | Yes, condensed (shares a ~45 min block with 03) | Full |
| 03 QC and trimming | Yes, condensed | Full |
| 04a STAR + featureCounts | Yes, ~75 min hands-on | n/a |
| 04b Kallisto | Concept only, ~10 min | Full |
| 05 DESeq2 (featureCounts) | Yes, ~75 min hands-on | n/a |
| 05b DESeq2 (Kallisto path) | No | Full |
| 06 Enrichment | Introduction plus pointers, ~30 min | Full |

This mirrors the 2026-01-22 schedule in `learners/setup.md`. Adjust after the timed dry run.

## Open findings

### Static review, 2026-09-25

- `[should-fix]` `learners/setup.md`: schedule heading still says 01/22/2026. Update to 10/06/2026 and confirm the instructor list (Tomas Ratkus, Rose Wilfong) for October.
- `[should-fix]` `config.yaml` `source:` points to `https://aseetharam.github.io/rcac_rnaseq_workshop`; the repo now lives at `rcac-bioinformatics/rnaseq-analysis`. Point to the canonical repo or the published site.
- `[nice]` `index.md` is the stock Workbench placeholder text. Replace with a short landing blurb: what the workshop is, audience, prerequisites, link to setup.
- `[nice]` `README.md` is one line. Add title, site link, license, delivery history.
- `[nice]` `instructor-notes/README.md` timing table lists episode 02 as 37 min; front matter says 20 + 15 = 35. Align.
- `[check]` Staged data at `/scratch/negishi/aseethar/rnaseq-workshop` and `..._results` dates to January; Negishi scratch purge has likely removed some or all of it. Re-stage from Fortress/Depot and verify before any execution testing.
- `[check]` `--account=rcac-rnaseq`: confirm the account still exists and whether October registrants must be added to it, or switch instructions to a workshop reservation or participants' own accounts.
- `[check]` Version drift to confirm during execution testing: GENCODE mouse release links (episode 02 pins vM38 on ftp.ebi.ac.uk; verify they still resolve), and R package API changes in the current `r-rnaseq` stack. Concrete case: episode 06 line ~621 calls `msigdbr(species = "Mus musculus", category = "H")`; newer msigdbr renamed `category`/`subcategory` to `collection`/`subcollection`, so check the installed version and update the call or pin the package.
- [blocker, fixed] tximport and msigdbr missing from the 4.4.0 host library.
- [should-fix] episode 06 msigdbr category= changed to collection=.
- [should-fix] biomaRt fallback to mart.tsv in episode 05.
- [ops] app forwarded host PATH/LD_LIBRARY_PATH; 4.6.0 library has 86 wrong-distro packages, rebuild pending.

### Setup page rewrite and package audit (prompt-03), 2026-10-01

Static findings surfaced while building the setup.md package list. Not yet fixed in the episodes.

- `[blocker]` LFC direction is inverted after shrinkage. `condition` levels sort alphabetically, so `WT_IR` is the reference and the apeglm coefficient is `condition_WT_mock_vs_WT_IR` (mock vs IR), while `results()` uses IR vs mock. Episode 05 then sets `res <- res_shrunk`, so the volcano plot, up/down labels, exported `DESeq2_results_joined.tsv`, and every direction-aware step in 06 (up/down ORA, GSEA NES sign) are flipped relative to the text ("right = higher in IR"). The shown outputs already reflect this: `summary(res_shrunk)` swaps the up/down counts (3260/3317 vs 3317/3260), and 06 GSEA reports "signal transduction in response to DNA damage" with negative NES. p53 targets (Cdkn1a, Mdm2, Bax, Bbc3) would appear as down. Files: `05-deseq2-expression-analyses.Rmd:686-697,736`; `05b-deseq2-expression-analyses-kallisto.Rmd:463-466,517-527`. Fix: `coldata$condition <- relevel(coldata$condition, ref = "WT_mock")` before building `dds` in both 05 and 05b, use `coef = "condition_WT_IR_vs_WT_mock"`, use contrast IR vs mock in 05b, then regenerate printed outputs and figures in 05, 05b, 06.
- `[blocker]` `05-deseq2-expression-analyses.Rmd:904` saves to `results/deseq2_results/dds_featurecounts.rds`, but the episode only creates `results/deseq2` (line 65). `saveRDS()` errors. Fix: save to `results/deseq2/dds_featurecounts.rds`.
- `[should-fix]` 05b contrast and output disagree: `05b:465` requests `c("condition", "WT_mock", "WT_IR")`, but the printed header at `05b:493-494` and the volcano callout at `05b:683` say WT_IR vs WT_mock. Resolved by the relevel fix above.
- `[should-fix]` 05b overwrites the shrunk table: `05b:560` builds `res_df` from `res_shrunk`, then `05b:601` rebuilds it from unshrunk `res`, so the volcano and saved tables use unshrunk LFCs. 05b also saves `res_df` (`05b:718-735`) instead of the annotated `res_annot`. Fix: build `res_annot` from the shrunk table and save that.
- `[should-fix]` Unused, mismatched packages are loaded, which forces learners to install them: `05:98-99` loads `EnsDb.Mmusculus.v79` (Ensembl 79, GRCm38, does not match GRCm39/vM38) and `ensembldb`; `05:106` loads `ComplexHeatmap`; `05b:100` loads `tximportData`. None is used. Fix: remove the four `library()` calls, then drop them from the setup.md install list.
- `[should-fix]` Episode 05 annotation spoiler (`05:152-206`) does not run: `rownames(counts)` should be `rownames(cts)`; `mart %>% select(...)` at line 199 uses `mart` before it exists (should be `annot`); the `write.table(` call is missing its closing parenthesis; the prose says `org.Mm.eg.db` but the code uses biomaRt; stray "0" at the end of line 152. `05b:581` says the annotation comes "from Episode 02"; it comes from this spoiler or the staged `data/mart.tsv`.
- `[should-fix]` Episode 05 OOD instructions are for Scholar: `05:47` and `05:75` say "Scholar", the form fields at `05:80-83` (`queue: rcac-rnaseq`, walltime, cores) do not match the Negishi form (Partition, Account, QoS, Wall Time, Cores, R version) now described in setup.md, and `fig/05_deseq/open-on-demand.png` is a Scholar screenshot. Narration (`instructor-notes/05-...-narration.md`) says `queue: workshop`. Fix: point 05 to setup.md, replace the screenshot with a Negishi one.
- `[should-fix]` `04a-genome-based-quantification.Rmd:64` uses `sinteractive -A rcac-workshop`; every other job uses `rcac-rnaseq`.
- `[should-fix]` `04a:193` and `04a:238` use `$SCRATCH/rcac_rnaseq/...`; the working directory is `rnaseq-workshop`.
- `[should-fix]` Episode 06 msigdbr call (extends the 2026-09-25 `[check]` entry): msigdbr 26.1.1 (current CRAN) still accepts `category = "H"` with a deprecation warning, but the `entrez_gene` column used at `06:645` is only added when the deprecated argument is used. Switching to `collection = "H"` without also changing `entrez_gene` to `ncbi_gene` breaks the episode. Change `06:621` and `06:645` together. msigdbr >= 24.1 also downloads gene set data on first use, so the OOD compute node needs internet egress.
- `[should-fix]` 06 depends on the 05 output path only: `06:161` reads `results/deseq2/DESeq2_results_joined.tsv`, but the prerequisites (`06:40`) say 05 or 05b. Add a line telling 05b learners to read `results/deseq2_kallisto/DESeq2_kallisto_results.tsv` instead.
- `[should-fix]` `04b:304-326` runs tximport in `R` on whatever node the learner is on, with no `sinteractive` line. Reading 8 `abundance.h5` files with 100 bootstraps (the shown `txi` includes `infReps`) is not a login-node task. Add the same `sinteractive` line used in 03.
- `[check]` Scratch variables: episodes use `$SCRATCH` (`02:254`, `03:63`, `04a`, `04b`, and `Sys.getenv("SCRATCH")` at `04b:313`), setup.md uses `$RCAC_SCRATCH`, and the R episodes hardcode `/scratch/negishi/$USER` (`05:110`, `05b:109`, `06:154`). Confirm `$SCRATCH` is set on Negishi (setup.md now asks learners to `echo` both); consider standardizing on one.
- `[check]` `01:91-97` asks learners to launch ConfoundingExplorer. It is GitHub-only (csoneson/ConfoundingExplorer 0.5.1), not on CRAN or Bioconductor 3.22, and the URL in `instructor-notes/README.md` is the package docs site, not a hosted app. Decide whether learners install it in OOD RStudio (`remotes::install_github`) or the instructor demos it.
- `[check]` Runtime packages not loaded by name: tximport needs `rhdf5` to read `abundance.h5` (04b), and vsn `meanSdPlot()` needs `hexbin` (05). Confirm both are in the OOD R stack.
- `[check]` OOD sizing in setup.md: 4 cores (~8 GB; the Negishi form states ~2 GB per core and a 4 h standby cap, per the scRNA-seq workshop's May 2026 form screenshot), 4 h, standby, R `4.4.0-bioconductor`. Record actual MaxRSS of 05/05b/06 in the runtime record and confirm the R version option loads every package.
- `[nice]` `02:119` says "four subdirectories" and lists three.
- `[nice]` `04a:588` keypoint says strandedness "should be inferred using aligned BAM files"; the episode uses Salmon on FASTQ. `04b:152` refers to "Step 0", which does not exist.
- `[nice]` `instructor-notes/README.md` lists `data/salmon_index_strand/`; 04a builds `data/salmon_index`.
- `[nice]` `06:225-230`: `bitr()` returns duplicate Ensembl-to-Entrez rows ("Mapped 11358 of 11330 genes (100.2%)"), so `universe_entrez` has duplicates. Wrap with `unique()`.
- `[nice]` CLAUDE.md says `updates.md` is at the repo root; it lives in `new-prompts/`.

### Execution testing (prompt-01-test-materials.md)

Append here: episode, finding, severity, file:line, proposed fix.

#### Local execution test, 2026-10-01 (prompt-00-local-test.md; full scale, staged 20 M-pair FASTQs, quay.io biocontainers, R 4.4.0 / Bioc 3.19 container; details in local-test/REPORT.md)

Static checks

- `[should-fix]` 02 subsampling spoiler (`02-downloading-and-organizing-files.Rmd:381`): `mv *_R1.fastq.gz *_R2.fastq.gz original_reads/` also matches the new `*_sub_R1.fastq.gz` / `*_sub_R2.fastq.gz`, so every subsampled file is moved into `original_reads/` and the rename loop at 384-388 fails with "No such file". Fix: write subsampled output to a separate directory (`mkdir -p sub; seqtk ... > sub/${sample}_R1.fastq.gz`), then move originals and `mv sub/* .`. Patch: `updates-patches/02-subsample-spoiler.diff`.
- `[nice]` Output and directory listings are fenced as `bash`/`r`, so copy-paste and syntax checkers treat them as code (shellcheck SC2287/SC2211, R `parse()` errors): `02:136-144` (prompt + `tree`), `02:412-435`, `03:72-105`, `04a:241-259`, `04b:330-346` (`> head(txi$counts)` with output in an `r` fence), `06:949-956` (`head(gene_list)` and its output in one `r` fence). Fix: fence output as `text`, keep the command in its own `r`/`bash` fence.
- `[nice]` 13 orphan figures under `episodes/fig/` (not referenced anywhere): `simplifyEnrichment.png`, `geneset.svg`, `msigdb.png`, `02_qc/fastqc-qual.png`, `01_fq2counts/{dir_org.png,geo-db_old.png,Screenshot_20260120_133452.png,srr-runs-old.png}`, `02_setup/{folders_working_dir,data_folder,rstudio_project,rstudio_console}.png`, `06-enrich/gsea-2.png`. Delete or use.
- `[nice]` Permanent (301) redirects in references: 7 `genomebiology/bmcbioinformatics.biomedcentral.com` links now go to `link.springer.com`; `www.ncbi.nlm.nih.gov/pmc/articles/PMC4766705/` to `pmc.ncbi.nlm.nih.gov`; `liebertpub.com/doi/10.1089/omi.2011.0118` to `journals.sagepub.com`. No 404s. 403/429 from academic.oup.com, annualreviews, bioRxiv, cshlp are bot blocking (`[check]` in a browser). `www.putty.org` timed out from here (`[check]`).
- `validate_lesson()` (sandpaper 0.17.1): clean, no output beyond section headers. All `include_graphics()` targets exist.
- shellcheck 0.10.0 on 34 bash blocks: beyond the items above, only style notes (SC2086 unquoted vars x49, SC2164 `cd` without `|| exit` x17). Not worth changing in teaching code.

Execution (shell track)

- `[should-fix]` 03 QC solution (`03-qc-and-trimming.Rmd:296-300`) says "All replicates show high, stable sequence quality" and "no samples show concerning patterns". Actual FastQC (0.12.1) on the staged data: per-base sequence quality FAIL in 12/16 files (lower quartile drops to Q2 by cycle 48), per-tile quality FAIL 12/16, per-base N content FAIL 7/16, per-sequence GC FAIL 8/16; adapter content PASS 16/16. Duplication differs by sample: IR_rep2 and mock_rep3 are ~96% unique vs 52-85% for the others (possible separate run/lane; worth a sentence). Fix: rewrite the solution around the real report (quality drop at read ends, N spikes from early cycles, adapters absent, trimming still not needed for STAR), and regenerate `fig/02_qc/*` and `multiqc-summary.png` from this dataset if they came from another one.
- `[nice]` 03 MultiQC block (`03:238`) starts with `cd rnaseq-workshop` while the learner is already in `$SCRATCH/rnaseq-workshop` from the FastQC block; it errors ("No such file or directory") before the corrective `cd $SCRATCH/rnaseq-workshop`. Delete the first `cd`.
- `[should-fix]` 04a strandedness spoiler (`04a:122`): the prose says ISF/ISR are "36.4M" and "35.2M"; the JSON directly above (which this run reproduced byte for byte: ISF 5,322,565, ISR 5,150,359, bias 0.508, `IU`) says 5.3M and 5.2M. Fix the prose. The strand call (`-s 0`, unstranded) is confirmed.
- `[should-fix]` 04a featureCounts (`04a:507`): `-p` without `--countReadPairs` counts reads, not fragments, in subread >= 2.0.2 (tested 2.0.6: totals ~46 M alignments per sample = 2 x reads; Pmaip1 IR_rep1 1701 reads vs 900 pairs). The prose (`04a:451`) says "how many read pairs map to each gene", and the printed DESeq2 outputs in 05/06 were made with fragment counts (published baseMean 614 vs 1166 here for ENSMUSG00000075122; 11330 vs 12058 genes pass the filter). Fix: add `--countReadPairs` if the Negishi `subread` default is >= 2.0.2 (it errors on older versions; check `module spider subread`). Patch: `updates-patches/04a-account-paths-countreadpairs.diff`.
- `[should-fix]` 04a alignment solution (`04a:437-441`) says "All samples have ~80% uniquely mapped reads. No sample is an outlier." Actual STAR 2.7.11b unique mapping: IR_rep1 61.6%, IR_rep2 63.5%, IR_rep3 77.4%, IR_rep4 78.3%, mock_rep1 70.9%, mock_rep2 68.1%, mock_rep3 68.4%, mock_rep4 83.3% (mean 71%; "too short" unmapped 1.4-4.6 M). Fix: quote the real range and use IR_rep1/2 as the discussion example.
- `[check]` 04a featureCounts assignment is 31-39% of alignments (multimapping 20-26%, no features 21-36%, unmapped counted because of `--outSAMunmapped Within`). The staged `README.md` QC checkpoint says ">60% assignment"; the 04a solution gives no number. Either explain the low figure (denominator includes unmapped and multimapped alignments) or drop the >60% claim from the staged README.
- `[check]` 04a STAR index job has no `--mem`. Measured MaxRSS 32.7 GB (genomeGenerate), 29.2 GB per mapping task. With Negishi's default memory per core, 20 cores should give ~40 GB, but add `#SBATCH --mem=40G` (or confirm the default) so a smaller default cannot OOM-kill the index build.
- `[should-fix]` Staged `scripts/*.sh` in `/scratch/negishi/aseethar/rnaseq-workshop` (copied to learners by setup.md) use `--account=workshop` and no `--qos`; the episodes use `--account=rcac-rnaseq --qos=standby`. A learner who submits the staged copy gets a rejected job. Re-stage with the episode versions.
- `[should-fix]` 04b kallisto runtime: the text (`04b:46`, `04b:101`, `04b:225`) says ~2-3 min per sample and index build 2-3 min, and the array requests `--time=1:00:00`. Measured here (kallisto 0.51.1, `-b 100 -t 16`): index 22 min (contended), quant 30-41 min per sample with 100 bootstraps on 16 threads (first two samples contended with STAR; uncontended times in REPORT.md). On standby a 1 h limit is at risk. Fix: use `-b 0` or `-b 10` for the workshop (bootstraps are only needed for sleuth, as the callout says), or raise `--time`; correct the runtime claims.
- `[should-fix]` 04b solution (`04b:406`) says "about 85 to 90 percent pseudoaligned". Actual `p_pseudoaligned`: see REPORT.md (mock_rep1 59.5%, IR_rep1 52.7%, IR_rep2 53.3%). Fix the number.
- `[should-fix]` 04b array script reads `samples.txt` by relative path (`04b:204`) and the episode never says where to save `quant_kallisto.sh` or where to run `sbatch` from; the preceding block leaves the learner in `data/`. Same pattern in 04a (`map_reads.sh`). Fix: "save as `$SCRATCH/rnaseq-workshop/scripts/quant_kallisto.sh` and submit from `scripts/`" (or use an absolute path to `samples.txt`).
- `[should-fix]` 04b `head(txi$counts)` output (`04b:330-346`) is not real output: columns are mock-first although `samples.txt` (from `ls`) is IR-first, and values such as `2.3456`, `4.5678`, `22.3456` are placeholders. Replace with the actual output (05b:136-150 shows real-looking IR-first output).
- `[check]` 04b tx2gene `awk` uses `match(..., arr)` (gawk only). Fine on Negishi (RHEL awk is gawk) and on this host; would fail with mawk on a learner's laptop.
- `[should-fix]` 04b tx2gene is built from the basic GTF (`04b:294-300`, 127,934 transcripts) but the Kallisto index is built from the comprehensive `gencode.vM38.transcripts-clean.fa` (278,396 transcripts). tximport reports "transcripts missing from tx2gene: 150462" and drops them: 12.4-12.6% of `est_counts` per sample are discarded, and multi-isoform genes are undercounted (Mdm2 baseMean 580 vs 1799 with a complete map; Bax 524 vs 947). Fix: build tx2gene from the same FASTA headers (`grep ">" gencode.vM38.transcripts.fa | cut -c2- | cut -d "|" -f 1,2 | tr "|" "\t" > tx2gene.tsv`); verified: no missing transcripts. Patch: `updates-patches/03-04b-cd-tx2gene-submitdir.diff`. The instructor notes and staged `data/tx2gene.tsv` need the same change.

Execution (R track; R 4.4.0 / Bioc 3.19 container, fresh session per episode)

- `[blocker]` confirmed by execution: LFC direction inverted (existing 2026-10-01 entry). As written, all six anchors come out DOWN in `results/deseq2/DESeq2_results_joined.tsv` (Cdkn1a -7.46, Mdm2 -2.84, Bax -3.12, Bbc3 -5.77, Pmaip1 -3.31, Gadd45a -0.51) and in the 05b table; 06 then puts every p53, DNA damage, and apoptosis term in the down-regulated ORA and gives them negative NES (-1.48 to -1.98). The committed figures `fig/05_deseq/volcano-plot.png` and `volcano-k.png` already show p53 targets on the "down" side under a title that says IR vs mock. Patch tested: `updates-patches/05-05b-lfc-direction-and-rds-path.diff` (relevel to `WT_mock`, `coef = "condition_WT_IR_vs_WT_mock"`, 05b contrast IR vs mock). After the patch all six anchors are UP in 05, five of six in 05b; 06 puts p53/DNA damage/apoptosis terms in the up-regulated ORA with positive NES (+1.47 to +1.99), KEGG p53 signaling and HALLMARK_P53_PATHWAY stay top hits, cell-cycle checkpoint terms are positive. Printed outputs and figures in 05, 05b, 06 must be regenerated.
- `[blocker]` confirmed by execution: `05:904` `saveRDS(dds, "results/deseq2_results/dds_featurecounts.rds")` errors "cannot open the connection" (directory does not exist). The patch above changes it to `results/deseq2/`; verified.
- `[check]` Anchor Gadd45a is weak: genome track padj 0.019 with shrunk LFC +0.51 (below the episode's log2(1.5) = 0.585 cutoff, so absent from `DESeq2_results_sig.tsv`); transcript track not significant (padj 0.146 as written, 0.076 with a complete tx2gene; LFC +0.48 to +0.59). Under the strict "missing anchor = blocker" rule this fails on the 05b track, but every other anchor is significant at padj < 1e-23 on both tracks and the direction is consistent, so this is effect size at 20 M reads, not a pipeline break. Decide: drop Gadd45a from the CLAUDE.md anchor list (or mark it "genome track only"), and soften the 05 solution (`05:923`) "Gadd45a should be upregulated".
- `[should-fix]` Printed counts and numbers in 05 and 06 do not match a fresh run with subread 2.0.6 (read counting): `dim(cts_coding)` 12058 vs printed 11330; summary table 12058 / 6324 / 1946 / 1937 vs 11330 / 5967 / 1802 / 1795; 06 `length(sig_genes)` 3883 vs 3597, `sig_entrez` 3895 vs 3608, up/down 1953/1942 vs 1810/1798. The top-gene list and PCA structure match. With `--countReadPairs` the counts halve and should reproduce the printed numbers. Regenerate after the 04a and 05 fixes.
- `[should-fix]` `05:561` prose says DE uses "the full count matrix (`cts`) without filtering by biotype"; the code (`05:565`) uses `cts_coding` (protein-coding, filtered). Align prose with code.
- `[should-fix]` Warnings in the current R stack, verbatim:
  - 06 block at `06:621`: "The `category` argument of `msigdbr()` is deprecated as of msigdbr 10.0.0. ℹ Please use the `collection` argument instead." (msigdbr 26.1.1). Tested fix: `collection = "H"` plus `ncbi_gene` instead of `entrez_gene` gives identical gene sets (50 sets, 7379 rows) and no warning. Patch: `updates-patches/06-msigdbr-collection.diff`. Requires msigdbr >= 10 on OOD; check the installed version first.
  - 06 `06:621`: "cannot rename file '/tmp/Rtmp.../file....zip' to '~/.cache/R/msigdbr/msigdb.2026.1.zip', reason 'Invalid cross-device link'". msigdbr >= 10 downloads MSigDB on first use and caches it; when `TMPDIR` and home are on different filesystems (as on many clusters) the cache write fails, and every session re-downloads. Hallmark still loaded. Needs compute-node internet; on OOD set `TMPDIR` under home or scratch, or pre-seed the cache.
  - 05 `05:407` `meanSdPlot()`: "`aes_string()` was deprecated in ggplot2 3.0.0" (from vsn). 06 `dotplot()`/`emapplot()`/`gseaplot2()`: "`aes_string()`/`aes_()` was deprecated in ggplot2 3.0.0", "Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0" (from enrichplot 1.24.4). Package-side; harmless, but learners will see them. Mention in a callout or move OOD to a newer enrichplot/vsn.
  - 06 `barplot(ego_bp, ...)`, `barplot(ekegg, ...)` (`06:423`, `06:579`, `06:1100`): "Arguments in `...` must be used. ✖ Problematic argument: • by = x". enrichplot 1.24 with ggplot2 4.x; harmless.
  - 06 `gseGO()` (`06:970`): "There are ties in the preranked stats (0.48% of the list)" and "For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero". Expected; explain in a callout or set `eps = 0`.
  - 06 `bitr()` (`06:225`): "0.32% of input gene IDs are fail to map..."; and "Mapped 12081 of 12058 genes (100.2%)" (duplicate mappings; existing `[nice]` entry about `unique()`).
  - All episodes: "package 'X' was built under R version 4.4.1" for most Bioconductor packages (Bioc 3.19 binaries built on R 4.4.1, run on R 4.4.0). If OOD's `4.4.0-bioconductor` library was installed the same way learners will see dozens of these on `library()`; `[check]` on OOD.
- `[nice]` 06 GSEA is not seeded: `gseGO()` returned 103 significant terms in one run and 71 in another on the same input, so learners' tables will differ from the printed output and from each other. Add `set.seed(42)` before `gseGO()` (and `seed = TRUE`).
- `[nice]` 06 captions say "for upregulated genes" (`06:409`, `06:427`) but the plotted object `ego_bp` uses all significant genes; `fig/06-enrich/go-bp-up-bar-simple.png` is shown (`06:452`) with no code that produces it. Fix captions, add `barplot(ego_bp_simple, showCategory = 15)` or drop the figure.
- `[check]` biomaRt could not be tested: on 2026-10-01 to 02 Ensembl BioMart was unavailable ("Your query has been redirected to https://status.ensembl.org", mirrors returned 502; `www.ensembl.org` redirects to `jun2026.archive.ensembl.org` which answers 403 to scripts). Both biomaRt 2.60.1 (Bioc 3.19) and 2.66.0 (Bioc 3.22) failed. The staged `data/mart.tsv` covers all 78,334 vM38 gene IDs (no gaps, no duplicates), so episode 05 is unaffected as written. Re-test the spoiler and the fallback before the delivery; do not demo biomaRt live without a fallback.
- `[should-fix]` 05 biomaRt spoiler (existing entry) re-confirmed: verbatim block fails `parse()` ("unexpected end of input" at the missing `)`). With the minimal fixes (closing paren, `annot %>%` instead of `mart %>%`, `counts` defined) it runs up to the Ensembl connection.
- 05, 05b, 06 run end to end in fresh sessions with no other errors. 06 runs unchanged from the 05 output and, with block `06:161` pointed at `results/deseq2_kallisto/DESeq2_kallisto_results.tsv`, from the 05b output (both as written and patched); this confirms the existing `[should-fix]` that 06 needs a line for 05b learners. KEGG REST (enrichKEGG) and MSigDB download worked from this machine.
- Resource use (OOD sizing check): peak RSS 05 1.5 GB, 04b tximport 3.4 GB, 05b 1.7 GB, 06 2.7 GB. All fit the 4-core (~8 GB) OOD request in setup.md. Wall time: 05 2.8 min, 05b 0.4 min, 06 6.3 min.
- `[check]` resolved: GENCODE vM38 links in episode 02 resolve (all three downloaded on 2026-10-01).

### Fixes on branch fix/oct6, 2026-10-02 (prompt-05 Part A; uncommitted, awaiting review)

Environment facts used (from `~/rnaseq_env_versions.txt`, Negishi): biocontainers modules fastqc 0.12.1, multiqc 1.23, fastp 0.23.2, star 2.7.11b, subread 2.0.1, kallisto 0.48.0, salmon 1.10.1, sra-tools 2.11.0-pl5262, seqtk 1.4. OOD R 4.4.0 on Rocky 8.10, Bioconductor 3.20 (DESeq2 1.46.0, tximport 1.34.0, clusterProfiler 4.14.6, enrichplot 1.26.6, biomaRt 2.62.1), msigdbr 26.1.1. OOD image `/depot/itap/aseethar/images/rstudio_bioc_ood_rocky8_r4.4.0_s2025.05.0-496_tex.sif`, host library `/apps/biocontainers/extras/r-ood/r4.4.0_s2025.05.0-496` bound at `/opt/R/host-site-library`. Line numbers below are in the edited files.

Placeholder convention: every output block, solution number, and timing claim that must come from the Negishi run is marked `<!-- NEGISHI:<id> -->`. A marker on its own line directly above a fence means "replace this block with the console output of the code block above it"; an inline marker follows a `TBD` in prose. Invented or known-wrong output was replaced by `TBD`; output that may be right but came from an unknown run keeps its text and the marker. 79 markers in episodes 02-06 plus 1 in `learners/setup.md`; `negishi-run/placeholders.tsv` lists them with their source. `grep -rn "NEGISHI:" episodes/ learners/` must return nothing after Part C.

Blockers
- `[blocker, fixed]` LFC direction. `05:106`, `05b:162` relevel `condition` to `WT_mock`; `05b:444` contrast now `c("condition", "WT_IR", "WT_mock")` (05 already used it); `lfcShrink(coef = "condition_WT_IR_vs_WT_mock")` at `05:719` and in 05b. New callout "Name your contrast" `05:658` explains explicit contrasts vs level order; 05b points to it. Verified locally (below).
- `[blocker, fixed]` `saveRDS` path: `05:930` saves to `results/deseq2/`; `05:91` and `05b:107` create the output directory in R with `dir.create(..., recursive = TRUE)`. The bare `mkdir -p results/deseq2` bash blocks (relative to an unknown working directory) are removed from both prereq boxes.

Should-fix
- `[should-fix, fixed]` 04a featureCounts counting mode. Negishi subread is 2.0.1, where `-p` alone counts fragments, so `--countReadPairs` is NOT added (it would error on 2.0.1). New callout "Counting reads or read pairs?" `04a:500` explains read vs fragment counting and the 2.0.2 behavior change. Prose `04a` Step 5 now says read pairs (fragments). The printed 05/06 numbers (11330 genes) were made with fragment counts, consistent with 2.0.1; local re-run with pair counts reproduced `dim(cts_coding)` = 11330.
- `[should-fix, fixed]` 04b tx2gene built from the transcript FASTA headers (`04b:300`), covering every indexed transcript; removes the gawk-only `match()` (`[check]` resolved). Prose explains why. Staged `data/tx2gene.tsv` and instructor notes still need the same change (regen_results.sh rebuilds it).
- `[should-fix, fixed]` 04b bootstraps: `-b 0` in the example (`04b:132`) and the array (`04b:216`); example uses `-t 4` to match its interactive job. Bootstrap callout rewritten (`04b:160`): when they matter (sleuth, fishpond) and when not (tximport + DESeq2). `abundance.h5` description and the h5-vs-tsv callout (`04b:368`) updated; notes the `rhdf5` dependency. 05b intro point 3 (`05b:48`) updated.
- `[should-fix, fixed]` 04b runtime claims (`04b:46`, `:101`, `:225` in the old file) removed and replaced by TBD markers.
- `[should-fix, fixed]` 04b fabricated `head(txi$counts)` output removed; command and output split into separate fences (`04b:339`). New `head -3 tx2gene.tsv` output block for the tx2gene format.
- `[should-fix, fixed]` 04b and 04a: where to save and submit array scripts. `04a:343` (map_reads.sh), `04a` count_features.sh, `04b:182` (quant_kallisto.sh); every `sbatch` block now starts with `cd $SCRATCH/rnaseq-workshop/scripts` or follows one.
- `[should-fix, fixed]` 04b heavy steps on a login node: `sinteractive` added before `kallisto index` (`04b:90`) and before the `r-rnaseq` R session (`04b:312`). "Step 0" reference fixed (`04b:154`).
- `[should-fix, fixed]` SLURM headers sized from measured runtimes plus margin, `--account=rcac-rnaseq --qos=standby --partition=cpu` unchanged: STAR index 20 CPU, `--mem=48G` (peak 32.7 GB), `--time=2:00:00` (41 min measured) at `04a:201`; STAR mapping `--mem=40G` (29.2 GB), `--time=1:00:00` (max 11.7 min/sample; was 8:00:00, which exceeds the 4 h standby cap) at `04a:350`; featureCounts `--mem=16G`, `--time=0:30:00` (2.9 min) at `04a:522`; kallisto quant `--mem=16G`, 1 h (`-b 0`). Salmon `sinteractive` raised to 2 h (`04a:64`; 28 min measured at 4 threads, contended). Staged copies: `negishi-run/staged-scripts/` (generated by build_kit.py).
- `[should-fix, fixed]` 04a account and paths: `sinteractive -A rcac-rnaseq` (`04a:64`), `$SCRATCH/rnaseq-workshop` instead of `rcac_rnaseq` (`04a:194`, star_index location).
- `[should-fix, fixed]` 04a strand prose: wrong 36.4M/35.2M replaced by TBD markers for ISF/ISR and bias; ISF/ISR definitions corrected.
- `[should-fix, fixed]` 04a mapping and counting solutions rewritten around the real metrics with TBD markers (`04a:463`, `04a:601`); assignment rate explained (denominator includes unmapped and multimapping pairs).
- `[should-fix, fixed]` 03 QC solution rewritten around the real FastQC findings with TBD markers (`03:298`).
- `[should-fix, fixed]` 05 biomaRt spoiler: `cts` not `counts`, `mart`/`annot` built in the right order, closing parenthesis, prose says biomaRt not org.Mm.eg.db, stray "0" removed (`05:134-190`). It now writes `data/mart_biomart.tsv`/`annot_biomart.tsv` so it cannot overwrite the staged files. New callout "If BioMart is unavailable" (`05:196`) with the `mart.tsv` fallback and the Ensembl-release caveat. The spoiler's `mart` has the same columns as `data/mart.tsv`, so the rest of 05 runs on either.
- `[should-fix, fixed]` 05 OOD instructions (Scholar, wrong form fields) replaced by a pointer to setup.md with the Negishi form values (`05:71`); Scholar screenshot `fig/05_deseq/open-on-demand.png` no longer referenced (orphan, see Part C).
- `[should-fix, fixed]` Unused packages removed: `EnsDb.Mmusculus.v79`, `ensembldb`, `ComplexHeatmap` (05), `tximportData` (05b). setup.md install list still lists them (not changed here; drop in the setup.md pass).
- `[should-fix, fixed]` Regression found and fixed while re-running 05: without `ensembldb` loaded first, `S4Vectors::rename` (attached by DESeq2) masks `dplyr::rename`, and the library-size block fails with "arguments in '...' must be character and not NA". Now `dplyr::rename` (`05:335`) and `dplyr::select` (`05:258`, defensive against `AnnotationDbi::select` when 06 runs in the same session).
- `[should-fix, fixed]` 05 prose said DE uses the unfiltered `cts`; now matches the code (`cts_coding`, `05:566`).
- `[should-fix, fixed]` 05b uses the shrunken table throughout: `res_df` from `res_shrunk` with `ensembl_gene_id_version` (`05b:543`), annotation joined to it, and `res_annot` saved (`05b:697`, `05b:707`). The file has the columns 06 needs. Annotation source sentence corrected (`05b:564`); dispersion plot caption fixed (`05b:421`).
- `[should-fix, fixed]` 06 msigdbr: `collection = "H"` and `ncbi_gene` (`06:623`, `06:631`); default `db_species` kept (human MSigDB mapped to mouse orthologs). Callout on the download, the cache warning, and the API change (`06:621`).
- `[should-fix, fixed]` 06 runs from either track: commented 05b line (`06:162`).
- `[should-fix, fixed]` Warning callouts: vsn/enrichplot deprecation warnings (`05:413`, `06:566`), gseGO ties and `eps` (`06:936`), msigdbr cache (`06:621`).
- `[should-fix, fixed]` 02 subsample spoiler (patch applied, `02:378`).
- `[should-fix, new, fixed]` 02 reference download block would hit `gunzip` "already exists" for learners who copied the staged data, which already has the uncompressed files; callout added (`02:252`).
- `[should-fix, new, fixed]` 02 said "100+ million reads" and "80-100+ million"; SRA run info gives 75,410,310 to 136,766,553 spots (read pairs) for the 8 runs. Both sentences now say 75 to 137 million read pairs (`02:359`, `02:371`); "20 million reads" corrected to read pairs.

Nice
- `[nice, fixed]` Output fenced as code: `02:136`, `02:412`, `03:73`, `04a` star_index listing, `04b` kallisto output tree, `04b:339`, `06` `head(gene_list)` now `text` fences separate from code. 02 tree listing order corrected (tree sorts `data, results, scripts`).
- `[nice, fixed]` 03 MultiQC stray `cd rnaseq-workshop` removed (`03:242`).
- `[nice, fixed]` 02 "four subdirectories" (`02:119`); 04a keypoint on strandedness (`04a:623`).
- `[nice, fixed]` 06 `bitr()` duplicates: `unique()` on all Entrez vectors (`06:264` and the up/down lists); mapped-genes message counts unique Ensembl IDs (`06:246`).
- `[nice, fixed]` 06 GSEA seeded: `set.seed(42)`, `seed = TRUE`, `eps = 0` (`06:912-920`).
- `[nice, fixed]` 06 captions say "all significant genes" (`06:398`, `06:416`); `barplot(ego_bp_simple)` added so `go-bp-up-bar-simple.png` has code (`06:443`).
- `[should-fix, fixed]` `learners/setup.md` package list: removed `EnsDb.Mmusculus.v79`, `ensembldb`, `ComplexHeatmap`, `tximportData`; added `rhdf5` (04b h5 import) and `hexbin` (`meanSdPlot()`); msigdbr version note now describes the `collection`/`ncbi_gene` code the episode uses; the "load-tested with R 4.5.0" claim replaced by a NEGISHI marker for the OOD package test.
- `[nice, fixed]` 301 redirects: `02:49` and 9 links in `learners/reference.md` point to the final URLs.

Local code-path check after the fixes (R 4.4.0 / Bioc 3.19 container, fresh session per episode; inputs: fragment counts from subread 2.0.6 `--countReadPairs`, complete-tx2gene `txi.rds`; harness `local-test/harness_fix/`). Code paths only; no local numbers went into the episodes.
- 05: no errors (after the `dplyr::rename` fix), 18 s, 1.1 GB. `resultsNames` = `condition_WT_IR_vs_WT_mock`; header "condition WT_IR vs WT_mock".
- 05b: no errors, 24 s, 1.7 GB. Output table has `ensembl_gene_id_version`, `pvalue`, `padj`, annotation.
- 06 genome track: no errors, 158 s, 3.2 GB; transcript track (05b line swapped in): no errors, 94 s, 2.9 GB. Both: KEGG p53 signaling and HALLMARK_P53_PATHWAY significant; "signal transduction by p53 class mediator" in the up-regulated ORA; GSEA DNA damage/p53/intrinsic apoptosis terms with positive NES (+1.78 to +1.99); no `eps` message; ties message remains (0.5 percent).
- Anchors, both tracks: all six have positive LFC. Cdkn1a, Mdm2, Bax, Bbc3, Pmaip1 padj < 1e-29. Gadd45a positive but not significant on either track (padj 0.068 genome, 0.076 transcript). Gadd45a decision deferred to Part C with the Negishi values.
- biomaRt spoiler: Ensembl still unreachable ("Unable to contact any Ensembl mirror"); the spoiler errors and the rest of 05 runs on `data/mart.tsv`. The successful-query path is untested; the Negishi run tests it.

Still open after Part A (not fixed here)
- `[should-fix]` Figures: all of `fig/05_deseq/*`, `fig/06-enrich/*`, `fig/02_qc/*` regenerate in Part C. `03:154` caption "Per tile sequence quality" shows `fastqc_per_base_sequence_content.png` (wrong image); fix when regenerating.
- `[nice]` 13 orphan figures plus `fig/05_deseq/open-on-demand.png`; delete or reuse in Part C.
- `[check]` `r-rnaseq` module (04b tximport) is not in `~/rnaseq_env_versions.txt`; preflight checks it loads and has tximport, readr, rhdf5.
- `[check]` `$SCRATCH` vs `$RCAC_SCRATCH`; preflight checks both are set and equal.
- `[check]` Staged `README.md` in the workshop data still claims ">60% assignment" and old tool versions; regen_results.sh does not rewrite it.
- `[nice]` 05b discussion table says Kallisto does "Sequence + GC bias" correction; Kallisto has only optional sequence-bias correction (`--bias`), GC bias is Salmon. Not changed.
- Static-review items unrelated to code (config.yaml `source:`, index.md, README.md, timing table, ConfoundingExplorer) unchanged.

### Negishi run kit, 2026-10-02 (prompt-05 Part B; `negishi-run/`, uncommitted)

Built to run every episode on Negishi as a new learner, fill the NEGISHI placeholders, and rebuild the completed-results copy. Usage and recovery in `negishi-run/README.md`; deviations from verbatim learner execution in `negishi-run/ADAPTATIONS.md`.

- `build_kit.py` extracts the episode blocks (same fence parsing as `local-test/harness/extract.py`), locates each by a content anchor (fails if an anchor matches 0 or 2+ blocks), and writes `generated/` (blocks, sessions, job scripts, R plans, step table, figure map, placeholder sources, episode checksums) and `staged-scripts/`. Deterministic: rebuilding gives byte-identical output (266 files); `--check` compares a fresh build with the files on disk and detects edited episodes.
- Every placeholder id has a source: 38 are the console output of the code block above the marker (the 17 in 06 from both tracks), 6 are files or listings written by the run, 37 are metrics computed by `py/summarize.py` (timings from `sacct`, FastQC module status, STAR/featureCounts/kallisto/Salmon outputs, anchors, PCA, track comparison). `build_kit.py` fails if a marker has no source.
- Safety: for `aseethar`, `$RCAC_SCRATCH/rnaseq-workshop` is the staged source. The kit points `SCRATCH`/`RCAC_SCRATCH` at a fresh test root and redirects the two literal `/scratch/negishi/$USER` paths (04a clean block, R `work_dir`) by setting `USER` to the test root's relative path; `lib.sh` refuses to run if the learner directory resolves inside the staged or results directories. The same `Sys.setenv(USER = ...)` step is in `ood_manual_check.md` so the interactive check cannot write into the staged data.
- `[should-fix, confirmed]` The January staged `scripts/map_reads.sh` and the other three staged scripts are rejected by SLURM validation (`--account=workshop`); the local shim also rejects an 8 h standby request. `negishi-run/staged-scripts/` has the episode versions; README lists the six files to replace.
- `[nice]` Ensembl BioMart answers HEAD with 405, so HTTP HEAD probes report it down when it is up; preflight and the smoke test use GET. A 200 from the registry URL did not mean `biomaRt::useMart()` worked (2026-10-02: registry 200, `useMart` "Unable to contact any Ensembl mirror").
- Verification done off-cluster: shellcheck 0.10.0 clean on all kit scripts and generated job scripts (learner code in the generated sessions and staged scripts keeps the teaching-code style notes SC2086/SC2164 already accepted above); `submit_all.sh --dry-run` through the local sbatch shim (now with `--test-only` checks of account, partition, QoS, time vs the 4 h standby cap): 18 steps, correct dependency graph, 0 failures; `--track`, `--from`, and bad `--from` behave; local run of steps 05, 05-biomart, 05b, 06-genome, 06-transcript, metrics in the kit harness (all PASS, 05-biomart WARN as designed while Ensembl is down), then `collect.sh`: 64 of 81 placeholders filled from local data (the 17 missing need the Negishi jobs and `sacct`), all 33 committed figure names regenerated.

## Runtime record

Append per-step wall time and MaxRSS from testing, plus tool versions and `sessionInfo()`. Compare totals against front matter teaching + exercises minutes and against the times claimed in the text (e.g., STAR index ~45 min).

#### Local execution test, 2026-10-01 (full scale, 22-thread laptop; full table in local-test/REPORT.md)

- 02 reference download + prep 6.5 min. 03 FastQC 16 files 25 min (contended), MultiQC 0.6 min.
- 04a salmon index + strand check 28 min (4 threads, contended); STAR index 41 min, 32.7 GB (claim 30-45 min: OK); STAR mapping 4.6-11.7 min per sample, 29.2 GB; featureCounts 2.9 min.
- 04b kallisto index 22 min (claim 2-3 min); quant `-b 100 -t 16` 20-23 min per sample uncontended (claim 2-3 min; array limit 1 h); tximport 0.4 min, 3.4 GB.
- 05 2.8 min / 1.5 GB; 05b 0.4 min / 1.7 GB; 06 6.3 min / 2.7 GB (all within the 4-core OOD request).
- Tool and package versions: local-test/versions.txt (quay.io biocontainers; R 4.4.0 / Bioc 3.19). Negishi module versions not verified.

## Revision opportunities

Non-blocking improvements worth doing before the next delivery cycle.

- Stage workshop data from Depot instead of Negishi scratch (2026-10-01). The scRNA-seq workshop stages from `/depot/workshop/data/scrna_workshop/`. Moving to e.g. `/depot/workshop/data/rnaseq_workshop/` removes the re-stage-before-every-delivery step. Learner-side paths (`$SCRATCH/rnaseq-workshop`) do not change; only the rsync sources do: `learners/setup.md` lines with `/scratch/negishi/aseethar/rnaseq-workshop` (copy command and both results lines), plus CLAUDE.md, `new-prompts/workshop-prep-todo.md`, and this file. No episode references the staged path.

## Done

Resolved items with date and commit.
