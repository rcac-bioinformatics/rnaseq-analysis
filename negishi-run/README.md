# Negishi run kit

Runs every episode of the lesson on Negishi as a new learner would, in the real environment (biocontainers modules, the `rcac-rnaseq` account on `standby`, and the Open OnDemand R 4.4.0 image), and brings back one tarball with every output, figure, timing, and number the lesson needs. It also rebuilds the precomputed "completed results" copy.

All learner commands are generated from the episodes by `build_kit.py`; see `docs/ADAPTATIONS.md` for the few places where unattended running differs from a learner. The staged data is only read: each run works in a fresh copy under `$KIT_SCRATCH_BASE` (default `/scratch/negishi/$USER/rnaseq-kit`).

## Command order

Run everything as `aseethar` from a Negishi login node. Use `tmux` or `screen` for step 2.

| Step | Command | Where | Expected wall time |
|---|---|---|---|
| 0 | send the kit (below) | laptop | 1 min |
| 0b | `bash ~/rnaseq-analysis/negishi-run/restage.sh` (only when the staged copy changes) | login node | 10 to 20 min |
| 1 | `bash ~/rnaseq-analysis/negishi-run/00_preflight.sh` | login node | 3 to 5 min |
| 2 | `bash ~/rnaseq-analysis/negishi-run/01_learner_setup.sh 2026-10-03` | login node | 10 to 20 min (copies ~20 GB) |
| 3 | `bash ~/rnaseq-analysis/negishi-run/submit_all.sh --run 2026-10-03 --dry-run` | login node | 1 min |
| 4 | `bash ~/rnaseq-analysis/negishi-run/submit_all.sh --run 2026-10-03` | login node | 1 min to submit; jobs run 2.5 to 4 h with standby waits |
| 5 | `docs/ood_manual_check.md`, once step `04a-post` is done | browser | 20 min |
| 6 | `bash ~/rnaseq-analysis/negishi-run/collect.sh --run 2026-10-03` | login node | 1 to 3 min |
| 7 | bring the tarball back (below) | laptop | 1 min |

Step 1 must show no FAIL before you continue. Fix FAIL rows first (most likely: staged data purged or unreadable, account membership, a changed module default). WARN rows are listed in the summary for review; the run can proceed.

Expected run time per stage after submission (compute only, both tracks in parallel; standby queue waits come on top and are what the run measures):

| Stage | Steps | Compute |
|---|---|---|
| Shared | 02 (fresh reference download, subsample test), 03 (FastQC 16 files, MultiQC) | 20 to 45 min |
| Genome track | 04a Salmon, STAR index, STAR mapping (8-task array), featureCounts, then 05, 06 | about 2.2 h on the critical path (index about 1 h) |
| Transcript track | 04b kallisto index (batch, 48 cores), quant (8-task array), tximport, then 05b, 06 | about 1.3 h |
| Kit metrics | after 05, 05b | 5 min |

Watch the queue with `squeue -u $USER -o '%.10i %.22j %.8T %.10M %R'`. Kit jobs are named `kit-<step>`; the learner scripts keep their own names (`star_index`, `read_mapping`, `featurecounts`, `kallisto_index`, `kallisto_quant`).

## Sending the kit and returning the results

From the laptop, send the episodes, the setup page, and the kit (the kit reads the episodes to check it was built from the same text, and for the front-matter times):

```bash
cd ~/svn/rnaseq-analysis
python3 negishi-run/build_kit.py --check          # must say "match a fresh build"
rsync -av --relative ./episodes ./learners ./negishi-run aseethar@negishi.rcac.purdue.edu:rnaseq-analysis/
```

After `collect.sh`, bring the tarball back (the command is printed by `collect.sh`):

```bash
rsync -avP aseethar@negishi.rcac.purdue.edu:/scratch/negishi/aseethar/rnaseq-kit/negishi-run-2026-10-03.tar.gz ~/svn/rnaseq-analysis/
```

Start with `SUMMARY.md` inside it: per-step and per-episode status with the first error of every failed step, the six anchor genes on both tracks, runtime against each episode's front matter, every placeholder value, and the R warnings learners will see.

## Options and recovery

`submit_all.sh --run ID [--track genome|transcript|both] [--qos standby|normal] [--from STEP] [--prebuilt-index] [--force] [--dry-run]`

- Reruns are safe. A step that already finished in this run (its marker in `runs/ID/markers/` says PASS or WARN, or its outputs exist for the learner scripts) is skipped; `--force` reruns everything.
- A failed step cancels its dependents (`--kill-on-invalid-dep=yes`). Read `runs/ID/logs/<step>.<jobid>.out` and `runs/ID/logs/session-<step>.log` (shell steps) or `runs/ID/logs/R-<step>.log` (R steps). Learner-script logs are in `runs/ID/scratch/rnaseq-workshop/scripts/cluster-*.out`. After fixing the cause, resubmit from that step: `submit_all.sh --run ID --from 04a-map`. Steps before it are not resubmitted; dependencies on them are dropped if they are done, and the command stops if one is not.
- Step names, in order (generated/steps.tsv): `02 03 04a-salmon 04a-index 04a-map 04a-mapqc 04a-count 04a-post 04b-index 04b-quant 04b-tximport 04b-qc 05 05-biomart 05b 06-genome 06-transcript metrics`.
- `--qos normal` uses the normal QoS (preflight shows whether `rcac-rnaseq` has it); default `standby` measures what learners experience.
- `--prebuilt-index` links the STAR index from the completed results copy, as the episode 04a callout tells learners to, instead of building it.
- If the episodes change after the kit was generated, preflight warns. Rebuild on Negishi with `python3 negishi-run/build_kit.py` (Python 3.6 is enough) and start a new run ID.

## Regenerating the completed results

```bash
bash ~/rnaseq-analysis/negishi-run/regen_results.sh --date 2026-10-03          # setup + both tracks, run ID regen-2026-10-03
bash ~/rnaseq-analysis/negishi-run/collect.sh --run regen-2026-10-03           # check its SUMMARY.md
bash ~/rnaseq-analysis/negishi-run/regen_results.sh --finalize --date 2026-10-03
```

`--finalize` refuses unless every step passed, checks that everything `instructor-notes/README.md` lists is present, copies the finished learner directory to a new `/scratch/negishi/aseethar/rnaseq-workshop_results.2026-10-03` (plus the transcript-track enrichment as `results/enrichment_kallisto/`), and writes the file-list diff against the current copy to `runs/regen-2026-10-03/records/regen-filelist-diff.txt`. It does not change the current `rnaseq-workshop_results` and does not change permissions.

### Swapping in the regenerated results (manual, after you have checked them)

```bash
cd /scratch/negishi/aseethar
chmod -R o+rX rnaseq-workshop_results.2026-10-03          # learners must be able to read it
mv rnaseq-workshop_results rnaseq-workshop_results.retired-2026-10-03
mv rnaseq-workshop_results.2026-10-03 rnaseq-workshop_results
ls -l rnaseq-workshop_results/data/star_index/SA rnaseq-workshop_results/results/counts/gene_counts_clean.txt
```

Rollback:

```bash
cd /scratch/negishi/aseethar
mv rnaseq-workshop_results rnaseq-workshop_results.2026-10-03
mv rnaseq-workshop_results.retired-2026-10-03 rnaseq-workshop_results
```

## Staging the learner copy on Depot

The staged copy learners rsync in setup lives on Depot at `/depot/DEPOT_PATH/rnaseq-workshop` (`STAGED` in `config.sh`). Replace `DEPOT_PATH` in `config.sh`, `learners/setup.md`, and `.claude/CLAUDE.md` first; the kit refuses to run until you do.

It must hold exactly what `staged_manifest.tsv` lists (31 files): the 16 FASTQ files, the four GENCODE vM38 references, `mart.tsv`, `annot.tsv`, `SRR_Acc_List.txt` (the 8 workshop runs; GSE71176 has 24, see Episode 02), `tx2gene.tsv` (built from the transcript FASTA, as Episode 04b does), and the seven files in `scripts/`. Nothing else: no indexes, `results/`, job logs, `README.md`, or `.ipynb_checkpoints/`, because learners would receive them and the kit would treat finished outputs as done.

`restage.sh` builds that copy from the old scratch staging (`STAGED_OLD`) and checks it:

```bash
bash ~/rnaseq-analysis/negishi-run/restage.sh              # asks before copying; --yes to skip
```

It copies only the manifest's data files (so the extras stay behind), regenerates `SRR_Acc_List.txt` and `tx2gene.tsv`, installs `staged-scripts/`, runs `chmod -R g+rX` on the Depot copy, and then runs `check_staged.sh --perm group --tree` (full md5 check, a few minutes). It never deletes anything; rerunning updates files in place. Send back the `data/` listing it prints so Episode 02 can be refilled.

The Depot copy is group-readable, so every learner must be a member of the Depot group (or use `STAGED_PERM=other` with a world-readable copy). Preflight checks the permissions you choose.

`check_staged.sh [DIR] [--quick] [--tree] [--perm group|other]` checks any copy on its own: presence, size, md5, permissions, and extra files. The manifest is generated on the laptop by `make_staged_manifest.py` from verified copies.

## What is not stored in git

- `generated/`: rebuild with `python3 negishi-run/build_kit.py` after cloning or pulling (Python 3.6 is enough). The rsync in "Sending the kit" copies it from the laptop.
- `staged-data/`: the regenerated `SRR_Acc_List.txt` and `tx2gene.tsv`. `make_staged_manifest.py` rebuilds them on the laptop and `restage.sh` regenerates them on Negishi; the manifest (tracked) holds their md5s.
- Result tarballs (`*.tar.gz`).

## What could not be verified before the first Negishi run

The kit was built and tested off-cluster: `build_kit.py` (deterministic, `--check`), shellcheck 0.10.0 on every kit script and generated job, a full `submit_all.sh --dry-run` through a local `sbatch` shim that enforces the account, partition, QoS, and the 4-hour standby limit, and a local run of the R steps (05, 05-biomart, 05b, 06 on both tracks, metrics), step 02, `summarize.py`, and `collect.sh` on real outputs. These parts depend on Negishi and are covered as follows:

| Not verifiable off-cluster | What covers it |
|---|---|
| Module names and default versions; `r-rnaseq` existing and having tximport, readr, rhdf5 | preflight loads each module, checks the default against `config.sh`, and runs R in `r-rnaseq` |
| `rcac-rnaseq` membership, standby/normal QoS, partition limits, memory requests | preflight checks `sacctmgr` and runs `sbatch --test-only` on every learner script and kit job |
| Staged data present, complete, world-readable, expected sizes | preflight checks `staged_manifest.tsv` (sizes, md5 of the annotation tables), path permissions, and unreadable files |
| Prebuilt indexes and the completed results copy | preflight reports each; the 04a callout path is checked specifically |
| OOD image, host library binding, R and Bioconductor versions, package set | preflight runs R in the image with `--cleanenv`, detects whether `R_LIBS_SITE` is needed, and loads every package the R episodes use |
| Compute-node internet (KEGG, BioMart, MSigDB via zenodo, ftp.ebi.ac.uk, NCBI SRA) | preflight runs a 5-minute standby job that requests each URL |
| Whether `module` is available in a non-login shell inside jobs | preflight checks; sessions fall back to Lmod's init file |
| `$SCRATCH` and `$RCAC_SCRATCH` equal `/scratch/negishi/$USER` (assumed by the USER redirection) | preflight checks |
| The OOD form, launch time, and the interactive experience | `docs/ood_manual_check.md` |
| Real runtimes, memory, and queue waits | the run itself (`sacct` in SUMMARY.md) |
