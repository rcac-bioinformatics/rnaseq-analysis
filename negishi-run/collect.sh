#!/bin/bash
# collect.sh: gather a kit run into one tarball to bring back (login node, 1 to 3 minutes).
# Usage: bash collect.sh --run RUN_ID
# Writes $KIT_SCRATCH_BASE/negishi-run-<RUN_ID>.tar.gz containing SUMMARY.md,
# placeholder_values.tsv, outputs/ (console output per NEGISHI placeholder), QC and
# quantification summaries, DE and enrichment tables, figures named as in episodes/fig/,
# sacct records, sessionInfo per R step, module versions, and all logs. Safe to rerun.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""
while [ $# -gt 0 ]; do
  case $1 in --run) RUN_ID=$2; shift 2;; *) kit_die "usage: collect.sh --run RUN_ID";; esac
done
[ -n "$RUN_ID" ] || kit_die "usage: collect.sh --run RUN_ID"
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN
kit_load_run
NAME="negishi-run-$RUN_ID"
C="$RUN/collect/$NAME"
mkdir -p "$C"/{qc,quant,analysis,figures,records,logs,checks,learner-logs}
kit_log "collecting $RUN into $C"

# ---- SLURM accounting for every job this run submitted (array tasks included)
ids=$(awk -F'\t' 'NR > 1 && $2 ~ /^[0-9]+$/ { printf "%s%s", sep, $2; sep = "," }' "$RUN/jobs.tsv" 2>/dev/null || true)
if [ -n "$ids" ]; then
  sacct -j "$ids" -P -n --format=JobID,JobName%40,Elapsed,MaxRSS,ReqMem,State,Submit,Start,End,AllocCPUS,Timelimit,ExitCode \
    > "$RUN/records/sacct.txt" 2>/dev/null || kit_log "sacct failed; runtimes will be missing"
fi

# ---- versions and preflight state
mkdir -p "$RUN/records/preflight"
latest=$(find "$KIT_SCRATCH_BASE/state" -maxdepth 1 -name 'preflight-*.txt' 2>/dev/null | sort | tail -1)
[ -n "$latest" ] && cp "$latest" "$RUN/records/preflight/"
[ -f "$KIT_SCRATCH_BASE/state/pkgcheck.txt" ] && cp "$KIT_SCRATCH_BASE/state/pkgcheck.txt" "$RUN/records/preflight/"
{
  for mv in $EXPECTED_MODULES r-rnaseq; do
    tool=${mv%%/*}
    printf '%s\t%s\n' "$tool" "$(bash -c "module load biocontainers >/dev/null 2>&1; module load $tool >/dev/null 2>&1; module -t list 2>&1 | grep -E '^$tool/' | head -1" 2>/dev/null || echo '?')"
  done
} > "$RUN/records/module_versions.tsv"
(cd "$KIT_DIR/.." && sha256sum -c "$KIT_GEN/episodes.sha256" > "$C/checks/episode_checksums.txt" 2>&1) || true

# ---- summary (also extracts per-placeholder outputs)
python3 "$KIT_DIR/py/summarize.py" "$RUN" "$KIT_GEN" "$C"

# ---- QC and quantification summaries from the learner directory
cpq() { local dst=$1; shift; mkdir -p "$dst"; local f; for f in "$@"; do [ -e "$f" ] && cp -r "$f" "$dst/"; done; return 0; }
cpq "$C/qc/fastqc" "$W"/results/qc_fastq/*_fastqc.zip
for d in qc_fastq qc_alignment qc_counts qc_kallisto; do
  cpq "$C/qc/$d" "$W/results/$d/multiqc_report.html" "$W/results/$d/multiqc_data"
done
cpq "$C/qc/multiqc_export_qc" "$RUN"/records/multiqc_export_qc/*
cpq "$C/qc/fastp" "$RUN"/records/fastp/*.json "$RUN"/records/fastp/*.html
cpq "$C/quant/star" "$W"/results/mapping/*Log.final.out
cpq "$C/quant/featurecounts" "$W/results/counts/gene_counts.txt.summary"
cpq "$C/quant/salmon" "$W/results/strand_check/lib_format_counts.json" "$W/results/strand_check/logs" "$W/results/strand_check/aux_info/meta_info.json"
for s in "$W"/results/kallisto_quant/*/; do
  [ -d "$s" ] || continue
  b=$(basename "$s"); cpq "$C/quant/kallisto/$b" "$s/run_info.json" "$s/$b.log"
done
cpq "$C/quant/listings" "$RUN"/records/listings/*

# ---- analysis outputs
cpq "$C/analysis/deseq2" "$W"/results/deseq2/*.tsv
cpq "$C/analysis/deseq2_kallisto" "$W"/results/deseq2_kallisto/*.tsv
cpq "$C/analysis/enrichment_genome" "$W"/results/enrichment/*.csv "$W"/results/enrichment/*.pdf
cpq "$C/analysis/enrichment_transcript" "$W_T"/results/enrichment/*.csv "$W_T"/results/enrichment/*.pdf
cpq "$C/analysis/biomart" "$W/data/mart_biomart.tsv"
cpq "$C/analysis/metrics" "$RUN"/records/metrics/*

# ---- figures: raw per block, and copies named as the committed episodes/fig/ files
cpq "$C/figures/raw" "$RUN"/records/figures/*
while IFS=$'\t' read -r step ep block k fig; do
  [ "$step" = step ] && continue
  src="$RUN/records/figures/$step/${ep}_${block}_$(printf '%02d' "$k").png"
  dst="$C/figures/episodes/$fig"
  [ "$step" = 06-transcript ] && dst="$C/figures/episodes_transcript_track/$fig"
  if [ -f "$src" ]; then mkdir -p "$(dirname "$dst")"; cp "$src" "$dst"
  else printf '%s\t%s\n' "$fig" "$src" >> "$C/figures/MISSING.tsv"; fi
done < "$KIT_GEN/figmap.tsv"

# ---- records and logs
cp -r "$RUN"/records/{console,events,sessioninfo,preflight} "$C/records/" 2>/dev/null || true
cp "$RUN"/records/*.tsv "$RUN"/records/*.txt "$C/records/" 2>/dev/null || true
cp "$RUN/kit.env" "$RUN"/jobs*.tsv "$C/records/" 2>/dev/null || true
cp -r "$RUN"/markers "$C/records/" 2>/dev/null || true
cp "$RUN"/logs/* "$C/logs/" 2>/dev/null || true
cp "$W"/scripts/cluster-* "$C/learner-logs/" 2>/dev/null || true
cp "$KIT_GEN/placeholders.tsv" "$KIT_GEN/steps.tsv" "$C/records/"

TAR="$KIT_SCRATCH_BASE/$NAME.tar.gz"
tar -czf "$TAR" -C "$RUN/collect" "$NAME"
kit_log "wrote $TAR ($(du -h "$TAR" | cut -f1)); summary: $C/SUMMARY.md"
kit_log "bring it back with: rsync -avP $(id -un)@negishi.rcac.purdue.edu:$TAR ."
