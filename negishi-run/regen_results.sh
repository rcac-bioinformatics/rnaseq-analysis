#!/bin/bash
# regen_results.sh: rebuild the precomputed "completed results" copy with the fixed lesson code.
# Two phases, so you can check the run before anything is copied:
#   bash regen_results.sh [--qos standby|normal] [--date YYYY-MM-DD]
#       runs 01_learner_setup.sh and submit_all.sh (both tracks) as run regen-<date>
#   bash regen_results.sh --finalize [--date YYYY-MM-DD]
#       after every step passed: copies the finished learner directory into a NEW directory
#       <dirname of STAGED_RESULTS>/rnaseq-workshop_results.<date>, leaving out the FASTQ and
#       reference files learners already have from the staged copy, and diffs its file list
#       against the current STAGED_RESULTS. The current copy is never touched; swapping
#       the new one into place is a separate manual step (README.md).
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

DATE=$(date +%F); QOS=$KIT_QOS_DEFAULT; FINAL=0
while [ $# -gt 0 ]; do
  case $1 in
    --date) DATE=$2; shift 2;;
    --qos) QOS=$2; shift 2;;
    --finalize) FINAL=1; shift;;
    *) kit_die "usage: regen_results.sh [--finalize] [--date YYYY-MM-DD] [--qos standby|normal]";;
  esac
done
RUN_ID="regen-$DATE"
RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID"
DEST="$(dirname "$STAGED_RESULTS")/$(basename "$STAGED_RESULTS").$DATE"

if [ $FINAL = 0 ]; then
  [ -e "$DEST" ] && kit_die "$DEST already exists; pick another --date"
  bash "$KIT_DIR/01_learner_setup.sh" "$RUN_ID"
  bash "$KIT_DIR/submit_all.sh" --run "$RUN_ID" --track both --qos "$QOS"
  kit_log "when all jobs finish: bash $KIT_DIR/collect.sh --run $RUN_ID, check SUMMARY.md, then bash $0 --finalize --date $DATE"
  exit 0
fi

export RUN
kit_load_run
[ -e "$DEST" ] && kit_die "$DEST already exists; refusing to overwrite"
bad=""
while IFS=$'\t' read -r step _; do
  [ "$step" = step ] && continue
  kit_step_passed "$step" || bad="$bad $step"
done < "$KIT_GEN/steps.tsv"
[ -z "$bad" ] || kit_die "not every step passed:$bad. Fix and rerun with submit_all.sh --run $RUN_ID --from <step>"

# what instructor-notes/README.md lists as precomputed backups
for p in data/star_index/SA data/salmon_index data/kallisto_index/transcripts.idx data/tx2gene.tsv \
         results/qc_fastq/multiqc_report.html results/qc_alignment/multiqc_report.html results/qc_counts/multiqc_report.html \
         results/qc_kallisto/multiqc_report.html results/counts/gene_counts_clean.txt results/kallisto_quant/txi.rds \
         results/deseq2/DESeq2_results_joined.tsv results/deseq2/dds_featurecounts.rds \
         results/deseq2_kallisto/DESeq2_kallisto_results.tsv results/enrichment/GO_BP_enrichment.csv \
         scripts/samples.txt scripts/samples.csv; do
  [ -e "$W/$p" ] || kit_die "expected $W/$p is missing"
done
[ "$(find "$W/results/mapping" -maxdepth 1 -name '*.bam' | wc -l)" = 8 ] || kit_die "expected 8 BAM files in $W/results/mapping"

kit_log "copying $W -> $DEST"
mkdir "$DEST"
# FASTQ and references are in the staged copy already (about 20 GB); learners copy single files from here
rsync -a --exclude='/data/*.fastq.gz' --exclude=/data/GRCm39.primary_assembly.genome.fa \
  --exclude=/data/gencode.vM38.primary_assembly.basic.annotation.gtf \
  --exclude=/data/gencode.vM38.transcripts.fa --exclude=/data/gencode.vM38.transcripts-clean.fa "$W/" "$DEST/"
# the transcript-track enrichment ran in its own directory; keep it next to the genome one
if [ -d "$W_T/results/enrichment" ]; then rsync -a "$W_T/results/enrichment/" "$DEST/results/enrichment_kallisto/"; fi

DIFF="$RUN/records/regen-filelist-diff.txt"
( cd "$STAGED_RESULTS" 2>/dev/null && find . -mindepth 1 | sort ) > "$RUN/records/regen-old-files.txt" || true
( cd "$DEST" && find . -mindepth 1 | sort ) > "$RUN/records/regen-new-files.txt"
{
  echo "# file list: $STAGED_RESULTS (old) vs $DEST (new)"
  echo "# only in old: $(comm -23 "$RUN/records/regen-old-files.txt" "$RUN/records/regen-new-files.txt" | wc -l)"
  comm -23 "$RUN/records/regen-old-files.txt" "$RUN/records/regen-new-files.txt" | sed 's/^/- /'
  echo "# only in new: $(comm -13 "$RUN/records/regen-old-files.txt" "$RUN/records/regen-new-files.txt" | wc -l)"
  comm -13 "$RUN/records/regen-old-files.txt" "$RUN/records/regen-new-files.txt" | sed 's/^/+ /'
} > "$DIFF"
kit_log "done: $DEST ($(du -sh "$DEST" | cut -f1)). File-list diff: $DIFF"
kit_log "not world-readable yet and not swapped in; see README.md, 'Swapping in the regenerated results'"
