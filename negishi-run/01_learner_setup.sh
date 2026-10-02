#!/bin/bash
# 01_learner_setup.sh: start a kit run by doing the learners/setup.md data copy exactly as
# written, into a fresh test directory, then the setup.md verification commands.
# Login node (this is where learners run it). Takes about 10 to 20 minutes for ~20 GB;
# run it inside tmux or screen.
# Usage: bash 01_learner_setup.sh [RUN_ID]     (default RUN_ID: YYYY-MM-DD-HHMM)
# Creates $KIT_SCRATCH_BASE/runs/<RUN_ID>/ with scratch/rnaseq-workshop as the learner copy.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=${1:-$(date +%Y-%m-%d-%H%M)}
RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID"
if [ -e "$RUN/markers/01.done" ] && kit_is_done 01 2>/dev/null; then
  kit_log "run $RUN_ID already set up ($(cat "$RUN/markers/01.done")); nothing to do"
  exit 0
fi
[ -d "$STAGED" ] || kit_die "staged data $STAGED not found; run 00_preflight.sh"
kit_init_run "$RUN"
export RUN
kit_load_run   # exports SCRATCH and RCAC_SCRATCH = $RUN/scratch (the test learner's scratch)
kit_marker 01 RUNNING
kit_log "run $RUN_ID: learner scratch is $RCAC_SCRATCH (staged source $STAGED is only read)"

# learner session: the three setup.md blocks, verbatim
S="$RUN/tmp/01-setup-session.sh"
{
  echo 'set -E'
  echo "trap 'echo \"### KIT-ERR rc=\$? line=\$LINENO cmd=\$BASH_COMMAND\" >&2' ERR"
  # shellcheck disable=SC2016  # written into the session file, expanded there
  echo 'cd "$RCAC_SCRATCH"'
  for f in 01-echo 02-rsync 03-verify; do
    echo "echo \"### KIT-BLOCK setup/$f begin \$(date +%s)\""
    cat "$KIT_GEN/setup/$f.sh"
    echo "echo \"### KIT-BLOCK setup/$f end \$(date +%s)\""
  done
} > "$S"
kit_run_session 01 "$S"

# checks on what the learner now has
LOG="$RUN/logs/session-01.log"
n=$(find "$W/data" -maxdepth 1 -name '*.fastq.gz' | wc -l)
[ "$n" = 16 ] || kit_die "learner copy has $n FASTQ files, expected 16"
# setup.md promises learners start without indexes, results, or job logs; a copy that
# already has them would let finished outputs stand in for the episode code
pre=$(kit_precomputed "$W")
if [ -n "$pre" ] && [ "${KIT_ALLOW_PRECOMPUTED:-0}" != 1 ]; then
  kit_marker 01 FAIL
  kit_die "the staged copy $STAGED contains precomputed outputs, so a learner (and this run) starts with finished results:
$pre
Remove them from the staged directory (see README.md), then start a new run ID. KIT_ALLOW_PRECOMPUTED=1 overrides."
fi
awk '/^### KIT-BLOCK setup\/03-verify begin/{f=1;next} /^### KIT-BLOCK setup\/03-verify end/{f=0} f' "$LOG" \
  | grep -E '[0-9.]+[KMGT]?\s+/' | tail -1 | awk '{print $1}' > "$RUN/records/listings/01-du.txt"
secs=$(awk '/KIT-BLOCK setup\/02-rsync begin/{b=$NF} /KIT-BLOCK setup\/02-rsync end/{e=$NF} END{print e-b}' "$LOG")
printf 'setup rsync\t%s\n' "$secs" > "$RUN/records/listings/01-rsync-seconds.tsv"
# listing for 02-out-data-tree (the data directory "after downloading")
( cd "$W" && { tree data 2>/dev/null || ls -1 data; } ) > "$RUN/records/listings/02-data-tree.txt"
kit_job_end 01 PASS
kit_log "setup done in ${secs}s; data size $(cat "$RUN/records/listings/01-du.txt"). Next: bash $KIT_DIR/submit_all.sh --run $RUN_ID"
