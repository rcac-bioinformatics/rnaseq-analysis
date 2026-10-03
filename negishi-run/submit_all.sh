#!/bin/bash
# submit_all.sh: submit the whole episode graph for one kit run (login node, under a minute).
# Usage:
#   bash submit_all.sh --run RUN_ID [--track genome|transcript|both] [--qos standby|normal]
#                      [--from STEP] [--prebuilt-index] [--force] [--dry-run]
# Steps and dependencies come from generated/steps.tsv (built from the episodes). Both
# tracks start in parallel after setup; 06 runs once per track. A step that is already
# done in this run (marker or outputs present) is skipped unless --force; --from STEP
# skips everything before STEP in the table. Job IDs go to $RUN/jobs.tsv.
# --dry-run submits nothing: it prints each sbatch command, checks every script with
# bash -n, and asks SLURM to validate the job headers with sbatch --test-only.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""; TRACK=both; QOS=$KIT_QOS_DEFAULT; FROM=""; PREBUILT=0; FORCE=0; DRY=0
usage() { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [ $# -gt 0 ]; do
  case $1 in
    --run) RUN_ID=$2; shift 2;;
    --track) TRACK=$2; shift 2;;
    --qos) QOS=$2; shift 2;;
    --from) FROM=$2; shift 2;;
    --prebuilt-index) PREBUILT=1; shift;;
    --force) FORCE=1; shift;;
    --dry-run) DRY=1; shift;;
    -h|--help) usage 0;;
    *) echo "unknown option $1" >&2; usage 1;;
  esac
done
case $TRACK in genome|transcript|both) ;; *) kit_die "--track must be genome, transcript, or both";; esac
case $QOS in standby|normal) ;; *) kit_die "--qos must be standby or normal";; esac
if [ -n "$FROM" ] && ! awk -F'\t' -v s="$FROM" '$1 == s { f = 1 } END { exit !f }' "$KIT_GEN/steps.tsv"; then
  kit_die "--from $FROM: no such step; steps are: $(awk -F'\t' 'NR > 1 { printf "%s ", $1 }' "$KIT_GEN/steps.tsv")"
fi
[ -n "$RUN_ID" ] || kit_die "--run RUN_ID is required (the ID given to 01_learner_setup.sh)"
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN

TMPRUN=0
if [ $DRY = 1 ] && [ ! -f "$RUN/kit.env" ]; then
  kit_log "dry run: $RUN has no kit.env; using a temporary run directory (nothing is written to $RUN)"
  RUN=$(mktemp -d "${TMPDIR:-/tmp}/kit-dryrun.XXXXXX")
  KIT_ALLOW_ANY_ROOT=1 kit_init_run "$RUN"
  TMPRUN=1
fi
if [ $TMPRUN = 1 ]; then
  # shellcheck source=/dev/null
  source "$RUN/kit.env"
  export SCRATCH="$KIT_LEARNER_SCRATCH" RCAC_SCRATCH="$KIT_LEARNER_SCRATCH" W W_T
else
  kit_load_run
  kit_is_done 01 || [ $DRY = 1 ] || kit_die "setup not done for $RUN_ID; run 01_learner_setup.sh $RUN_ID first"
fi
if [ -f "$KIT_SCRATCH_BASE/state/prebuilt.sh" ]; then
  # shellcheck source=/dev/null
  source "$KIT_SCRATCH_BASE/state/prebuilt.sh"
fi

# metrics only summarizes other steps' outputs and is cheap, so it is never treated as done
step_done() { [ $FORCE = 1 ] && return 1; [ "$1" = metrics ] && return 1; kit_step_passed "$1"; }
step_track() { awk -F'\t' -v s="$1" '$1 == s { print $3 }' "$KIT_GEN/steps.tsv"; }
in_track() { local t; t=$(step_track "$1"); [ "$TRACK" = both ] || [ "$t" = common ] || [ "$t" = "$TRACK" ]; }

# refuse while an earlier submission for this run is still queued or running: two copies
# of a step would work in the same learner directory at the same time
if [ $DRY = 0 ] && [ -f "$RUN/jobs.tsv" ]; then
  prev=$(awk -F'\t' 'NR > 1 && $2 ~ /^[0-9]+$/ { printf "%s%s", s, $2; s = "," }' "$RUN/jobs.tsv")
  if [ -n "$prev" ]; then
    live=$(squeue -h -j "$prev" -o '%i %j %T' 2>/dev/null || true)
    [ -z "$live" ] || kit_die "jobs from an earlier submission of this run are still queued or running; wait for them or scancel them first:
$live"
  fi
fi

declare -A JOB
SUBMITTED=0; FAKE=9000000
JOBS_TSV="$RUN/jobs.tsv"
[ $DRY = 1 ] && JOBS_TSV="$RUN/jobs-dryrun.tsv"
[ -f "$JOBS_TSV" ] || printf 'step\tjobid\tsubmitted\tdeps\tqos\n' > "$JOBS_TSV"

sbatch_or_dry() {  # step workdir args... ; prints the job id; returns 1 if submission (or test-only) failed
  local step=$1 wd=$2; shift 2
  if [ $DRY = 1 ]; then
    local args=() a
    for a in "$@"; do case $a in --dependency=*|--parsable|--kill-on-invalid-dep=*) ;; *) args+=("$a");; esac; done
    printf '  (cd %s && sbatch %s)\n' "$wd" "$*" >&2
    local r
    r=$(cd "$wd" 2>/dev/null || cd "$RUN"; sbatch --test-only "${args[@]}" 2>&1 | tail -1) || true
    echo "$((FAKE + SUBMITTED + 1))"   # fake id; the parent shell counts submissions
    case "$r" in *"to start at"*) printf '  test-only: %s\n' "$r" >&2;;
      *) printf '  test-only FAILED: %s\n' "$r" >&2; return 1;; esac
  else
    (cd "$wd" && sbatch "$@")
  fi
}
submitted_or_fail() {  # called with the exit status of sbatch_or_dry
  if [ "$1" -ne 0 ]; then
    if [ $DRY = 1 ]; then DRYFAIL=$((DRYFAIL + 1)); else kit_die "sbatch failed for $step"; fi
  fi
}

DRYFAIL=0
started=0
[ -z "$FROM" ] && started=1
while IFS=$'\t' read -r step kind _ deps dep_type _ _ _ _ _ script desc; do
  [ "$step" = step ] && continue
  if [ $started = 0 ]; then
    if [ "$step" = "$FROM" ]; then started=1; else continue; fi
  fi
  in_track "$step" || continue
  if step_done "$step"; then kit_log "skip $step (already done)"; continue; fi

  # dependencies: only on steps submitted in this invocation; others must be done already
  dl=""
  for d in ${deps//,/ }; do
    [ "$d" = - ] && continue
    if [ -n "${JOB[$d]:-}" ]; then dl="$dl:${JOB[$d]}"
    elif step_done "$d" || [ "$dep_type" = afterany ] || ! in_track "$d"; then :
    elif [ $DRY = 1 ]; then kit_log "dry run: $step depends on $d, which is not done"
    else kit_die "$step depends on $d, which is neither done nor submitted; use --from $d"; fi
  done
  dep=(); [ -n "$dl" ] && dep=("--dependency=$dep_type$dl")
  common=(--parsable --kill-on-invalid-dep=yes "--qos=$QOS")
  [ ${#dep[@]} -gt 0 ] && common+=("${dep[@]}")
  kit_log "submit $step ($desc)${dl:+ after${dl//:/ }}"

  case $kind in
    kitjob|ood)
      bash -n "$KIT_GEN/jobs/$step.sh"
      [ -f "$KIT_GEN/sessions/$step.sh" ] && bash -n "$KIT_GEN/sessions/$step.sh"
      rc=0
      jid=$(sbatch_or_dry "$step" "$RUN" "${common[@]}" -o "$RUN/logs/$step.%j.out" -e "$RUN/logs/$step.%j.out" \
            "--export=ALL,RUN=$RUN,KIT_DIR=$KIT_DIR,KIT_GEN=$KIT_GEN" "$KIT_GEN/jobs/$step.sh") || rc=$?
      submitted_or_fail $rc
      ;;
    learner)
      if [ "$step" = 04a-index ] && [ $PREBUILT = 1 ]; then
        [ -n "${PREBUILT_STAR:-}" ] || kit_die "--prebuilt-index: no prebuilt STAR index (see 00_preflight.sh)"
        kit_log "04a-index: using the episode's prebuilt-index callout instead of building"
        if [ $DRY = 0 ]; then kit_run_session 04a-index-prebuilt "$KIT_GEN/submit-time/04a-index-prebuilt.sh"; kit_marker 04a-index PASS; fi
        continue
      fi
      # the learner saves the script as the episode says, runs any login-node commands, then submits
      bash -n "$KIT_GEN/learner-scripts/$script"
      if [ $DRY = 0 ]; then
        cp "$KIT_GEN/learner-scripts/$script" "$W/scripts/$script"
        [ -f "$KIT_GEN/submit-time/$step.sh" ] && kit_run_session "$step-submit" "$KIT_GEN/submit-time/$step.sh"
      else
        printf '  cp generated/learner-scripts/%s %s/scripts/\n' "$script" "$W" >&2
        [ -f "$KIT_GEN/submit-time/$step.sh" ] && { bash -n "$KIT_GEN/submit-time/$step.sh"; printf '  run generated/submit-time/%s.sh\n' "$step" >&2; }
      fi
      wd="$W/scripts"; [ $DRY = 1 ] && wd="$KIT_GEN/learner-scripts"
      rc=0
      jid=$(sbatch_or_dry "$step" "$wd" "${common[@]}" "$script") || rc=$?
      submitted_or_fail $rc
      ;;
    *) kit_die "unknown step kind $kind";;
  esac
  jid=${jid%%;*}
  JOB[$step]=$jid
  SUBMITTED=$((SUBMITTED + 1))
  printf '%s\t%s\t%s\t%s\t%s\n' "$step" "$jid" "$(kit_ts)" "${dl#:}" "$QOS" >> "$JOBS_TSV"
done < "$KIT_GEN/steps.tsv"

if [ $DRY = 1 ]; then
  kit_log "dry run: $SUBMITTED steps checked, $DRYFAIL sbatch --test-only failure(s); nothing submitted"
  [ $DRYFAIL = 0 ]
else
  kit_log "submitted $SUBMITTED step(s); job IDs in $JOBS_TSV"
  kit_log "watch: squeue -u $(id -un) -o '%.10i %.22j %.8T %.10M %R'   then: bash $KIT_DIR/collect.sh --run $(basename "$RUN")"
fi
