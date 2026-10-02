# negishi-run/lib.sh: functions shared by the kit scripts and the generated jobs.
# Sourced by scripts that already run under `set -euo pipefail`.
# shellcheck shell=bash

KIT_DIR=${KIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
export KIT_DIR
# shellcheck source=negishi-run/config.sh
source "$KIT_DIR/config.sh"
KIT_GEN=${KIT_GEN:-$KIT_DIR/$KIT_GEN_DIRNAME}
export KIT_GEN

kit_ts() { date '+%Y-%m-%dT%H:%M:%S'; }
kit_log() { printf '[%s] %s\n' "$(kit_ts)" "$*" >&2; }
kit_die() { kit_log "ERROR: $*"; exit 1; }

kit_apptainer() {
  if command -v apptainer >/dev/null 2>&1; then command -v apptainer
  elif command -v singularity >/dev/null 2>&1; then command -v singularity
  else return 1; fi
}

# Refuse to run if the learner working directory could be the staged source.
kit_guard() {
  local w=${W:?W not set}
  # KIT_ALLOW_ANY_ROOT=1 is only for testing the kit off-cluster
  case "$w" in /scratch/negishi/*) ;; *) [ "${KIT_ALLOW_ANY_ROOT:-0}" = 1 ] || kit_die "W=$w is not under /scratch/negishi/";; esac
  case "$(readlink -m "$w")/" in
    "$(readlink -m "$STAGED")"/*|"$(readlink -m "$STAGED_RESULTS")"/*)
      kit_die "W=$w resolves inside the staged data; refusing to run";;
  esac
  [ "$(readlink -m "$w")" != "$(readlink -m "$STAGED")" ] || kit_die "W is the staged directory"
  case "$(readlink -m "$w")" in "$(readlink -m "$KIT_SCRATCH_BASE")"/*) ;;
    *) kit_die "W=$w is not under KIT_SCRATCH_BASE=$KIT_SCRATCH_BASE";; esac
}

# Per-run settings frozen by submit_all.sh (or 01_learner_setup.sh) in $RUN/kit.env
kit_load_run() {
  : "${RUN:?RUN not set}"
  [ -f "$RUN/kit.env" ] || kit_die "$RUN/kit.env missing; start with 01_learner_setup.sh"
  # shellcheck source=/dev/null
  source "$RUN/kit.env"
  export RUN RUN_ID W W_T KIT_REAL_USER KIT_LEARNER_SCRATCH KIT_USER_OVERRIDE KIT_USER_OVERRIDE_T STAGED STAGED_RESULTS
  export SCRATCH="$KIT_LEARNER_SCRATCH" RCAC_SCRATCH="$KIT_LEARNER_SCRATCH"
  kit_guard
}

# Write $RUN/kit.env for a new run directory
kit_init_run() {
  local run=$1
  case "$run" in /scratch/negishi/*) ;; *) [ "${KIT_ALLOW_ANY_ROOT:-0}" = 1 ] || kit_die "run directory $run must be under /scratch/negishi/";; esac
  mkdir -p "$run"/{scratch,scratch_t,logs,records/listings,records/console,records/events,records/figures,records/sessioninfo,markers,tmp,ood}
  {
    printf 'RUN=%q\n' "$run"
    printf 'RUN_ID=%q\n' "$(basename "$run")"
    printf 'W=%q\n' "$run/scratch/rnaseq-workshop"
    printf 'W_T=%q\n' "$run/scratch_t/rnaseq-workshop"
    printf 'KIT_REAL_USER=%q\n' "$(id -un)"
    printf 'KIT_LEARNER_SCRATCH=%q\n' "$run/scratch"
    printf 'KIT_USER_OVERRIDE=%q\n' "${run#/scratch/negishi/}/scratch"
    printf 'KIT_USER_OVERRIDE_T=%q\n' "${run#/scratch/negishi/}/scratch_t"
    printf 'KIT_DIR=%q\n' "$KIT_DIR"
    printf 'KIT_GEN=%q\n' "$KIT_GEN"
    printf 'KIT_BUILT_FROM=%q\n' "$(sha256sum "$KIT_GEN/episodes.sha256" | cut -c1-16)"
  } > "$run/kit.env"
}

kit_marker() { printf '%s %s\n' "$2" "$(kit_ts)" > "$RUN/markers/$1.done"; }
# done = marker says PASS or WARN (a step that is running or failed says RUNNING)
kit_is_done() { [ -f "$RUN/markers/$1.done" ] && grep -qE '^(PASS|WARN) ' "$RUN/markers/$1.done"; }

kit_on_exit() {
  local rc=$?
  if [ "$rc" -ne 0 ] && [ -n "${KIT_STEP:-}" ] && [ -n "${RUN:-}" ]; then
    kit_marker "$KIT_STEP" FAIL
    printf '%s\t%s\t%s\t%s\n' "$KIT_STEP" "${SLURM_JOB_ID:-none}" "$(kit_ts)" "FAIL" >> "$RUN/records/step-end.tsv"
  fi
}

# Precomputed outputs in a learner copy: listed on stdout (empty = clean)
kit_precomputed() {
  local w=$1 p
  for p in data/star_index data/salmon_index data/kallisto_index; do [ -e "$w/$p" ] && echo "$p"; done
  [ -d "$w/results" ] && find "$w/results" -mindepth 1 -maxdepth 1 2>/dev/null | sed "s#^$w/##" | head -20
  find "$w/scripts" -maxdepth 1 -name 'cluster-*' 2>/dev/null | head -3 | sed "s#^$w/##"
  return 0
}

kit_job_begin() {
  local name=$1
  kit_load_run
  export KIT_STEP=$name
  trap kit_on_exit EXIT
  kit_log "step $name begin: job ${SLURM_JOB_ID:-none} on $(hostname) cpus=${SLURM_CPUS_ON_NODE:-?}"
  kit_marker "$name" RUNNING
  printf '%s\t%s\t%s\n' "$name" "${SLURM_JOB_ID:-none}" "$(kit_ts)" >> "$RUN/records/step-start.tsv"
}

kit_job_end() {
  local name=$1 status=${2:-PASS}
  kit_marker "$name" "$status"
  printf '%s\t%s\t%s\t%s\n' "$name" "${SLURM_JOB_ID:-none}" "$(kit_ts)" "$status" >> "$RUN/records/step-end.tsv"
  kit_log "step $name end: $status"
}

# Run a generated learner session in a child bash without -e/-u/pipefail (as a learner's
# shell behaves). Any command that fails is logged as "### KIT-ERR" by the session's ERR
# trap; the step fails if there is at least one.
kit_run_session() {
  # one log per attempt (session-<step>.<jobid>.log); session-<step>.log is a copy of the latest
  local name=$1 file=$2 log="$RUN/logs/session-$1.${SLURM_JOB_ID:-$$}.log" rc=0
  kit_log "session $name: $file -> $log"
  ( cd "$W" 2>/dev/null || cd "$KIT_LEARNER_SCRATCH"; bash --noprofile --norc "$file" ) > "$log" 2>&1 || rc=$?
  cp "$log" "$RUN/logs/session-$name.log"
  cp "$log" "$RUN/records/console/session-$name.log"
  local nerr
  nerr=$(grep -c '^### KIT-ERR' "$log" || true)
  if [ "$rc" -ne 0 ] || [ "$nerr" -gt 0 ]; then
    kit_log "session $name: exit $rc, $nerr failed command(s); first:"
    grep -m1 '^### KIT-ERR' "$log" >&2 || true
    exit 1
  fi
}

# OOD R episode: fresh R session in the OOD image, as the RStudio app runs it
kit_ood_run() {
  local name=$1 plan=$2 user_path=$3 script=${4:-$KIT_DIR/R/run_blocks.R}
  local app home rlib rc=0 extra=""
  app=$(kit_apptainer) || kit_die "apptainer/singularity not found"
  # fresh, empty home and user library for every attempt, like a new learner
  home="$RUN/ood/home/$name.${SLURM_JOB_ID:-$$}"; rlib="$RUN/ood/rlib/$name.${SLURM_JOB_ID:-$$}"
  mkdir -p "$home" "$rlib"
  if [ -f "$KIT_SCRATCH_BASE/state/ood_env.sh" ]; then
    # shellcheck source=/dev/null
    source "$KIT_SCRATCH_BASE/state/ood_env.sh"   # sets OOD_R_LIBS_SITE if preflight found it necessary
  fi
  [ -n "${OOD_R_LIBS_SITE:-}" ] && extra=",R_LIBS_SITE=$OOD_R_LIBS_SITE"
  kit_log "R step $name: $app exec --cleanenv $OOD_SIF (USER=$user_path)"
  local binds=(--bind "$OOD_HOST_LIB:$OOD_HOST_LIB_MOUNT:ro" --bind "$KIT_SCRATCH_BASE" --bind "$KIT_DIR")
  local x; for x in ${KIT_EXTRA_BINDS:-}; do binds+=(--bind "$x"); done   # off-cluster testing only
  "$app" exec --cleanenv --home "$home" "${binds[@]}" \
    --env "USER=$user_path,R_LIBS_USER=$rlib,KIT_RECORDS=$RUN/records,KIT_STEP=$name,KIT_GEN=$KIT_GEN,W=/scratch/negishi/$user_path/rnaseq-workshop,KIT_RUN_START=$(stat -c %Y "$RUN/kit.env")$extra" \
    "$OOD_SIF" Rscript "$script" "$plan" > "$RUN/logs/R-$name.${SLURM_JOB_ID:-$$}.log" 2>&1 || rc=$?
  cp "$RUN/logs/R-$name.${SLURM_JOB_ID:-$$}.log" "$RUN/logs/R-$name.log"
  cp "$RUN/logs/R-$name.log" "$RUN/records/console/R-$name.log"
  case $rc in
    0) kit_job_end "$name" PASS ;;
    3) kit_job_end "$name" WARN ;;   # only optional (spoiler) blocks failed
    *) kit_log "R step $name failed (exit $rc); see $RUN/logs/R-$name.log"; exit 1 ;;
  esac
}

# ---- helpers used by "KIT EXTRA" lines in the generated sessions (exported below)
kit_compare_refs() {  # fresh_dir staged_dir -> TSV on stdout
  local fresh=$1 staged=$2 f
  printf 'file\tfresh_bytes\tstaged_bytes\tfresh_md5\tstaged_md5\tmatch\n'
  for f in GRCm39.primary_assembly.genome.fa gencode.vM38.primary_assembly.basic.annotation.gtf \
           gencode.vM38.transcripts.fa gencode.vM38.transcripts-clean.fa; do
    local a=- b=- ma=- mb=-
    [ -f "$fresh/$f" ] && a=$(stat -c %s "$fresh/$f") && ma=$(md5sum < "$fresh/$f" | cut -c1-32)
    [ -f "$staged/$f" ] && b=$(stat -c %s "$staged/$f") && mb=$(md5sum < "$staged/$f" | cut -c1-32)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$a" "$b" "$ma" "$mb" "$([ "$ma" = "$mb" ] && [ "$ma" != - ] && echo yes || echo NO)"
  done
}

kit_subsample_inputs() {  # dir: tiny copies (10,000 reads) of the staged FASTQs
  local d=$1 s r
  mkdir -p "$d"
  for s in WT_Bcell_mock_rep1 WT_Bcell_mock_rep2 WT_Bcell_mock_rep3 WT_Bcell_mock_rep4 \
           WT_Bcell_IR_rep1 WT_Bcell_IR_rep2 WT_Bcell_IR_rep3 WT_Bcell_IR_rep4; do
    for r in R1 R2; do
      { zcat "$STAGED/data/${s}_${r}.fastq.gz" || true; } | head -n 40000 | gzip > "$d/${s}_${r}.fastq.gz"
    done
  done
}

kit_check_subsample() {  # dir: read counts and pairing after the spoiler
  local d=$1 s n1 n2
  printf 'sample\tR1_reads\tR2_reads\tpaired_names_match\toriginals_archived\n'
  for s in WT_Bcell_mock_rep1 WT_Bcell_mock_rep2 WT_Bcell_mock_rep3 WT_Bcell_mock_rep4 \
           WT_Bcell_IR_rep1 WT_Bcell_IR_rep2 WT_Bcell_IR_rep3 WT_Bcell_IR_rep4; do
    n1=$(( $(zcat "$d/${s}_R1.fastq.gz" 2>/dev/null | wc -l) / 4 ))
    n2=$(( $(zcat "$d/${s}_R2.fastq.gz" 2>/dev/null | wc -l) / 4 ))
    local same=NO
    if cmp -s <(zcat "$d/${s}_R1.fastq.gz" 2>/dev/null | awk 'NR%4==1{print $1}' | sed 's#/1$##') \
              <(zcat "$d/${s}_R2.fastq.gz" 2>/dev/null | awk 'NR%4==1{print $1}' | sed 's#/2$##'); then same=yes; fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$n1" "$n2" "$same" \
      "$([ -f "$d/original_reads/${s}_R1.fastq.gz" ] && echo yes || echo NO)"
  done
}

kit_multiqc_export() {  # input_dir out_dir: MultiQC with --export into a kit-owned directory
  local in=$1 out=$2
  mkdir -p "$out"
  module load biocontainers
  module load multiqc
  multiqc --export --force "$in" -o "$out" > "$out/multiqc.log" 2>&1 || echo "### KIT-ERR multiqc export failed" >&2
}

kit_fastp_test() {
  module load biocontainers
  bash --noprofile --norc "$KIT_GEN/sessions/03-fastp.sh" > "$RUN/logs/fastp.log" 2>&1 \
    || echo "### KIT-ERR fastp spoiler test failed (see logs/fastp.log)" >&2
}

export -f kit_ts kit_log kit_compare_refs kit_subsample_inputs kit_check_subsample kit_multiqc_export kit_fastp_test

# 06 on the transcript track runs in its own light learner directory, so its
# results/enrichment does not collide with the genome-track run
kit_transcript_dir() {
  mkdir -p "$W_T/results" "$W_T/data" "$W_T/scripts"
  ln -sfn "$W/results/deseq2_kallisto" "$W_T/results/deseq2_kallisto"
}
