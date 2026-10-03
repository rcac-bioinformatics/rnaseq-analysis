#!/bin/bash
# 00_preflight.sh: check Negishi before a kit run. Login node, about 3 to 5 minutes
# (most of it waiting for a 5-minute standby job for the network test).
# Usage: bash 00_preflight.sh [--no-srun]
# Writes a PASS/WARN/FAIL table to the screen and to $KIT_SCRATCH_BASE/state/preflight-<time>.txt,
# plus state/ood_env.sh (R library setting the R jobs use) and state/pkgcheck.txt.
# Exit status 1 if any check FAILs.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

DO_SRUN=1
[ "${1:-}" = "--no-srun" ] && DO_SRUN=0
STATE="$KIT_SCRATCH_BASE/state"
mkdir -p "$STATE"
OUT="$STATE/preflight-$(date +%Y%m%d-%H%M%S).txt"
NFAIL=0; NWARN=0
row() {  # status check detail
  local s=$1; shift
  [ "$s" = FAIL ] && NFAIL=$((NFAIL + 1))
  [ "$s" = WARN ] && NWARN=$((NWARN + 1))
  printf '%-5s %-34s %s\n' "$s" "$1" "${2:-}" | tee -a "$OUT"
}
section() { printf '\n== %s ==\n' "$1" | tee -a "$OUT"; }
T0=$(date +%s)
printf 'Negishi run kit preflight, %s, user %s, host %s\n' "$(kit_ts)" "$(id -un)" "$(hostname)" | tee "$OUT"

# ---------------------------------------------------------------- kit and tools
section "Kit"
if (cd "$KIT_DIR/.." && sha256sum --quiet -c "$KIT_GEN/episodes.sha256" >/dev/null 2>&1); then
  row PASS "episodes match the kit build" "$(wc -l < "$KIT_GEN/episodes.sha256") files"
else
  row WARN "episodes differ from the kit build" "rerun: python3 $KIT_DIR/build_kit.py (episodes changed after the kit was generated)"
fi
for t in python3 sha256sum sbatch sacct sacctmgr srun rsync; do
  if command -v "$t" >/dev/null 2>&1; then row PASS "command $t" "$(command -v "$t")"; else row FAIL "command $t" "not found"; fi
done
if command -v python3 >/dev/null 2>&1; then
  pv=$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')
  if python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)'; then row PASS "python3 >= 3.6" "$pv"; else row FAIL "python3 >= 3.6" "$pv"; fi
fi
if APP=$(kit_apptainer); then row PASS "apptainer/singularity" "$APP: $("$APP" --version 2>&1 | head -1)"; else row FAIL "apptainer/singularity" "not found"; APP=; fi
if bash --noprofile --norc -c 'type module' >/dev/null 2>&1; then
  row PASS "module function exported" "available in a non-login child shell"
elif [ -n "${LMOD_PKG:-}" ] && [ -f "$LMOD_PKG/init/bash" ]; then
  row PASS "module function" "not exported; sessions source \$LMOD_PKG/init/bash"
else
  row WARN "module function" "not exported and LMOD_PKG unset; sessions fall back to /etc/profile.d/lmod.sh"
fi
for v in SCRATCH RCAC_SCRATCH; do
  if [ "${!v:-}" = "/scratch/negishi/$(id -un)" ]; then row PASS "\$$v" "${!v}"
  else row WARN "\$$v" "'${!v:-unset}' (episodes assume /scratch/negishi/\$USER)"; fi
done
case "$KIT_SCRATCH_BASE" in /scratch/negishi/*) row PASS "KIT_SCRATCH_BASE" "$KIT_SCRATCH_BASE";;
  *) row FAIL "KIT_SCRATCH_BASE" "$KIT_SCRATCH_BASE must be under /scratch/negishi/";; esac

# ---------------------------------------------------------------- modules
section "Modules (module load biocontainers, then the tool; default version must match)"
declare -A VCMD=( [fastqc]="fastqc --version" [multiqc]="multiqc --version" [fastp]="fastp --version"
  [star]="STAR --version" [subread]="featureCounts -v" [kallisto]="kallisto version"
  [salmon]="salmon --version" [sra-tools]="fasterq-dump --version" [seqtk]="true" )
for mv in $EXPECTED_MODULES; do
  tool=${mv%%/*}
  res=$(timeout 120 bash -c "module load biocontainers >/dev/null 2>&1; module load $tool >/dev/null 2>&1 || exit 3; \
        module -t list 2>&1 | grep -E '^$tool/' | head -1; ${VCMD[$tool]} 2>&1 | grep -v -i warning | head -2 | tr '\n' ' '" 2>&1) || {
    row FAIL "module $tool" "did not load"; continue; }
  loaded=$(printf '%s\n' "$res" | head -1)
  if [ "$loaded" = "$mv" ]; then row PASS "module $tool" "$loaded; $(printf '%s\n' "$res" | sed -n 2p)"
  else row FAIL "module $tool" "default is '$loaded', expected $mv (episodes load the default)"; fi
done
res=$(timeout 300 bash -c 'module load biocontainers >/dev/null 2>&1; module load r-rnaseq >/dev/null 2>&1 || exit 3
  module -t list 2>&1 | grep -E "^r-rnaseq/" | head -1
  R --no-save --no-restore --quiet -e "for (p in readLines(\"'"$KIT_GEN"'/r-rnaseq_packages.txt\")) cat(p, tryCatch(as.character(packageVersion(p)), error = function(e) \"MISSING\"), \"\n\"); cat(R.version.string, \"\n\")" 2>&1 | grep -E "^[A-Za-z0-9.]+ ([0-9.-]+|MISSING) *$|^R version"' 2>&1) || res="LOADFAIL"
if [ "$res" = LOADFAIL ]; then row FAIL "module r-rnaseq (04b tximport)" "did not load"
elif printf '%s' "$res" | grep -q MISSING; then row FAIL "module r-rnaseq (04b tximport)" "$(printf '%s' "$res" | tr '\n' ' ')"
else row PASS "module r-rnaseq (04b tximport)" "$(printf '%s' "$res" | tr '\n' ' ')"; fi

# ---------------------------------------------------------------- SLURM
section "SLURM account and QoS"
assoc=$(sacctmgr -nP show assoc where user="$(id -un)" format=account,qos 2>/dev/null | grep -i "^$KIT_ACCOUNT|" || true)
if [ -n "$assoc" ]; then row PASS "member of $KIT_ACCOUNT" "qos: ${assoc#*|}"; else row FAIL "member of $KIT_ACCOUNT" "no association for $(id -un)"; fi
for q in standby normal; do
  if printf '%s' "$assoc" | grep -qw "$q"; then row PASS "QoS $q on $KIT_ACCOUNT" "$(sacctmgr -nP show qos "$q" format=name,maxwall 2>/dev/null | head -1)"
  elif [ "$q" = standby ]; then row FAIL "QoS standby on $KIT_ACCOUNT" "not in the association"
  else row WARN "QoS normal on $KIT_ACCOUNT" "not available; use --qos standby"; fi
done
TMPJ=$(mktemp -d "${TMPDIR:-/tmp}/kitpf.XXXXXX")
for f in "$KIT_GEN"/learner-scripts/*.sh "$KIT_GEN"/jobs/*.sh; do
  r=$(cd "$TMPJ" && sbatch --test-only "$f" 2>&1 | tail -1) || true
  case "$r" in *"to start at"*) row PASS "sbatch --test-only $(basename "$f")" "${r#sbatch: }";;
    *) row FAIL "sbatch --test-only $(basename "$f")" "$r";; esac
done
rmdir "$TMPJ" 2>/dev/null || true

# ---------------------------------------------------------------- scratch space
section "Scratch space"
staged_gb=$(du -s --block-size=1G "$STAGED" 2>/dev/null | cut -f1 || echo 0)
need=$(( staged_gb + 60 ))
avail=$(df --output=avail -B1G "$KIT_SCRATCH_BASE" 2>/dev/null | tail -1 | tr -d ' ' || echo 0)
detail="staged ${staged_gb} GB; one learner run needs about ${need} GB (copy + indexes + BAMs); regen needs as much again plus the copy into the dated results directory"
if command -v myquota >/dev/null 2>&1; then
  row INFO "myquota" "$(myquota 2>/dev/null | grep -i scratch | head -1 | tr -s ' ')"
fi
if [ "${avail:-0}" -ge $(( need * 2 )) ]; then row PASS "filesystem free space" "${avail} GB free; $detail"
else row WARN "filesystem free space" "${avail} GB free; $detail (check your quota with myquota)"; fi

# ---------------------------------------------------------------- staged data
section "Staged data ($STAGED)"
case $STAGED in *DEPOT_PATH*) row FAIL "staged path" "STAGED still contains the placeholder DEPOT_PATH; set it in config.sh, learners/setup.md, .claude/CLAUDE.md";; esac
if [ "$STAGED_PERM" = group ]; then RPOS=4; XPOS=6; else RPOS=7; XPOS=9; fi
path_ok=1
d=$STAGED
while [ "$d" != / ]; do
  perm=$(stat -c %A "$d" 2>/dev/null || echo "----------")
  case "${perm:$XPOS:1}" in x|s|t) ;; *) row FAIL "path traversable by $STAGED_PERM" "$d is $perm"; path_ok=0;; esac
  d=$(dirname "$d")
done
[ $path_ok = 1 ] && row PASS "path traversable by $STAGED_PERM" "$STAGED"
nbad=0
while IFS=$'\t' read -r p typ bytes md5 need src; do
  case "$p" in path|\#*) continue;; esac
  f="$STAGED/$p"
  if [ "$typ" = dir ]; then
    if [ -d "$f" ]; then row PASS "dir $p" "$(stat -c %A "$f")"
    elif [ "$need" = optional ]; then row WARN "dir $p" "absent ($src)"
    else row FAIL "dir $p" "missing"; nbad=$((nbad + 1)); fi
    continue
  fi
  if [ ! -f "$f" ]; then row FAIL "file $p" "missing"; nbad=$((nbad + 1)); continue; fi
  perm=$(stat -c %A "$f")
  [ "${perm:$RPOS:1}" = r ] || { row FAIL "file $p" "not readable by $STAGED_PERM ($perm)"; nbad=$((nbad + 1)); continue; }
  sz=$(stat -c %s "$f")
  if [ "$bytes" != - ] && [ "$sz" != "$bytes" ]; then row FAIL "file $p" "size $sz, expected $bytes ($src)"; nbad=$((nbad + 1)); continue; fi
  # md5 only for small files here (run check_staged.sh for the full check)
  if [ "$md5" != - ] && [ "$bytes" -lt 50000000 ] && [ "$(md5sum < "$f" | cut -c1-32)" != "$md5" ]; then
    row FAIL "file $p" "md5 differs from the manifest ($src)"; nbad=$((nbad + 1)); continue; fi
  if [ "$p" = data/SRR_Acc_List.txt ]; then
    n=$(grep -c . "$f"); [ "$n" = 8 ] || { row FAIL "file $p" "$n lines, expected 8"; continue; }
  fi
done < "$KIT_DIR/staged_manifest.tsv"
[ $nbad = 0 ] && row PASS "staged manifest" "all required files present, readable, expected size; small files md5-checked (full check: check_staged.sh)"
pre=$(kit_precomputed "$STAGED" | tr '\n' ' ')
if [ -z "$pre" ]; then row PASS "no precomputed outputs staged" "learners start from FASTQ and references, as setup.md says"
else row FAIL "precomputed outputs in staged copy" "$pre(learners would start with finished results; the kit would skip steps)"; fi
if [ "$STAGED_PERM" = group ]; then unread=$(find "$STAGED" ! -perm -g+r 2>/dev/null | head -5 | tr '\n' ' ')
else unread=$(find "$STAGED" ! -perm -o+r 2>/dev/null | head -5 | tr '\n' ' '); fi
if [ -z "$unread" ]; then row PASS "everything readable by $STAGED_PERM" "learners must be members of the Depot group"
else row FAIL "unreadable by $STAGED_PERM" "$unread (learner rsync fails on these)"; fi
for f in "$KIT_DIR"/staged-scripts/*; do
  b=$(basename "$f")
  if [ ! -f "$STAGED/scripts/$b" ]; then row WARN "staged scripts/$b" "absent; copy negishi-run/staged-scripts/$b"
  elif cmp -s "$f" "$STAGED/scripts/$b"; then row PASS "staged scripts/$b" "matches the episodes"
  else row WARN "staged scripts/$b" "differs from the episodes; replace with negishi-run/staged-scripts/$b"; fi
done
if [ -f "$STAGED/data/tx2gene.tsv" ]; then
  ntx=$(grep -c '>' "$STAGED/data/gencode.vM38.transcripts.fa" 2>/dev/null || echo 0)
  nmap=$(wc -l < "$STAGED/data/tx2gene.tsv")
  if [ "$ntx" = "$nmap" ]; then row PASS "staged data/tx2gene.tsv" "$nmap transcripts (complete)"
  else row WARN "staged data/tx2gene.tsv" "$nmap rows vs $ntx transcripts: old GTF-based map; replace (episode 04b now builds it from the FASTA)"; fi
fi

section "Completed results ($STAGED_RESULTS)"
if [ -d "$STAGED_RESULTS" ]; then
  for p in data/star_index/SA data/star_index/Genome data/salmon_index data/kallisto_index/transcripts.idx \
           results/qc_fastq results/counts/gene_counts_clean.txt results/kallisto_quant/txi.rds \
           results/deseq2/DESeq2_results_joined.tsv; do
    if [ -e "$STAGED_RESULTS/$p" ]; then row PASS "results $p" ""
    else row WARN "results $p" "missing (instructor-notes/README.md lists it; regen_results.sh rebuilds it)"; fi
  done
  nb=$(find "$STAGED_RESULTS/results/mapping" -maxdepth 1 -name '*.bam' 2>/dev/null | wc -l)
  if [ "$nb" = 8 ]; then row PASS "results mapping BAMs" "8"; else row WARN "results mapping BAMs" "$nb of 8"; fi
  [ -e "$STAGED_RESULTS/data/star_index/SA" ] || row WARN "04a callout prebuilt index" "episode 04a tells learners to ln -s $STAGED_RESULTS/data/star_index; it is missing"
  if [ "$STAGED_PERM" = group ]; then unread=$(find "$STAGED_RESULTS" ! -perm -g+r 2>/dev/null | head -5 | tr '\n' ' ')
  else unread=$(find "$STAGED_RESULTS" ! -perm -o+r 2>/dev/null | head -5 | tr '\n' ' '); fi
  if [ -z "$unread" ]; then row PASS "results readable by $STAGED_PERM" ""
  else row WARN "results unreadable by $STAGED_PERM" "$unread (learners copy single files and ln -s the STAR index)"; fi
else
  row WARN "completed results directory" "missing; setup.md and episode 04a point learners to it"
fi
PB_STAR=""; PB_KAL=""
if [ -e "$STAGED_RESULTS/data/star_index/SA" ]; then PB_STAR="$STAGED_RESULTS/data/star_index"; fi
if [ -e "$STAGED_RESULTS/data/kallisto_index/transcripts.idx" ]; then PB_KAL="$STAGED_RESULTS/data/kallisto_index"; fi
printf 'PREBUILT_STAR=%q\nPREBUILT_KALLISTO=%q\n' "$PB_STAR" "$PB_KAL" > "$STATE/prebuilt.sh"

# ---------------------------------------------------------------- OOD image
section "Open OnDemand R image"
if [ -r "$OOD_SIF" ]; then row PASS "SIF readable" "$OOD_SIF ($(du -h "$OOD_SIF" | cut -f1))"; else row FAIL "SIF readable" "$OOD_SIF"; fi
if [ -d "$OOD_HOST_LIB" ] && [ -r "$OOD_HOST_LIB" ]; then
  row PASS "host library readable" "$OOD_HOST_LIB ($(find "$OOD_HOST_LIB" -mindepth 1 -maxdepth 1 -type d | wc -l) packages)"
else row FAIL "host library readable" "$OOD_HOST_LIB"; fi
: > "$STATE/ood_env.sh"
if [ -n "$APP" ] && [ -r "$OOD_SIF" ]; then
  RCHK='pk <- readLines(Sys.getenv("KIT_PK")); ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
cat("R", as.character(getRversion()), "Bioc", tryCatch(as.character(packageVersion("BiocVersion")), error = function(e) "?"), "\n")
cat("HOSTLIB", any(grepl("host-site-library", .libPaths())), "\n")
cat("LOADED", sum(ok), "of", length(ok), "\n"); cat("MISSING", paste(pk[!ok], collapse = ","), "\n")
for (p in c("DESeq2", "clusterProfiler", "msigdbr", "biomaRt", "enrichplot")) cat("VERSION", p, tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"), "\n")
u <- tryCatch({ b <- paste(deparse(get("check_cache", asNamespace("msigdbr"))), collapse = ""); regmatches(b, regexpr("https?://[^\"]+", b)) }, error = function(e) "")
cat("MSIGDB_URL", if (length(u)) u else "", "\n")'
  run_r() {  # $1: extra --env entries
    timeout 300 "$APP" exec --cleanenv --home "$STATE/ood-home" --bind "$OOD_HOST_LIB:$OOD_HOST_LIB_MOUNT:ro" --bind "$KIT_DIR" \
      --env "KIT_PK=$KIT_GEN/ood_packages.txt,R_LIBS_USER=$STATE/empty-rlib$1" "$OOD_SIF" Rscript -e "$RCHK" 2>/dev/null
  }
  mkdir -p "$STATE/ood-home" "$STATE/empty-rlib"
  r1=$(run_r '' || true)
  if ! printf '%s' "$r1" | grep -q "^HOSTLIB TRUE"; then
    r2=$(run_r ",R_LIBS_SITE=$OOD_HOST_LIB_MOUNT" || true)
    if printf '%s' "$r2" | grep -q "^HOSTLIB TRUE"; then
      echo "OOD_R_LIBS_SITE=$OOD_HOST_LIB_MOUNT" > "$STATE/ood_env.sh"
      row WARN "host library on .libPaths()" "only with R_LIBS_SITE=$OOD_HOST_LIB_MOUNT; R jobs will set it (check the OOD app does the same)"
      r1=$r2
    else
      row FAIL "host library on .libPaths()" "$OOD_HOST_LIB_MOUNT not used by R in the image"
    fi
  else
    row PASS "host library on .libPaths()" "with --cleanenv and the bind only"
  fi
  rv=$(printf '%s\n' "$r1" | grep '^R ' || true)
  case "$rv" in "R $EXPECTED_R_VERSION Bioc $EXPECTED_BIOC_VERSION"*) row PASS "R and Bioconductor" "$rv";;
    *) row WARN "R and Bioconductor" "'$rv', expected R $EXPECTED_R_VERSION Bioc $EXPECTED_BIOC_VERSION";; esac
  ld=$(printf '%s\n' "$r1" | grep '^LOADED' || echo "LOADED ?"); ms=$(printf '%s\n' "$r1" | grep '^MISSING' | cut -c9- || true)
  if [ -z "$(printf '%s' "$ms" | tr -d ' ')" ] && printf '%s' "$ld" | grep -q 'LOADED'; then row PASS "episode packages" "${ld#LOADED }: $(tr '\n' ' ' < "$KIT_GEN/ood_packages.txt")"
  else row FAIL "episode packages" "${ld#LOADED }; missing: $ms"; fi
  row INFO "package versions" "$(printf '%s\n' "$r1" | grep '^VERSION' | cut -c9- | tr '\n' ';')"
  printf '%s: %s; R/Bioc: %s; missing: %s\n' "$(date +%F)" "${ld#LOADED }" "${rv#R }" "${ms:-none}" > "$STATE/pkgcheck.txt"
  MSIGDB_URL=$(printf '%s\n' "$r1" | grep '^MSIGDB_URL' | cut -d' ' -f2 || true)
fi
MSIGDB_URL=${MSIGDB_URL:-https://zenodo.org/records/18968178/files/msigdb.2026.1.zip}

# ---------------------------------------------------------------- compute-node egress
section "Network from a compute node (5-minute standby job; waits at most 4 minutes to start)"
if [ $DO_SRUN = 1 ]; then
  # shellcheck disable=SC2016  # expanded by the bash on the compute node, not here
  EG='for u in "$@"; do printf "%s %s\n" "$(curl -s -o /dev/null -r 0-0 -L --max-time 30 -w "%{http_code}" "$u" || echo 000)" "$u"; done'
  eg=$(timeout 240 srun -A "$KIT_ACCOUNT" -q "$KIT_QOS_DEFAULT" -p "$KIT_PARTITION" -N1 -n1 -t 00:05:00 --job-name=kit-egress \
       bash -c "$EG" _ https://rest.kegg.jp/info/kegg \
       "https://www.ensembl.org/biomart/martservice?type=registry" \
       "https://useast.ensembl.org/biomart/martservice?type=registry" \
       "$MSIGDB_URL" \
       https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/ \
       "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRR2121778" 2>&1) || true
  if [ -z "$eg" ] || ! printf '%s' "$eg" | grep -qE '^[0-9]{3} '; then
    row WARN "compute-node network test" "job did not start or failed within 4 min (standby queue?): $(printf '%s' "$eg" | tail -1)"
  else
    while read -r code url; do
      case "$code" in 2??|3??) row PASS "egress $(printf '%s' "$url" | cut -d/ -f3)" "HTTP $code";;
        *) row FAIL "egress $(printf '%s' "$url" | cut -d/ -f3)" "HTTP $code $url";; esac
    done < <(printf '%s\n' "$eg" | grep -E '^[0-9]{3} ')
  fi
else
  row WARN "compute-node network test" "skipped (--no-srun)"
fi

# ---------------------------------------------------------------- estimates
section "Estimates (est_min = local measurement x1.5; queue wait not included)"
awk -F'\t' 'NR>1 { ch += $6 * $7 * $10 / 60.0; printf "      %-14s %3d cpu x %d task(s) x %3d min\n", $1, $6, $7, $10 } END { printf "      total about %.0f core-hours per full run (both tracks); regen_results.sh is the same again\n", ch }' \
  "$KIT_GEN/steps.tsv" | tee -a "$OUT"
cat <<'EOF' | tee -a "$OUT"
      critical path (genome track): setup copy ~15 min, STAR index ~65, mapping ~18, featureCounts ~6,
      clean ~5, episode 05 ~10, episode 06 ~15: about 2.2 h of run time plus standby queue waits.
      Transcript track in parallel: kallisto index ~35, quant ~12, tximport ~6, 05b ~6, 06 ~15: about 1.3 h.
EOF

section "Result"
printf 'FAIL %d, WARN %d, %d s. Saved to %s\n' "$NFAIL" "$NWARN" "$(( $(date +%s) - T0 ))" "$OUT" | tee -a "$OUT"
[ "$NFAIL" -eq 0 ]
