#!/bin/bash
# check_staged.sh: check a staged learner copy against staged_manifest.tsv.
# Usage: bash check_staged.sh [DIR] [--quick] [--tree] [--perm group|other]
#   DIR      staged copy to check (default: STAGED from config.sh)
#   --perm   who must be able to read it: group (Depot, default from STAGED_PERM) or other
#   --quick  md5 only for files under 50 MB (sizes are always checked); the full check
#            reads about 20 GB and takes a few minutes
#   --tree   also print the data/ listing that Episode 02 shows
# Checks: every expected file present with the expected size and md5, readable by the group
# (or by others), directories traversable likewise, and nothing extra (indexes, results/,
# job logs, README.md, .ipynb_checkpoints).
# Exit status 1 if anything fails. Read only: changes nothing.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

DIR=$STAGED; QUICK=0; TREE=0; PERM=$STAGED_PERM; next_perm=0
for a in "$@"; do
  if [ $next_perm = 1 ]; then PERM=$a; next_perm=0; continue; fi
  case $a in
    --quick) QUICK=1;;
    --tree) TREE=1;;
    --perm) next_perm=1;;
    -h|--help) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) DIR=$a;;
  esac
done
MAN="$KIT_DIR/staged_manifest.tsv"
case $PERM in group) RPOS=4; XPOS=6;; other) RPOS=7; XPOS=9;; *) kit_die "--perm must be group or other";; esac
case $DIR in *DEPOT_PATH*) kit_die "replace DEPOT_PATH in config.sh (STAGED) with the real Depot path first";; esac
[ -d "$DIR" ] || kit_die "$DIR is not a directory"
bad=0
fail() { printf 'FAIL  %s\n' "$*"; bad=$((bad + 1)); }
ok() { printf 'ok    %s\n' "$*"; }

expected=$(mktemp)
while IFS=$'\t' read -r p typ bytes md5 _need _src; do
  case "$p" in path|\#*) continue;; esac
  echo "$p" >> "$expected"
  f="$DIR/$p"
  if [ "$typ" = dir ]; then
    if [ -d "$f" ]; then ok "$p/"; else fail "$p/ missing"; fi
    continue
  fi
  [ -f "$f" ] || { fail "$p missing"; continue; }
  perm=$(stat -L -c %A "$f")
  [ "${perm:$RPOS:1}" = r ] || { fail "$p not readable by $PERM ($perm)"; continue; }
  sz=$(stat -L -c %s "$f")
  [ "$sz" = "$bytes" ] || { fail "$p size $sz, expected $bytes"; continue; }
  if [ $QUICK = 1 ] && [ "$bytes" -gt 50000000 ]; then ok "$p (size only)"; continue; fi
  m=$(md5sum < "$f" | cut -c1-32)
  if [ "$m" = "$md5" ]; then ok "$p"; else fail "$p md5 $m, expected $md5"; fi
done < "$MAN"

# directories inside the copy, and on the path to it, must be traversable
while IFS= read -r d; do
  perm=$(stat -c %A "$d")
  case "${perm:$XPOS:1}" in x|s|t) ;; *) fail "directory $d not traversable by $PERM ($perm)";; esac
done < <(find "$DIR" -type d)
d=$(readlink -f "$DIR")
while [ "$d" != / ]; do
  perm=$(stat -c %A "$d")
  case "${perm:$XPOS:1}" in x|s|t) ;; *) fail "directory $d not traversable by $PERM ($perm)";; esac
  d=$(dirname "$d")
done
# anything not in the manifest
extra=$(cd "$DIR" && find . -mindepth 1 \( -type f -o -type l \) | sed 's#^\./##' | sort | comm -23 - <(sort "$expected"))
if [ -n "$extra" ]; then
  printf 'FAIL  not in the manifest (learners would receive these):\n'; printf '%s\n' "$extra" | head -40 | sed 's/^/        /'
  bad=$((bad + 1))
fi
rm -f "$expected"

if [ $TREE = 1 ]; then
  printf '\n# Episode 02 listing of %s/data\n' "$DIR"
  (cd "$DIR" && { tree data 2>/dev/null || ls -1 data; })
fi
printf '\n%s: %d problem(s)\n' "$DIR" "$bad"
[ "$bad" -eq 0 ]
