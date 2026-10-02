#!/bin/bash
# restage.sh: build the staged learner copy on Depot from the current staged data, then check it.
# Usage: bash restage.sh [--from SRC] [--to DEST] [--yes]
#   SRC   current staged copy (default STAGED_OLD in config.sh, the old scratch location)
#   DEST  new staged copy on Depot (default STAGED in config.sh)
#   --yes skip the confirmation prompt
# Copies only the files staged_manifest.tsv lists (so README.md, .ipynb_checkpoints/, indexes,
# results/, and job logs are left behind), regenerates SRR_Acc_List.txt (the 8 workshop runs)
# and tx2gene.tsv (Episode 04b's command on the transcript FASTA), installs the episode
# scripts from staged-scripts/, makes everything group-readable, and runs check_staged.sh.
# Never deletes anything; rerunning updates files in place.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

SRC=$STAGED_OLD; DEST=$STAGED; YES=0
while [ $# -gt 0 ]; do
  case $1 in
    --from) SRC=$2; shift 2;;
    --to) DEST=$2; shift 2;;
    --yes) YES=1; shift;;
    -h|--help) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) kit_die "unknown option $1";;
  esac
done
case $DEST in *DEPOT_PATH*) kit_die "replace DEPOT_PATH in config.sh (STAGED) with the real Depot path first";; esac
[ -d "$SRC/data" ] || kit_die "source $SRC/data not found"
[ "$(readlink -m "$SRC")" != "$(readlink -m "$DEST")" ] || kit_die "source and destination are the same"

MAN="$KIT_DIR/staged_manifest.tsv"
LIST=$(mktemp)
awk -F'\t' '$2 == "file" && $1 ~ /^data\// && $1 != "data/SRR_Acc_List.txt" && $1 != "data/tx2gene.tsv" { print $1 }' "$MAN" > "$LIST"
kit_log "restage: $(wc -l < "$LIST") data files from $SRC"
kit_log "         to $DEST (group-readable), plus SRR_Acc_List.txt, tx2gene.tsv, and $(find "$KIT_DIR/staged-scripts" -type f | wc -l) scripts"
missing=$(while read -r p; do [ -f "$SRC/$p" ] || echo "$p"; done < "$LIST")
[ -z "$missing" ] || kit_die "missing in the source: $missing"
if [ $YES = 0 ]; then
  read -r -p "Proceed? [y/N] " ans
  case $ans in y|Y|yes) ;; *) kit_die "cancelled";; esac
fi

mkdir -p "$DEST/data" "$DEST/scripts"
rsync -a --info=progress2 --files-from="$LIST" "$SRC/" "$DEST/"
rm -f "$LIST"
printf '%s\n' SRR2121786 SRR2121787 SRR2121788 SRR2121789 SRR2121778 SRR2121779 SRR2121780 SRR2121781 \
  > "$DEST/data/SRR_Acc_List.txt"
# Episode 04b command (field 1 transcript ID, field 2 gene ID of the GENCODE FASTA headers)
grep ">" "$DEST/data/gencode.vM38.transcripts.fa" | cut -c2- | cut -d "|" -f 1,2 | tr "|" "\t" > "$DEST/data/tx2gene.tsv"
cp "$KIT_DIR"/staged-scripts/* "$DEST/scripts/"
chmod -R g+rX "$DEST"
kit_log "restage: copied; checking against the manifest (full md5, a few minutes)"
bash "$KIT_DIR/check_staged.sh" "$DEST" --perm group --tree
