#!/bin/bash
# Run validate_annotation.pl over every genome of a project directory (the annotator's
# own output in work/, before any hand edit) and score the verdicts against the hand
# edits in manual_edits.tsv: every hand edit should be found with the right verdict
# (rename-to-pseudogene and add-pseudogene are REAL; delete and move are ARTIFACT),
# and every verdict that is not PASS should correspond to a hand edit.
# USAGE: validate_truth.sh <project_dir> --family F --order O [--out DIR]
set -u; HERE="$(cd "$(dirname "$0")" && pwd)"; PROJ="${1:?usage: validate_truth.sh <project_dir> --family F --order O [--out DIR]}"; shift
FAM=""; ORD=""; OUT="$PROJ/validation"
while [ $# -gt 0 ]; do case "$1" in --family) FAM="$2"; shift 2;; --order) ORD="$2"; shift 2;; --out) OUT="$2"; shift 2;; *) echo "unknown option $1" >&2; exit 1;; esac; done
mkdir -p "$OUT"; cd "$PROJ"
for d in annotations/*_annobtd; do p=${d#annotations/}; p=${p%_annobtd}; a=$d/work/${p}_VERDANT_cleaned_annotation.txt
    [ -s "$a" ] && [ -d $d/work/files ] || { echo "skip $p (no pipeline output or guides)" >&2; continue; }
    f=$(ls ${p}.fsa ${p}_FULLCP.fsa 2>/dev/null | head -1); g=$(awk -F'\t' -v p="$p" '$1==p{print $3; exit}' annotations/status.tsv 2>/dev/null); [ "$g" = NA ] && g=""
    perl "$HERE/../tools/validate_annotation.pl" "$a" "$f" --guides-dir $d/work/files --guide-seqs $d/work/sequenceFiles ${FAM:+--family "$FAM"} ${ORD:+--order "$ORD"} ${g:+--genus "$g"} --out "$OUT/$p.tsv" 2>> "$OUT/run.log"
done
perl "$HERE/score_validation.pl" "$OUT" "$PROJ/manual_edits.tsv" | tee "$OUT/score.txt"
