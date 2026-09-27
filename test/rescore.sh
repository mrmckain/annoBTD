#!/bin/bash
# Re-score an existing results directory with the current comparator, without
# re-running the pipeline. Useful when the comparator itself changes.
# USAGE: ./rescore.sh <results_dir>
set -uo pipefail
cd "$(dirname "$0")"
D="$1"
[ -d "$D" ] || { echo "no such directory: $D" >&2; exit 1; }
SUM="$D/summary.tsv"
printf "genome\ttruth\tpredicted\tmatched\trecall%%\tboth_exact%%\t5p_exact%%\t3p_exact%%\t5p_median\t3p_median\tmissing\tspurious\tseconds\n" > "$SUM"
for w in "$D"/work_*; do
    [ -d "$w" ] || continue
    g="$(basename "$w" | sed 's/^work_//')"
    r="$w/${g}_VERDANT_cleaned_annotation.txt"
    [ -s "$r" ] || continue
    perl compare_annotations.pl "data/$g.truth.tsv" "$r" \
        --label "$g" --detail "$D/$g.detail.tsv" > "$D/$g.compare.txt"
    perl summarize.pl "$D/$g.compare.txt" "$SUM" "$g" 0
done
