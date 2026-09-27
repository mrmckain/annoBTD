#!/bin/bash
# Download the guide genomes named in a select_guides.pl ranking and lay them out
# the way run_annotation_pipeline.sh expects.
#
# USAGE: fetch_guides.sh <ranking.tsv> <seqdir> <filedir> [cachedir]
#
# ranking.tsv is the output of select_guides.pl: accession, organism, family,
# order, jaccard, rank. For each accession this writes:
#   <seqdir>/<name>   the plastome FASTA
#   <filedir>/<name>  its annotation, in the tab-separated format
#                     get_annotated_regions_fromverdant.pl reads
# and prints the comma-separated guide list on stdout, ready to hand to the
# pipeline as its second argument.
#
# Downloads are cached by accession in <cachedir> (default ./gb_cache), so
# re-running for another target only fetches what is new.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

RANKING="${1:?usage: fetch_guides.sh <ranking.tsv> <seqdir> <filedir> [cachedir]}"
SEQDIR="${2:?}"
FILEDIR="${3:?}"
CACHE="${4:-./gb_cache}"

# gb_to_truth.pl converts a GenBank flatfile to the annotation format used here.
GB2TRUTH=""
for c in "$HERE/gb_to_truth.pl" "$HERE/test/gb_to_truth.pl"; do
    [ -f "$c" ] && GB2TRUTH="$c" && break
done
[ -n "$GB2TRUTH" ] || { echo "ERROR: cannot find gb_to_truth.pl" >&2; exit 1; }

mkdir -p "$SEQDIR" "$FILEDIR" "$CACHE"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

guides=""
while IFS=$'\t' read -r acc org fam ord jac rank; do
    [ -n "${acc:-}" ] || continue
    # A filesystem- and pipeline-safe name: the pipeline splits its guide list on
    # commas and looks the names up as filenames.
    name="$(echo "$org" | tr ' ' '_' | tr -cd '[:alnum:]_.-')"
    [ -n "$name" ] || name="$acc"
    name="${name}_${acc%%.*}"

    gb="$CACHE/${acc}.gb"
    fa="$CACHE/${acc}.fsa"
    if [ ! -s "$gb" ]; then
        curl -sS -m 180 --retry 2 --retry-delay 3 \
            "${EUTILS}?db=nuccore&id=${acc}&rettype=gbwithparts&retmode=text" -o "$gb"
        sleep 1
    fi
    if [ ! -s "$fa" ]; then
        curl -sS -m 120 --retry 2 --retry-delay 3 \
            "${EUTILS}?db=nuccore&id=${acc}&rettype=fasta&retmode=text" -o "$fa"
        sleep 1
    fi
    if [ ! -s "$gb" ] || [ ! -s "$fa" ]; then
        echo "  WARNING: download failed for $acc ($org); skipping" >&2
        rm -f "$gb" "$fa"
        continue
    fi

    if ! perl "$GB2TRUTH" "$gb" > "$FILEDIR/$name" 2>/dev/null || [ ! -s "$FILEDIR/$name" ]; then
        echo "  WARNING: no usable annotation parsed from $acc ($org); skipping" >&2
        rm -f "$FILEDIR/$name"
        continue
    fi
    cp "$fa" "$SEQDIR/$name"
    guides="${guides:+$guides,}$name"
    echo "  guide $rank: $name ($fam, jaccard $jac)" >&2
done < "$RANKING"

[ -n "$guides" ] || { echo "ERROR: no guides could be prepared" >&2; exit 1; }
echo "$guides"
