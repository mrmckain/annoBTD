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
for c in "$HERE/../test/gb_to_truth.pl" "$HERE/gb_to_truth.pl"; do
    [ -f "$c" ] && GB2TRUTH="$c" && break
done
[ -n "$GB2TRUTH" ] || { echo "ERROR: cannot find gb_to_truth.pl" >&2; exit 1; }

mkdir -p "$SEQDIR" "$FILEDIR" "$CACHE"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# NCBI allows 3 requests/s without a key, 10 with one; batches of jobs sharing a cache can hit the limit
KEY="${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}"

# A record is usable only when complete: GenBank from LOCUS to the closing //, FASTA from a > header.
# A cache file that fails this (an error page, or a download another job was killed in) is fetched again.
gb_ok() { [ -s "$1" ] && head -c 5 "$1" | grep -q '^LOCUS' && tail -n 3 "$1" | grep -q '^//'; }
fa_ok() { [ -s "$1" ] && head -c 1 "$1" | grep -q '^>'; }
# fetch <url> <check> <cache path>: download to a temporary file and move it into the cache only once it
# passes <check>, so jobs sharing the cache never read a partial record; prints where the record is.
# NCBI answers bursts with 429, so a failed try waits a random, growing interval (up to ~3 min in all)
# and looks in the cache again first: another job may have fetched the record meanwhile. An unwritable
# cache keeps the record in the run's scratch space instead. Never fails the script.
fetch() {
    local url="$1" check="$2" dest="$3" tmp try code scratch="$SEQDIR/.${3##*/}"
    for try in 1 2 3 4 5 6; do
        "$check" "$dest" && { echo "$dest"; return 0; }
        tmp=$(mktemp "$dest.part.XXXXXX" 2>/dev/null) || tmp=$(mktemp "$SEQDIR/.part.XXXXXX") || return 1
        code=$(curl -sS -m 180 -w '%{http_code}' -o "$tmp" "$url$KEY" 2>/dev/null) || code=000   # timeouts and cut connections count as transient
        if [ "$code" = 200 ] && "$check" "$tmp"; then
            chmod a+r "$tmp" 2>/dev/null
            if mv -f "$tmp" "$dest" 2>/dev/null; then echo "$dest"; else mv -f "$tmp" "$scratch" && echo "$scratch"; fi
            sleep 1; return 0
        fi
        rm -f "$tmp"
        case "$code" in 000|408|429|5??) sleep $(( try * 5 + RANDOM % (try * 10) )) ;; *) return 1 ;; esac   # only transient failures are retried
    done
    return 1
}

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
    # the cache may be gzipped: expand into a scratch file beside the guides, not into the cache
    if ! gb_ok "$gb" && [ -s "$gb.gz" ]; then gzip -dc "$CACHE/${acc}.gb.gz" > "$SEQDIR/.${acc}.gb" && gb="$SEQDIR/.${acc}.gb"; fi
    gb_ok "$gb" || gb=$(fetch "${EUTILS}?db=nuccore&id=${acc}&rettype=gbwithparts&retmode=text" gb_ok "$CACHE/${acc}.gb") || gb=""
    fa_ok "$fa" || fa=$(fetch "${EUTILS}?db=nuccore&id=${acc}&rettype=fasta&retmode=text" fa_ok "$CACHE/${acc}.fsa") || fa=""
    if [ -z "$gb" ] || [ -z "$fa" ]; then
        echo "  WARNING: download failed for $acc ($org); skipping" >&2
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

rm -f "$SEQDIR"/.*.gb "$SEQDIR"/.*.fsa "$SEQDIR"/.part.*
# guide genes that GenBank names wrongly (e.g. psbJ annotated as a second psaJ) take the name the other guides give them
perl "$HERE/relabel_guides.pl" "$FILEDIR" "$SEQDIR" >&2 || echo "  WARNING: relabel_guides.pl failed; guide names not checked" >&2
[ -n "$guides" ] || { echo "ERROR: no guides could be prepared" >&2; exit 1; }
echo "$guides"
