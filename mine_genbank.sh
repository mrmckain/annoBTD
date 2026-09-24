#!/bin/bash
# Fetch the full GenBank flat file for every accession in a taxonomy table into a
# cache directory, skipping any already present. Resumable: stop and restart at
# will. ~0.5 records/s against NCBI's rate limit.
#
# USAGE: mine_genbank.sh <taxonomy.tsv> <cache_dir>
set -u
TAX="$1"; CACHE="$2"; mkdir -p "$CACHE"
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
total=$(awk 'NR>1' "$TAX" | wc -l | tr -d ' ')
have=$(ls "$CACHE"/*.gb 2>/dev/null | wc -l | tr -d ' ')
echo "resuming: $have of $total already cached" >&2
n=0
awk -F'\t' 'NR>1{print $1}' "$TAX" | while read -r acc; do
    [ -s "$CACHE/$acc.gb" ] && continue
    curl -sS -m 180 --retry 3 "${E}?db=nuccore&id=${acc}&rettype=gbwithparts&retmode=text" -o "$CACHE/$acc.gb.part" \
      && grep -q '^ORIGIN' "$CACHE/$acc.gb.part" && mv "$CACHE/$acc.gb.part" "$CACHE/$acc.gb" \
      || { rm -f "$CACHE/$acc.gb.part"; echo "failed: $acc" >&2; }
    n=$((n+1)); [ $((n % 100)) -eq 0 ] && echo "fetched $n" >&2
    sleep 0.4
done
echo "done" >&2
