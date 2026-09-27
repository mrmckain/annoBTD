#!/bin/bash
# Download the regression test set from NCBI: GenBank flatfile + FASTA per genome,
# then build a ground-truth table from each flatfile.
#
# USAGE: ./fetch_test_set.sh [outdir]      (default: ./data)
#
# The set is deliberately spread across land plants so a change that helps one
# lineage and hurts another shows up. It carries three grasses (Oryza, Sorghum,
# Zea) because exon structure is lineage specific - grasses differ from tobacco in
# the exon counts of several genes - and a single grass cannot show whether
# reference selection is picking a relative or a stranger. Pinus has a drastically
# reduced IR and Marchantia is a bryophyte, so both stress IR detection and the
# guide scoring.
#
# Accessions are checked against their GenBank DEFINITION line after download:
# NC_008325 is Daucus carota, not a grass, and was mislabelled here once already.

set -euo pipefail
cd "$(dirname "$0")"
OUT="${1:-data}"
mkdir -p "$OUT"

EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# accession<TAB>name
GENOMES=$(cat <<'EOF'
NC_001879.2	Nicotiana_tabacum
NC_000932.1	Arabidopsis_thaliana
NC_002202.1	Spinacia_oleracea
NC_008325.1	Daucus_carota
NC_001320.1	Oryza_sativa
NC_008602.1	Sorghum_bicolor
NC_001666.2	Zea_mays
NC_001631.1	Pinus_thunbergii
NC_001319.1	Marchantia_paleacea
EOF
)

fetch () {
    local acc="$1" name="$2" type="$3" ext="$4"
    local dest="$OUT/${name}.${ext}"
    if [ -s "$dest" ]; then
        echo "  have  $dest"
        return
    fi
    echo "  fetch $dest"
    curl -sS -m 180 --retry 3 --retry-delay 3 \
        "${EUTILS}?db=nuccore&id=${acc}&rettype=${type}&retmode=text" -o "$dest"
    if [ ! -s "$dest" ]; then
        echo "ERROR: empty download for $acc ($type)" >&2
        rm -f "$dest"
        exit 1
    fi
    sleep 1   # stay under the NCBI rate limit
}

echo "$GENOMES" | while IFS=$'\t' read -r acc name; do
    [ -n "$acc" ] || continue
    echo "$name ($acc)"
    fetch "$acc" "$name" gbwithparts gb
    fetch "$acc" "$name" fasta       fsa

    # Guard against mislabelled accessions: the name must appear in the record.
    genus="${name%%_*}"
    if ! grep -m1 '^DEFINITION' "$OUT/${name}.gb" | grep -qi "$genus"; then
        echo "ERROR: $OUT/${name}.gb is not $name --" >&2
        grep -m1 '^DEFINITION' "$OUT/${name}.gb" >&2
        exit 1
    fi
    perl gb_to_truth.pl "$OUT/${name}.gb" > "$OUT/${name}.truth.tsv"
done

echo
echo "Test set ready in $OUT/:"
for f in "$OUT"/*.truth.tsv; do
    printf "  %-40s %4d features\n" "$(basename "$f")" "$(wc -l < "$f")"
done
