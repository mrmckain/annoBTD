#!/bin/bash
# Fetch (or build) the annoBTD reference database bundle.
#
# USAGE: fetch_refdb.sh [--version V] [--url BASE_URL] [--dir DIR]      download + install
#        fetch_refdb.sh --build V [--dir DIR] [--out DIR]                 pack DIR into a bundle
#
# The reference database is too large for git (1.5 GB: two MinHash sketch databases,
# the taxonomy of every complete plastome in GenBank, and the profile tables derived
# from 67,461 records). It is published as one gzip tarball per release, with a
# manifest and a SHA-256 checksum:
#   annoBTD_refdb_<V>.tar.gz          the files, relative to the refdb directory
#   annoBTD_refdb_<V>.tar.gz.sha256   checksum of the tarball
# Default source: the GitHub release "refdb-<V>" of mrmckain/annoBTD. --url points at
# another directory of the same two files (a mirror, an institutional server, a local
# path as file:///...). The bundle is unpacked into --dir (default ANNOBTD_REFDB, else
# refdb/ beside this script), replacing files of the same name; the manifest is
# verified after unpacking.
#
# --build packs the current --dir into a bundle in --out (default .), writing the
# manifest (file, bytes, rows, sha256) and VERSION (date, source snapshot) into it.
# Publishing is then: create a release tagged refdb-<V> and attach both files.
set -uo pipefail
ANNOBTD_DIR="$(cd "$(dirname "$0")" && pwd)"
VERSION="2026-09"; URL=""; DIR="${ANNOBTD_REFDB:-$ANNOBTD_DIR/refdb}"; BUILD=""; OUT="."
while [ $# -gt 0 ]; do
    case "$1" in
        --version) VERSION="$2"; shift 2 ;; --url) URL="$2"; shift 2 ;; --dir) DIR="$2"; shift 2 ;;
        --build) BUILD=1; VERSION="$2"; shift 2 ;; --out) OUT="$2"; shift 2 ;;
        *) echo "fetch_refdb.sh: unknown option $1" >&2; exit 1 ;;
    esac
done
NAME="annoBTD_refdb_$VERSION.tar.gz"
sha() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
# the files a run needs (record_meta.tsv and species_consensus.tsv are refresh-only audit tables and stay out)
FILES="reference_sketches.tsv reference_sketches_full.tsv taxonomy_all.tsv gene_expect.tsv gene_profile.tsv species_units.tsv species_lengths.tsv all_lengths.tsv taxonomy_consistency.tsv lineage_species_counts.tsv lineage_assignments.tsv"

if [ -n "$BUILD" ]; then
    cd "$DIR" || exit 1
    for f in $FILES; do [ -s "$f" ] || { echo "fetch_refdb.sh: $DIR/$f missing" >&2; exit 1; }; done
    { printf "file\tbytes\trows\tsha256\n"; for f in $FILES; do printf "%s\t%s\t%s\t%s\n" "$f" "$(stat -f %z "$f" 2>/dev/null || stat -c %s "$f")" "$(($(wc -l < "$f")))" "$(sha "$f")"; done; } > MANIFEST.tsv
    { echo "version	$VERSION"; echo "built	$(date -u '+%Y-%m-%dT%H:%M:%SZ')"; echo "records	$(($(wc -l < taxonomy_all.tsv) - 1)) complete plastomes in taxonomy_all.tsv"; echo "sketches	$(wc -l < reference_sketches_full.tsv | tr -d ' ') in reference_sketches_full.tsv, $(wc -l < reference_sketches.tsv | tr -d ' ') curated in reference_sketches.tsv"; echo "species	$(awk -F'\t' 'NR>1 && $8==1' species_units.tsv | wc -l | tr -d ' ') species representatives in species_units.tsv"; } > VERSION
    mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
    echo "packing $DIR -> $OUT/$NAME" >&2
    tar -czf "$OUT/$NAME" MANIFEST.tsv VERSION $FILES || exit 1
    ( cd "$OUT" && sha "$NAME" > "$NAME.sha256" && ls -la "$NAME" "$NAME.sha256" >&2 )
    echo "publish: create a GitHub release tagged refdb-$VERSION on mrmckain/annoBTD and attach $NAME and $NAME.sha256" >&2
    exit 0
fi

URL="${URL:-https://github.com/mrmckain/annoBTD/releases/download/refdb-$VERSION}"
mkdir -p "$DIR"; TMP="$DIR/.download"; mkdir -p "$TMP"
echo "fetching $URL/$NAME" >&2
curl -fL --retry 3 -o "$TMP/$NAME.sha256" "$URL/$NAME.sha256" || { echo "fetch_refdb.sh: cannot fetch the checksum from $URL" >&2; exit 1; }
curl -fL --retry 3 -C - -o "$TMP/$NAME" "$URL/$NAME" || { echo "fetch_refdb.sh: download failed" >&2; exit 1; }
want=$(cut -d' ' -f1 "$TMP/$NAME.sha256"); have=$(sha "$TMP/$NAME")
[ "$want" = "$have" ] || { echo "fetch_refdb.sh: checksum mismatch (expected $want, got $have); not installing" >&2; exit 1; }
echo "checksum ok; unpacking into $DIR" >&2
tar -xzf "$TMP/$NAME" -C "$DIR" || exit 1
rm -rf "$TMP"
bad=0; while IFS=$'\t' read -r f bytes rows sum; do [ "$f" = file ] && continue; [ "$(sha "$DIR/$f")" = "$sum" ] || { echo "  $f: checksum differs from the manifest" >&2; bad=1; }; done < "$DIR/MANIFEST.tsv"
[ $bad = 0 ] && { echo "reference database $VERSION installed in $DIR:" >&2; sed 's/^/  /' "$DIR/VERSION" >&2; } || exit 1
