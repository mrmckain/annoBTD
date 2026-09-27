#!/bin/bash
# annoBTD: annotate a plastid genome from GenBank guide genomes chosen for it.
#
# USAGE: annobtd.sh <plastome.fasta> <prefix> [--family F --order O [--genus G]]
#            [--db curated|full|hybrid] [-n 8] [--outdir DIR] [--refdb DIR]
#            [--exclude ACC[,ACC..]] [--allow-congeners]
#
#   <plastome.fasta>  one circular plastome (several records are concatenated, with a
#                     warning); <prefix> names the output files
#   --family/--order  the target's lineage, which chooses the length expectations and
#                     the guide filters; when omitted, taken from the three most similar
#                     genomes in the database (printed, so you can check)
#   --genus           optional; enables the genus stratum of the profile
#   --db              curated: 7,912 structurally verified plastomes (default: hybrid)
#                     full:    every complete plastome in GenBank
#                     hybrid:  curated when it holds 3 same-family guides, else full
#   -n                guides per genome (default 8)
#   --exclude         accessions never used as guides (e.g. the genome's own record)
#   --outdir          default ./<prefix>_annobtd
#   --refdb           reference database directory (default ANNOBTD_REFDB, else refdb/
#                     beside this script; see fetch_refdb.sh)
#
# OUTPUT (in --outdir):
#   <prefix>_annotation.txt   gene<TAB>start<TAB>end<TAB>strand, 1-based, inclusive;
#                             exons as gene_exonN, introns as gene_intronN, spacers as
#                             geneA~geneB, plus LSC/SSC/IRA/IRB/FULL region rows
#   <prefix>_guides.tsv       the guides used and why each was chosen
#   <prefix>_post_filter.txt  every post-annotation decision (duplicates, names,
#                             strands, junctions, starts)
#   <prefix>_checks/          check_annotation.tsv and check_splice.tsv
#   work/                     the full pipeline run (BLAST, ORFs, logs)
#
# Needs: perl 5, BLAST+ (blastn, tblastx, makeblastdb on PATH or BLAST_BIN), curl for
# fetching guide records from NCBI, and bin/aragorn (built from third_party/aragorn) for
# tRNA boundaries; without ARAGORN the tRNA steps that need it are skipped.
set -uo pipefail
ANNOBTD_DIR="$(cd "$(dirname "$0")" && pwd)"
FASTA="${1:?usage: annobtd.sh <plastome.fasta> <prefix> [options]}"; PREFIX="${2:?prefix}"; shift 2
FAM=""; ORD=""; GEN="NA"; DBMODE="hybrid"; NGUIDES=8; OUTDIR=""; REFDB="${ANNOBTD_REFDB:-$ANNOBTD_DIR/refdb}"; EXCL=""; CONGENERS=""
while [ $# -gt 0 ]; do
    case "$1" in
        --family) FAM="$2"; shift 2 ;; --order) ORD="$2"; shift 2 ;; --genus) GEN="$2"; shift 2 ;;
        --db) DBMODE="$2"; shift 2 ;; -n) NGUIDES="$2"; shift 2 ;; --outdir) OUTDIR="$2"; shift 2 ;;
        --refdb) REFDB="$2"; shift 2 ;; --exclude) EXCL="$2"; shift 2 ;; --allow-congeners) CONGENERS="--allow-congeners"; shift ;;
        *) echo "annobtd.sh: unknown option $1" >&2; exit 1 ;;
    esac
done
[ -s "$FASTA" ] || { echo "annobtd.sh: cannot read $FASTA" >&2; exit 1; }
FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"
OUTDIR="${OUTDIR:-$(pwd)/${PREFIX}_annobtd}"; WORK="$OUTDIR/work"; mkdir -p "$WORK/sequenceFiles" "$WORK/files" "$OUTDIR/${PREFIX}_checks"
# reference database and tools
for f in reference_sketches.tsv gene_expect.tsv species_units.tsv species_lengths.tsv all_lengths.tsv taxonomy_consistency.tsv lineage_species_counts.tsv; do
    [ -s "$REFDB/$f" ] || { echo "annobtd.sh: $REFDB/$f is missing. Run $ANNOBTD_DIR/fetch_refdb.sh (or test/refresh_references.sh --install)." >&2; exit 1; }
done
if [ "$DBMODE" != curated ]; then [ -s "$REFDB/reference_sketches_full.tsv" ] || { echo "annobtd.sh: $REFDB/reference_sketches_full.tsv is missing (needed for --db $DBMODE)." >&2; exit 1; }; fi
if [ -n "${BLAST_BIN:-}" ]; then [ -x "$BLAST_BIN/tblastx" ] || { echo "annobtd.sh: no tblastx in BLAST_BIN=$BLAST_BIN" >&2; exit 1; }
else command -v tblastx >/dev/null && command -v blastn >/dev/null && command -v makeblastdb >/dev/null || { echo "annobtd.sh: BLAST+ (blastn, tblastx, makeblastdb) not on PATH; install it or set BLAST_BIN" >&2; exit 1; }; fi
command -v curl >/dev/null || { echo "annobtd.sh: curl is needed to fetch guide records from NCBI" >&2; exit 1; }
[ -x "${ANNOBTD_ARAGORN:-$ANNOBTD_DIR/bin/aragorn}" ] || echo "annobtd.sh: WARNING bin/aragorn not found; tRNA boundary refinement and strand checks are skipped (see third_party/aragorn/README.txt)" >&2
case "$DBMODE" in
    curated) DB="$REFDB/reference_sketches.tsv"; DB2="" ;;
    full)    DB="$REFDB/reference_sketches_full.tsv"; DB2="" ;;
    hybrid)  DB="$REFDB/reference_sketches.tsv"; DB2="$REFDB/reference_sketches_full.tsv" ;;
    *) echo "annobtd.sh: --db must be curated, full or hybrid" >&2; exit 1 ;;
esac
# lineage: given, or the majority family/order of the three nearest genomes
if [ -z "$FAM" ] || [ -z "$ORD" ]; then
    NEAR="${DB2:-$DB}"
    perl "$ANNOBTD_DIR/select_guides.pl" "$FASTA" "$NEAR" -n 3 --one-per-genus ${EXCL:+--exclude "$EXCL"} > "$WORK/nearest.tsv" 2>/dev/null
    [ -s "$WORK/nearest.tsv" ] || { echo "annobtd.sh: could not rank the genome against $NEAR" >&2; exit 1; }
    [ -z "$FAM" ] && FAM=$(cut -f3 "$WORK/nearest.tsv" | sort | uniq -c | sort -rn | awk 'NR==1{print $2}')
    [ -z "$ORD" ] && ORD=$(cut -f4 "$WORK/nearest.tsv" | sort | uniq -c | sort -rn | awk 'NR==1{print $2}')
    echo "lineage taken from the nearest genomes: family $FAM, order $ORD ($(awk -F'\t' '{printf "%s %.3f; ", $2, $5}' "$WORK/nearest.tsv"))" >&2
fi
export ANNOBTD_FAMILY="$FAM" ANNOBTD_ORDER="$ORD" ANNOBTD_GENUS="$GEN"
export ANNOBTD_PROFILE="$REFDB/gene_expect.tsv" ANNOBTD_LINEAGE_COUNTS="$REFDB/lineage_species_counts.tsv"
export ANNOBTD_SPECIES_UNITS="$REFDB/species_units.tsv" ANNOBTD_SPECIES_GENES="$REFDB/species_lengths.tsv"
export ANNOBTD_TRNA_LIBRARY="${ANNOBTD_TRNA_LIBRARY:-$ANNOBTD_DIR/trna_library.fasta}"
[ -n "$EXCL" ] && export ANNOBTD_TRNA_LIBRARY_EXCLUDE="$EXCL"
# guides
hyb=$("$ANNOBTD_DIR/choose_guides.sh" "$FASTA" "$WORK" --db "$DB" ${DB2:+--db2 "$DB2"} -n "$NGUIDES" --family "$FAM" --order "$ORD" --genus "$GEN" ${EXCL:+--exclude "$EXCL"} $CONGENERS --refdb "$REFDB") || { echo "annobtd.sh: no guides could be chosen" >&2; exit 2; }
[ -n "$hyb" ] && echo "$hyb" >&2
CACHE="$REFDB/gb_cache"; [ -d "$CACHE" ] || { CACHE="$ANNOBTD_DIR/test/gb_cache"; [ -d "$CACHE" ] || CACHE="$OUTDIR/gb_cache"; }
guides=$("$ANNOBTD_DIR/fetch_guides.sh" "$WORK/ranking.tsv" "$WORK/sequenceFiles" "$WORK/files" "$CACHE" 2>>"$WORK/guides.log")
[ -n "$guides" ] || { echo "annobtd.sh: could not fetch the guide records (see $WORK/guides.log)" >&2; exit 2; }
echo "guides: $(echo "$guides" | tr ',' '\n' | wc -l | tr -d ' ') ($(awk -F'\t' '{printf "%s ", $2}' "$WORK/ranking.tsv"| sed 's/ $//'))" >&2
cp "$FASTA" "$WORK/target.fsa"
# annotate
( cd "$WORK" && ANNOBTD_DIR="$ANNOBTD_DIR" SEQDIR=./sequenceFiles FILEDIR=./files "$ANNOBTD_DIR/run_annotation_pipeline.sh" target.fsa "$guides" "$PREFIX" ) > "$WORK/pipeline.log" 2>&1
RESULT="$WORK/${PREFIX}_VERDANT_cleaned_annotation.txt"
[ -s "$RESULT" ] || { echo "annobtd.sh: the pipeline produced no annotation; see $WORK/pipeline.log" >&2; tail -15 "$WORK/pipeline.log" >&2; exit 3; }
cp "$RESULT" "$OUTDIR/${PREFIX}_annotation.txt"
{ printf "accession\tfamily\torder\tjaccard\tcompared\tmismatches\tstatus\tspecies\trepresentative\tspecies_mismatch\tnote\n"; cat "$WORK/guide_picks.tsv"; } > "$OUTDIR/${PREFIX}_guides.tsv"
[ -s "$WORK/${PREFIX}_post_filter.txt" ] && cp "$WORK/${PREFIX}_post_filter.txt" "$OUTDIR/"
for c in check_annotation check_splice; do [ -s "$WORK/${PREFIX}_$c.tsv" ] && cp "$WORK/${PREFIX}_$c.tsv" "$OUTDIR/${PREFIX}_checks/"; done
awk -F'\t' '$1 !~ /~|^(LSC|SSC|IRA|IRB|FULL)$/ && $1 !~ /_intron/ { n=$1; sub(/_exon[0-9]+$/,"",n); t=(n ~ /^trn/)?"tRNA":(n ~ /^rrn/)?"rRNA":"CDS"; g[t][n]=1 } END { for (t in g) { c=0; for (k in g[t]) c++; printf "  %s genes: %d\n", t, c } }' "$RESULT" 2>/dev/null || true
echo "annotation -> $OUTDIR/${PREFIX}_annotation.txt" >&2
