#!/bin/bash
# Leave-one-out regression over the test set: annotate each genome using all the
# others as guide references, then score the result against its GenBank truth.
#
# USAGE: ./run_regression.sh [--only Genus_species] [--outdir DIR] [--self]
#
# --self annotates each genome using ONLY its own annotation as the guide. That
# removes reference divergence entirely: the pipeline is given the exact answer,
# so any error is the algorithm's own. It is the sharpest diagnostic here.
#
# Needs BLAST+ on PATH or BLAST_BIN set (see ../run_multiblast.sh).
#
# Reads the reference database in ../refdb (ANNOBTD_REFDB to point elsewhere;
# filled by ../fetch_refdb.sh or built by ./refresh_references.sh --install):
# gene_expect.tsv, species_units.tsv, species_lengths.tsv, lineage_species_counts.tsv,
# taxonomy_consistency.tsv, all_lengths.tsv, taxonomy_all.tsv, and the guide
# databases reference_sketches.tsv (curated) and reference_sketches_full.tsv (every
# cached GenBank plastome). Guides are chosen by ../choose_guides.sh.
# Results land in <outdir>/ : one .compare.txt summary and one .detail.tsv per genome,
# plus summary.tsv collecting the headline numbers across the whole set.

set -uo pipefail
cd "$(dirname "$0")"

ANNOBTD_DIR="$(cd .. && pwd)"
REFDB="${ANNOBTD_REFDB:-$ANNOBTD_DIR/refdb}"
DATA="$(pwd)/data"
OUTDIR="results"
ONLY=""
SELF=""
DB=""
DB2=""
EXGENUS=""
NGUIDES=3
GENUSOPT="--one-per-genus"

ORIG_ARGS="$*"

while [ $# -gt 0 ]; do
    case "$1" in
        --only)   ONLY="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --self)   SELF=1; shift ;;
        --db)     DB="$2"; shift 2 ;;
        --db2)    DB2="$2"; shift 2 ;;   # secondary database: hybrid choice per target (see below)
        --exclude-genus) EXGENUS=1; shift ;;   # benchmark: no guide from the target's genus (novel-lineage simulation)
        --nguides) NGUIDES="$2"; shift 2 ;;
        --allow-congeners) GENUSOPT=""; shift ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

[ -d "$DATA" ] || { echo "No test data. Run ./fetch_test_set.sh first." >&2; exit 1; }

GENOMES=()
for f in "$DATA"/*.truth.tsv; do
    n="$(basename "$f" .truth.tsv)"
    [ -n "$ONLY" ] && [ "$n" != "$ONLY" ] && continue
    GENOMES+=("$n")
done
[ ${#GENOMES[@]} -gt 0 ] || { echo "No genomes selected." >&2; exit 1; }

ALL=()
for f in "$DATA"/*.truth.tsv; do ALL+=("$(basename "$f" .truth.tsv)"); done

mkdir -p "$OUTDIR"

# Record what this run actually was, into the results directory.
#
# Without this a results directory cannot be reproduced or honestly compared
# against another. That is not hypothetical: results_v18 was produced with
# ANNOBTD_SCREEN_REFS=0, nothing recorded it, and a later comparison against it
# read as a 10-feature regression in the annotator when the only difference was
# that one variable. Every ANNOBTD_* knob is written with its effective value -
# including the ones left unset, so a default that changes later is still visible
# - together with the command line and a checksum of each script on the annotation
# path, so a code change between runs is detectable too.
CONFIG="$OUTDIR/run_config.txt"
{
    echo "# annoBTD regression run"
    echo "date	$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "host	$(hostname)"
    echo "command	$0 $ORIG_ARGS"
    echo "annobtd_dir	$ANNOBTD_DIR"
    echo "mode	$([ -n "$SELF" ] && echo self || echo leave-one-out)"
    echo "nguides	$NGUIDES"
    echo "db	${DB:-(none)}"
    echo "db2	${DB2:-(none)}"
    echo "exclude_genus	$([ -n "$EXGENUS" ] && echo yes || echo no)"
    echo
    echo "# environment (value, or the default that applied when unset)"
    for v in ANNOBTD_MAX_OVERLAP:0.9 ANNOBTD_MIN_ORF_NT:19 ANNOBTD_MIN_LEN_RATIO:0.6 \
             ANNOBTD_SCREEN_REFS:1 ANNOBTD_DISABLE_RANKING:unset ANNOBTD_REF_RANKING:unset \
             ANNOBTD_SCORE_K:2 ANNOBTD_RANK_K:3 ANNOBTD_SCORE_MIN:0.70 ANNOBTD_EXPECT_MIN_CONSENSUS:0.90 ANNOBTD_EXPECT_WEAK_CONSENSUS:1.01 ANNOBTD_EXPECT_MIN_N:10 ANNOBTD_EXPECT_MIN_N_GENUS:10 \
             ANNOBTD_SPLICE_THRESHOLD:0.5 ANNOBTD_SPLICE_REFINE:unset \
             ANNOBTD_KMER_K:unset ANNOBTD_SKETCH_SIZE:unset ANNOBTD_DENOVO_NODES:120000 \
             ANNOBTD_SKIP_CHECKS:unset ANNOBTD_TRACE:unset ANNOBTD_PROFILE:unset ANNOBTD_FAMILY:unset ANNOBTD_ORDER:unset ANNOBTD_GUIDE_QUALITY:unset ANNOBTD_MAX_OUTLIERS:3 ANNOBTD_GUIDE_LENGTHS:all_lengths.tsv ANNOBTD_MAX_MISMATCH:2 ANNOBTD_CONSISTENCY_EXON_LENGTH:unset ANNOBTD_CONSISTENCY_EXON_COUNT:unset ANNOBTD_NO_CONSISTENCY:unset ANNOBTD_SPECIES_UNITS:species_units.tsv ANNOBTD_MAX_SPECIES_MISMATCH:unset ANNOBTD_NO_SPECIES_COLLAPSE:unset ANNOBTD_CANDIDATES_MULT:6 ANNOBTD_LINEAGE_COUNTS:lineage_species_counts.tsv ANNOBTD_MIN_PRESENCE:0.5 ANNOBTD_PRESENCE_MODE:family ANNOBTD_SPECIES_GENES:species_lengths.tsv ANNOBTD_SKIP_POSTFILTER:unset ANNOBTD_TAXONOMY_CHECK:taxonomy_consistency.tsv ANNOBTD_NO_TAXONOMY_CHECK:unset ANNOBTD_HYBRID_MIN_FAMILY:3 ANNOBTD_SPECIES_MISMATCH_WEIGHT:0 ANNOBTD_ARAGORN:bin/aragorn ANNOBTD_TRNA_CONVENTION:trna_window_residuals.tsv ANNOBTD_NO_ARAGORN:unset ANNOBTD_TRNA_LIBRARY:trna_library.fasta ANNOBTD_NO_TRNA_LIBRARY:unset \
             ANNOBTD_DEBUG_COORDS:unset ANNOBTD_DEBUG_SPLICE:unset ; do
        name="${v%%:*}"; def="${v##*:}"
        if [ -n "${!name+x}" ]; then echo "$name	${!name}	set"
        else echo "$name	$def	default"; fi
    done
    echo
    echo "# md5 of the scripts on the annotation path"
    for f in "$ANNOBTD_DIR"/*.pl "$ANNOBTD_DIR"/*.sh; do
        [ -f "$f" ] || continue
        echo "$(basename "$f")	$(md5 -q "$f" 2>/dev/null || md5sum "$f" | cut -d" " -f1)"
    done
} > "$CONFIG"
echo "run config -> $CONFIG"

SUMMARY="$OUTDIR/summary.tsv"
printf "genome\ttruth\tpredicted\tmatched\trecall%%\tboth_exact%%\t5p_exact%%\t3p_exact%%\t5p_median\t3p_median\tmissing\tspurious\tseconds\n" > "$SUMMARY"

for target in "${GENOMES[@]}"; do
    echo "=============================================================="
    echo "  $target"
    echo "=============================================================="

    WORK="$OUTDIR/work_$target"
    rm -rf "$WORK"; mkdir -p "$WORK/sequenceFiles" "$WORK/files"

    # Third scoring signal: the mined per-gene expectation, stratified by the
    # target's own lineage (looked up by accession). Resolved before guide selection
    # because the lineage-consistency filter below needs it.
    TAXTABLE="$REFDB/taxonomy_all.tsv"
    tacc=$(grep -m1 '^>' "$DATA/$target.fsa" | sed 's/^>\([^ .]*\).*/\1/')
    read -r ANNOBTD_FAMILY ANNOBTD_ORDER <<< "$(awk -F'\t' -v a="$tacc" 'NR>1{b=$1; sub(/\.[0-9]+$/,"",b); if(b==a){print $6" "$7; exit}}' "$TAXTABLE")"
    export ANNOBTD_FAMILY="${ANNOBTD_FAMILY:-NA}" ANNOBTD_ORDER="${ANNOBTD_ORDER:-NA}"
    # the tRNA sequence library must not name a gene from the target's own record
    export ANNOBTD_TRNA_LIBRARY_EXCLUDE="$target,$tacc"
    # genus: first word of the organism name in the taxonomy table (hybrids "x Genus" allowed)
    ANNOBTD_GENUS=$(awk -F'\t' -v a="$tacc" 'NR>1{b=$1; sub(/\.[0-9]+$/,"",b); if(b==a){n=split($2,w," "); g=w[1]; if(g=="x" && n>1) g=w[2]; if(g ~ /^[A-Z][a-z]+$/) print g; exit}}' "$TAXTABLE"); export ANNOBTD_GENUS="${ANNOBTD_GENUS:-NA}"
    PROFILE_FILE="${ANNOBTD_PROFILE_FILE:-$REFDB/gene_expect.tsv}"
    if [ -z "${ANNOBTD_NO_PROFILE:-}" ] && [ -s "$PROFILE_FILE" ]; then export ANNOBTD_PROFILE="$PROFILE_FILE"; else unset ANNOBTD_PROFILE; fi
    COUNTS_FILE="${ANNOBTD_LINEAGE_COUNTS_FILE:-$REFDB/lineage_species_counts.tsv}"
    GENES_FILE="${ANNOBTD_SPECIES_GENES_FILE:-$REFDB/species_lengths.tsv}"
    if [ -s "$GENES_FILE" ]; then export ANNOBTD_SPECIES_GENES="$GENES_FILE"; else unset ANNOBTD_SPECIES_GENES; fi
    export ANNOBTD_SPECIES_UNITS="${ANNOBTD_SPECIES_UNITS:-$REFDB/species_units.tsv}"
    if [ -s "$COUNTS_FILE" ]; then export ANNOBTD_LINEAGE_COUNTS="$COUNTS_FILE"; else unset ANNOBTD_LINEAGE_COUNTS; fi
    printf "target\t%s\tANNOBTD_PROFILE=%s\tANNOBTD_FAMILY=%s\tANNOBTD_ORDER=%s\tANNOBTD_GENUS=%s\n" "$target" "${ANNOBTD_PROFILE:-unset}" "$ANNOBTD_FAMILY" "$ANNOBTD_ORDER" "$ANNOBTD_GENUS" >> "$CONFIG"

    guides=""
    if [ -n "$DB" ]; then
        # Guides come from the reference database, chosen by k-mer similarity to
        # this genome rather than from the handful of other test genomes. The
        # target's own accession is excluded so leave-one-out stays honest.
        # The target's own record (and the other records of its genome unit, e.g. a
        # RefSeq copy and its INSDC source) are excluded so leave-one-out stays honest.
        # choose_guides.sh applies the lineage-consistency, species-collapse and
        # taxonomy-check filters and the hybrid rule; see its header.
        tacc_full=$(grep -m1 '^>' "$DATA/$target.fsa" | sed 's/^>\([^ ]*\).*/\1/')
        hyb=$("$ANNOBTD_DIR/choose_guides.sh" "$DATA/$target.fsa" "$WORK" --db "$DB" ${DB2:+--db2 "$DB2"} -n "$NGUIDES" \
                --family "$ANNOBTD_FAMILY" --order "$ANNOBTD_ORDER" --genus "$ANNOBTD_GENUS" --exclude "$tacc_full" \
                ${EXGENUS:+--exclude-genus} $( [ -z "$GENUSOPT" ] && echo --allow-congeners ) --refdb "$REFDB") \
            || { echo "  FAILED: no guides selected from $DB" >&2; continue; }
        [ -n "$hyb" ] && printf "hybrid\t%s\t%s\n" "$target" "${hyb#hybrid	}" >> "$CONFIG"
        guides=$("$ANNOBTD_DIR/fetch_guides.sh" "$WORK/ranking.tsv" \
                 "$WORK/sequenceFiles" "$WORK/files" "$(pwd)/gb_cache" 2>>"$WORK/guides.log")
        if [ -z "$guides" ]; then
            echo "  FAILED: no guides selected from $DB" >&2
            continue
        fi
    elif [ -n "$SELF" ]; then
        # Self-annotation: the genome is its own and only guide. The pipeline is
        # handed the exact annotation it should reproduce, so anything less than a
        # perfect result is an algorithmic defect, not reference divergence.
        cp "$DATA/$target.fsa"       "$WORK/sequenceFiles/$target"
        cp "$DATA/$target.truth.tsv" "$WORK/files/$target"
        guides="$target"
    else
        # Every other genome in the set acts as a guide reference.
        for g in "${ALL[@]}"; do
            [ "$g" = "$target" ] && continue
            cp "$DATA/$g.fsa"       "$WORK/sequenceFiles/$g"
            cp "$DATA/$g.truth.tsv" "$WORK/files/$g"
            guides="${guides:+$guides,}$g"
        done
    fi
    cp "$DATA/$target.fsa" "$WORK/target.fsa"

    start=$(date +%s)
    ( cd "$WORK" && \
      ANNOBTD_DIR="$ANNOBTD_DIR" SEQDIR=./sequenceFiles FILEDIR=./files \
      "$ANNOBTD_DIR/run_annotation_pipeline.sh" target.fsa "$guides" "$target" \
    ) > "$WORK/pipeline.log" 2>&1
    rc=$?
    elapsed=$(( $(date +%s) - start ))

    RESULT="$WORK/${target}_VERDANT_cleaned_annotation.txt"
    if [ $rc -ne 0 ] || [ ! -s "$RESULT" ]; then
        echo "  FAILED (exit $rc, ${elapsed}s) - see $WORK/pipeline.log"
        tail -15 "$WORK/pipeline.log" | sed 's/^/    /'
        printf "%s\tFAILED\t\t\t\t\t\t\t\t\t\t\t%d\n" "$target" "$elapsed" >> "$SUMMARY"
        continue
    fi

    perl compare_annotations.pl "$DATA/$target.truth.tsv" "$RESULT" \
        --label "$target" --detail "$OUTDIR/$target.detail.tsv" \
        | tee "$OUTDIR/$target.compare.txt"
    echo "  (${elapsed}s)"

    perl summarize.pl "$OUTDIR/$target.compare.txt" "$SUMMARY" "$target" "$elapsed"

    # BLAST indexes are pure derived data and dominate the on-disk size.
    rm -f "$WORK"/*.n?? "$WORK"/*.njs 2>/dev/null

    echo
done

echo "=============================================================="
echo "  SUMMARY  ($SUMMARY)"
echo "=============================================================="
column -t -s $'\t' "$SUMMARY"
