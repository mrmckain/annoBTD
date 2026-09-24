#!/bin/bash
# USAGE: run_annotation_pipeline.sh <plastome.fsa> <guide1,guide2,...> <species_prefix>
#
# Guide names are looked up as $SEQDIR/<name> (sequence) and $FILEDIR/<name> (annotation).
#
# ANNOBTD_DIR - where these scripts live      (default: annoBTD_multiref)
# SEQDIR      - guide sequences               (default: ./sequenceFiles)
# FILEDIR     - guide annotations             (default: ./files)
# See run_multiblast.sh for BLAST_BIN / BLAST_THREADS / BLAST_EVALUE.

set -euo pipefail

ANNOBTD_DIR="${ANNOBTD_DIR:-annoBTD_multiref}"
SEQDIR="${SEQDIR:-./sequenceFiles}"
FILEDIR="${FILEDIR:-./files}"

PLASTOME="$1"
GUIDES="$2"
SP="$3"

[ -s "$PLASTOME" ] || { echo "ERROR: no plastome at $PLASTOME" >&2; exit 1; }

# get_annotated_regions_fromverdant.pl appends, so stale output from an earlier run
# would silently contaminate this one.
rm -f "${SP}_annotated_regions_fromverdant_genes.fsa" \
      "${SP}_annotated_regions_fromverdant_trnas.fsa" \
      "${SP}_annotated_regions_fromverdant_rrnas.fsa"

# Rank the guides by k-mer similarity to this genome, so downstream reference
# choice can prefer a close relative. Exon structure is lineage specific, and a
# guide picked for convenience transfers the wrong splicing pattern.
if [ -n "${ANNOBTD_DISABLE_RANKING:-}" ]; then
    echo "reference ranking disabled" >&2
    export ANNOBTD_REF_RANKING=""
else
RANKING="${SP}_reference_ranking.txt"
: > "$RANKING"
IDX=0
for file in $(echo "$GUIDES" | tr "," " "); do
    IDX=$((IDX + 1))
    echo -e "${IDX}\t${file}"
done > "${RANKING}.idx"
if perl "$ANNOBTD_DIR/rank_references_by_kmer.pl" "$PLASTOME" \
        $(awk -F'\t' -v d="$SEQDIR" '{print d"/"$2}' "${RANKING}.idx") \
        > "${RANKING}.raw" 2>/dev/null; then
    awk -F'\t' 'NR==FNR{idx[$2]=$1; next} {print idx[$1]"\t"$1"\t"$4}' \
        "${RANKING}.idx" "${RANKING}.raw" > "$RANKING"
    echo "reference ranking (closest first):" >&2
    sort -t$'\t' -k3 -n "$RANKING" | head -3 | sed 's/^/  /' >&2

    # Splice-site refinement helps only when the guide is distant enough that its
    # alignment cannot be trusted for the exon junction. With a close relative it
    # moves correct boundaries - a real intron lacking a GT donor scores 0 and any
    # neighbour looks better - so gate it on the closest guide's k-mer containment.
    if [ -z "${ANNOBTD_SPLICE_REFINE:-}" ]; then
        TOPC=$(sort -t$'\t' -k3 -rn "${RANKING}.raw" 2>/dev/null | head -1 | cut -f3)
        if [ -n "$TOPC" ] && awk -v c="$TOPC" -v t="${ANNOBTD_SPLICE_THRESHOLD:-0.5}" \
               'BEGIN{exit !(c < t)}'; then
            export ANNOBTD_SPLICE_REFINE=1
            echo "closest guide containment $TOPC < ${ANNOBTD_SPLICE_THRESHOLD:-0.5}: splice refinement ON" >&2
        else
            echo "closest guide containment $TOPC: splice refinement off" >&2
        fi
    fi
else
    echo "WARNING: k-mer ranking failed; reference choice falls back to score only" >&2
fi
export ANNOBTD_REF_RANKING="$RANKING"
fi

perl "$ANNOBTD_DIR/annotate_plastome.pl" "$PLASTOME" "$SP" "${ANNOBTD_MAX_OVERLAP:-0.9}" "${ANNOBTD_MIN_ORF_NT:-19}" "${ANNOBTD_MIN_LEN_RATIO:-0.6}"

COUNT=0
for file in $(echo "$GUIDES" | tr "," " "); do
    COUNT=$((COUNT + 1))
    perl "$ANNOBTD_DIR/get_annotated_regions_fromverdant.pl" \
        "$SEQDIR/$file" "$FILEDIR/$file" "$SP" "$COUNT"
done

# Screen the guide references for bad annotations before they are used. GenBank
# holds real errors, and annoBTD transfers whatever the winning reference says, so
# a reference whose own CDS does not translate is a direct route to a bad call.
# ANNOBTD_SCREEN_REFS=0 disables it.
if [ "${ANNOBTD_SCREEN_REFS:-1}" != "0" ]; then
    GF="${SP}_annotated_regions_fromverdant_genes.fsa"
    # ANNOBTD_SCREEN_COMPARATIVE=1 also enables the peer-comparison checks (length
    # and 5-mer outliers). They are off by default because they cost recall against
    # a small set of curated guides, but a large set of RAW GenBank guides is the
    # case they were designed for - see check_references.pl.
    COMPOPT=""
    [ -n "${ANNOBTD_SCREEN_COMPARATIVE:-}" ] && COMPOPT="--comparative"
    if perl "$ANNOBTD_DIR/check_references.pl" "$GF" $COMPOPT \
            --tsv "${SP}_reference_problems.tsv" --drop "${GF}.screened" 2>/dev/null; then
        if [ -s "${GF}.screened" ]; then
            mv "${GF}.screened" "$GF"
        fi
    fi
fi

"$ANNOBTD_DIR/run_multiblast.sh" \
    "${SP}_annotated_regions_fromverdant_genes.fsa" \
    "${SP}_annotated_regions_fromverdant_trnas.fsa" \
    "${SP}_annotated_regions_fromverdant_rrnas.fsa" \
    "${SP}_orffinder_seqs.fsa" \
    "$PLASTOME"

perl "$ANNOBTD_DIR/identify_best_ref_for_orf_withscoring_v0.7.pl" \
    "${SP}_orffinder_seqs.fsa_trnas.blastn" \
    "${SP}_annotated_regions_fromverdant_trnas.fsa" \
    "$PLASTOME" tRNA "$SP" "${SP}_orffinder_coordinates.txt"

# Was hardcoded to Schizachyrium.fsa, which fed the wrong sequence to every
# non-grass genome. The query for the rRNA blastn is the plastome being annotated.
perl "$ANNOBTD_DIR/identify_best_ref_for_orf_withscoring_v0.7.pl" \
    "${SP}_orffinder_seqs.fsa_rrnas.blastn" \
    "${SP}_annotated_regions_fromverdant_rrnas.fsa" \
    "$PLASTOME" rRNA "$SP" "${SP}_orffinder_coordinates.txt"

perl "$ANNOBTD_DIR/identify_best_ref_for_orf_withscoring_v0.7.pl" \
    "${SP}_orffinder_seqs.fsa_genes.tblastx" \
    "${SP}_annotated_regions_fromverdant_genes.fsa" \
    "${SP}_orffinder_seqs.fsa" protein "$SP" "${SP}_orffinder_coordinates.txt"

perl "$ANNOBTD_DIR/match_orfs_to_blast_v2.2.pl" \
    "${SP}_annotated_regions_fromverdant_genes.fsa" \
    "${SP}_orffinder_seqs.fsa" \
    "${SP}_best_orfs_for_refs_SCORE.txt" \
    "${SP}_orffinder_coordinates.txt" \
    "$PLASTOME" \
    "${SP}_annotated_regions_fromverdant_trnas.fsa" \
    "${SP}_annotated_regions_fromverdant_rrnas.fsa" \
    "${SP}_orffinder_seqs.fsa_trnas.blastn" \
    "${SP}_orffinder_seqs.fsa_rrnas.blastn" \
    "$SP" \
    "${SP}_best_tRNA_ref_SCORE.txt" \
    "${SP}_best_rRNA_ref_SCORE.txt"

# ---------------------------------------------------------------------------
# Post-annotation checks. Neither is required to produce the annotation; both
# report on it. Set ANNOBTD_SKIP_CHECKS=1 to omit them.
#
# check_annotation.pl translates what was called - frames, starts, stops, length.
# check_splice.pl reads the intron junctions, which is the only thing that catches
# a boundary drawn a few nucleotides off when the displaced version happens to
# encode the same protein. That case passes every protein-level test, and it is
# the one that gets published and then propagated into everything annotated from
# it. Exit status is not checked: these are advisory.
ANNOTATION="${SP}_VERDANT_cleaned_annotation.txt"
# Post-annotation filter (post_filter_annotation.pl): collapse duplicate tRNA/rRNA
# calls that several guides transferred under different names or shifted windows
# (coordinates decided by acceptor-stem pairing), and drop CDS genes annotated in
# fewer than ANNOBTD_MIN_PRESENCE (0.5) of the target lineage's species (needs
# ANNOBTD_PROFILE, ANNOBTD_LINEAGE_COUNTS, ANNOBTD_FAMILY/ORDER). ANNOBTD_SKIP_POSTFILTER
# disables it. tRNA names are settled by sequence against trna_library.fasta
# (build_trna_library.pl; ANNOBTD_TRNA_LIBRARY to point elsewhere, ANNOBTD_NO_TRNA_LIBRARY
# to skip, ANNOBTD_TRNA_LIBRARY_EXCLUDE=name,acc,... to leave the genome's own entries out).
if [ -z "${ANNOBTD_SKIP_POSTFILTER:-}" ] && [ -s "$ANNOTATION" ]; then
    perl "$ANNOBTD_DIR/post_filter_annotation.pl" "$PLASTOME" "$ANNOTATION" \
        --guides-dir "${FILEDIR:-./files}" \
        ${ANNOBTD_PROFILE:+--profile "$ANNOBTD_PROFILE"} ${ANNOBTD_LINEAGE_COUNTS:+--lineage-counts "$ANNOBTD_LINEAGE_COUNTS"} \
        --family "${ANNOBTD_FAMILY:-NA}" --order "${ANNOBTD_ORDER:-NA}" --min-presence "${ANNOBTD_MIN_PRESENCE:-0.5}" \
        $( [ -z "${ANNOBTD_NO_ARAGORN:-}" ] && [ -x "${ANNOBTD_ARAGORN:-$ANNOBTD_DIR/bin/aragorn}" ] && printf -- "--aragorn %s --trna-convention %s" "${ANNOBTD_ARAGORN:-$ANNOBTD_DIR/bin/aragorn}" "${ANNOBTD_TRNA_CONVENTION:-$ANNOBTD_DIR/trna_window_residuals.tsv}" ) \
        ${ANNOBTD_SPECIES_UNITS:+--species-units "$ANNOBTD_SPECIES_UNITS"} ${ANNOBTD_SPECIES_GENES:+--species-genes "$ANNOBTD_SPECIES_GENES"} --presence-mode "${ANNOBTD_PRESENCE_MODE:-family}" \
        $( [ -z "${ANNOBTD_NO_TRNA_LIBRARY:-}" ] && [ -s "${ANNOBTD_TRNA_LIBRARY:-$ANNOBTD_DIR/trna_library.fasta}" ] && printf -- "--trna-library %s" "${ANNOBTD_TRNA_LIBRARY:-$ANNOBTD_DIR/trna_library.fasta}" ) ${ANNOBTD_TRNA_LIBRARY_EXCLUDE:+--trna-library-exclude "$ANNOBTD_TRNA_LIBRARY_EXCLUDE"} ${BLAST_BIN:+--blast-bin "$BLAST_BIN"} \
        --log "${SP}_post_filter.txt" || true
fi
if [ -z "${ANNOBTD_SKIP_CHECKS:-}" ] && [ -s "$ANNOTATION" ]; then
    perl "$ANNOBTD_DIR/check_annotation.pl" "$PLASTOME" "$ANNOTATION" \
        --refs "${SP}_annotated_regions_fromverdant_genes.fsa" \
        --tsv "${SP}_check_annotation.tsv" || true
    perl "$ANNOBTD_DIR/check_splice.pl" "$PLASTOME" "$ANNOTATION" \
        --tsv "${SP}_check_splice.tsv" --fix "${SP}_splice_fixes.txt" || true
fi
