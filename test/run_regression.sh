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
# Reads the derived reference tables in this directory (built and installed by
# ./refresh_references.sh): gene_expect.tsv, species_units.tsv, species_lengths.tsv,
# lineage_species_counts.tsv, taxonomy_consistency.tsv, all_lengths.tsv, and
# ../../taxonomy_all.tsv; guide databases are ../reference_sketches.tsv (curated)
# and ../reference_sketches_full.tsv (every cached GenBank plastome).
# Results land in <outdir>/ : one .compare.txt summary and one .detail.tsv per genome,
# plus summary.tsv collecting the headline numbers across the whole set.

set -uo pipefail
cd "$(dirname "$0")"

ANNOBTD_DIR="$(cd .. && pwd)"
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
             ANNOBTD_SKIP_CHECKS:unset ANNOBTD_TRACE:unset ANNOBTD_PROFILE:unset ANNOBTD_FAMILY:unset ANNOBTD_ORDER:unset ANNOBTD_GUIDE_QUALITY:unset ANNOBTD_MAX_OUTLIERS:3 ANNOBTD_GUIDE_LENGTHS:all_lengths.tsv ANNOBTD_MAX_MISMATCH:2 ANNOBTD_CONSISTENCY_EXON_LENGTH:unset ANNOBTD_CONSISTENCY_EXON_COUNT:unset ANNOBTD_NO_CONSISTENCY:unset ANNOBTD_SPECIES_UNITS:species_units.tsv ANNOBTD_MAX_SPECIES_MISMATCH:unset ANNOBTD_NO_SPECIES_COLLAPSE:unset ANNOBTD_CANDIDATES_MULT:6 ANNOBTD_LINEAGE_COUNTS:lineage_species_counts.tsv ANNOBTD_MIN_PRESENCE:0.5 ANNOBTD_PRESENCE_MODE:family ANNOBTD_SPECIES_GENES:species_lengths.tsv ANNOBTD_SKIP_POSTFILTER:unset ANNOBTD_TAXONOMY_CHECK:taxonomy_consistency.tsv ANNOBTD_NO_TAXONOMY_CHECK:unset ANNOBTD_HYBRID_MIN_FAMILY:3 ANNOBTD_SPECIES_MISMATCH_WEIGHT:0 ANNOBTD_ARAGORN:bin/aragorn ANNOBTD_TRNA_CONVENTION:trna_window_residuals.tsv ANNOBTD_NO_ARAGORN:unset \
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
    TAXTABLE="$(cd .. && cd .. && pwd)/taxonomy_all.tsv"
    tacc=$(grep -m1 '^>' "$DATA/$target.fsa" | sed 's/^>\([^ .]*\).*/\1/')
    read -r ANNOBTD_FAMILY ANNOBTD_ORDER <<< "$(awk -F'\t' -v a="$tacc" 'NR>1{b=$1; sub(/\.[0-9]+$/,"",b); if(b==a){print $6" "$7; exit}}' "$TAXTABLE")"
    export ANNOBTD_FAMILY="${ANNOBTD_FAMILY:-NA}" ANNOBTD_ORDER="${ANNOBTD_ORDER:-NA}"
    # genus: first word of the organism name in the taxonomy table (hybrids "x Genus" allowed)
    ANNOBTD_GENUS=$(awk -F'\t' -v a="$tacc" 'NR>1{b=$1; sub(/\.[0-9]+$/,"",b); if(b==a){n=split($2,w," "); g=w[1]; if(g=="x" && n>1) g=w[2]; if(g ~ /^[A-Z][a-z]+$/) print g; exit}}' "$TAXTABLE"); export ANNOBTD_GENUS="${ANNOBTD_GENUS:-NA}"
    PROFILE_FILE="${ANNOBTD_PROFILE_FILE:-$(pwd)/gene_expect.tsv}"
    if [ -z "${ANNOBTD_NO_PROFILE:-}" ] && [ -s "$PROFILE_FILE" ]; then export ANNOBTD_PROFILE="$PROFILE_FILE"; else unset ANNOBTD_PROFILE; fi
    COUNTS_FILE="${ANNOBTD_LINEAGE_COUNTS_FILE:-$(pwd)/lineage_species_counts.tsv}"
    GENES_FILE="${ANNOBTD_SPECIES_GENES_FILE:-$(pwd)/species_lengths.tsv}"
    if [ -s "$GENES_FILE" ]; then export ANNOBTD_SPECIES_GENES="$GENES_FILE"; else unset ANNOBTD_SPECIES_GENES; fi
    export ANNOBTD_SPECIES_UNITS="${ANNOBTD_SPECIES_UNITS:-$(pwd)/species_units.tsv}"
    if [ -s "$COUNTS_FILE" ]; then export ANNOBTD_LINEAGE_COUNTS="$COUNTS_FILE"; else unset ANNOBTD_LINEAGE_COUNTS; fi
    printf "target\t%s\tANNOBTD_PROFILE=%s\tANNOBTD_FAMILY=%s\tANNOBTD_ORDER=%s\tANNOBTD_GENUS=%s\n" "$target" "${ANNOBTD_PROFILE:-unset}" "$ANNOBTD_FAMILY" "$ANNOBTD_ORDER" "$ANNOBTD_GENUS" >> "$CONFIG"

    guides=""
    if [ -n "$DB" ]; then
        # Guides come from the reference database, chosen by k-mer similarity to
        # this genome rather than from the handful of other test genomes. The
        # target's own accession is excluded so leave-one-out stays honest.
        acc=$(grep -m1 '^>' "$DATA/$target.fsa" | sed 's/^>\([^ .]*\).*/\1/')
        # Leave-one-out must also leave out the target's OTHER accessions: a RefSeq
        # record and the INSDC record it was derived from are the same genome, and
        # the full database holds both (Pinus thunbergii NC_001631 / D17510).
        SU_EX="${ANNOBTD_SPECIES_UNITS:-$(pwd)/species_units.tsv}"
        if [ -s "$SU_EX" ]; then acc=$(awk -F'\t' -v b="$acc" 'NR>1{a=$1; sub(/\.[0-9]+$/,"",a); hit=(a==b); if(!hit && $9!="-"){n=split($9,m,","); for(i=1;i<=n;i++){x=m[i]; sub(/\.[0-9]+$/,"",x); if(x==b) hit=1}} if(hit){out=$1; if($9!="-") out=out","$9; print out; exit}}' "$SU_EX"); acc="${acc:-$(grep -m1 '^>' "$DATA/$target.fsa" | sed 's/^>\([^ .]*\).*/\1/')}"; fi
        # Two filters sit between the similarity ranking and the guide list, each
        # asked for 3x the guides so the top N survive:
        #   ANNOBTD_GUIDE_QUALITY  accession -> outlier-CDS count (flag_profile_outliers.pl);
        #                          above ANNOBTD_MAX_OUTLIERS (3) the record is misannotated.
        #   ANNOBTD_GUIDE_LENGTHS  accession, gene, exons, prot_len (from cds_rows); each
        #                          guide's CDS lengths are compared to the TARGET lineage's
        #                          settled expectations (score_guide_consistency.pl). More
        #                          than ANNOBTD_MAX_MISMATCH (2) genes off means the guide
        #                          follows another lineage's annotation conventions and
        #                          would transfer the wrong boundaries even if correct.
        #   ANNOBTD_SPECIES_UNITS  build_species_consensus.pl units table. Candidates collapse
        #                          to ONE per species (the species' representative if present,
        #                          else the unit agreeing best with the multi-group consensus),
        #                          so a well-sampled species cannot fill the list from one lab;
        #                          RefSeq/INSDC copies of one genome count once. With
        #                          ANNOBTD_MAX_SPECIES_MISMATCH set, a unit disagreeing with its
        #                          species' consensus (>= 2 groups) on more features than that
        #                          is treated as failing.
        # Guides passing all are taken in similarity order; if fewer than N pass, the
        # rest are backfilled by fewest mismatches. Unprofiled guides pass.
        #   ANNOBTD_TAXONOMY_CHECK check_taxonomy_consistency.pl output; accessions whose
        #                          sequence places them in another order than their declared
        #                          family (verdict MISLABELED) are never guides. REVIEW rows
        #                          (same order, family disputed) pass.
        SPECIES_UNITS="${ANNOBTD_SPECIES_UNITS:-$(pwd)/species_units.tsv}"
        TAX_CHECK="${ANNOBTD_TAXONOMY_CHECK:-$(pwd)/taxonomy_consistency.tsv}"
        [ -n "${ANNOBTD_NO_TAXONOMY_CHECK:-}" ] && TAX_CHECK=""
        [ -n "${ANNOBTD_NO_SPECIES_COLLAPSE:-}" ] && SPECIES_UNITS=""
        GUIDE_LENGTHS="${ANNOBTD_GUIDE_LENGTHS:-$(pwd)/all_lengths.tsv}"
        [ -n "${ANNOBTD_NO_CONSISTENCY:-}" ] && GUIDE_LENGTHS=""
        # Select guides from one database into ranking$SUF.tsv / guide_picks$SUF.tsv
        select_from_db() { local DBSEL="$1" SUF="$2"
        if { [ -n "${ANNOBTD_GUIDE_QUALITY:-}" ] && [ -s "$ANNOBTD_GUIDE_QUALITY" ]; } || { [ -s "$GUIDE_LENGTHS" ] && [ -n "${ANNOBTD_PROFILE:-}" ]; }; then
            perl "$ANNOBTD_DIR/select_guides.pl" "$DATA/$target.fsa" "$DBSEL" -n $((NGUIDES*${ANNOBTD_CANDIDATES_MULT:-6})) $GENUSOPT --exclude "$acc" ${EXGENUS:+--exclude-genus "$ANNOBTD_GENUS"} > "$WORK/ranking_raw$SUF.tsv" 2>/dev/null
            if [ -s "$GUIDE_LENGTHS" ] && [ -n "${ANNOBTD_PROFILE:-}" ]; then
                perl "$ANNOBTD_DIR/score_guide_consistency.pl" "$WORK/ranking_raw$SUF.tsv" "$ANNOBTD_FAMILY" "$ANNOBTD_ORDER" \
                     "$ANNOBTD_PROFILE" "$GUIDE_LENGTHS" "${ANNOBTD_EXPECT_MIN_N:-10}" "${ANNOBTD_EXPECT_MIN_CONSENSUS:-0.90}" > "$WORK/ranking_scored$SUF.tsv"
            else
                awk -F'\t' '{print $0"\t0\t0\tNA\t-"}' "$WORK/ranking_raw$SUF.tsv" > "$WORK/ranking_scored$SUF.tsv"
            fi
            awk -F'\t' -v Q="${ANNOBTD_GUIDE_QUALITY:-/dev/null}" -v M="${ANNOBTD_MAX_OUTLIERS:-3}" -v C="${ANNOBTD_MAX_MISMATCH:-2}" -v N="$NGUIDES" \
                -v SU="${SPECIES_UNITS:-/dev/null}" -v SM="${ANNOBTD_MAX_SPECIES_MISMATCH:-}" -v TX="${TAX_CHECK:-/dev/null}" -v W="${ANNOBTD_SPECIES_MISMATCH_WEIGHT:-0}" '
                BEGIN{while((getline l<Q)>0){split(l,f,"\t"); q[f[1]]=f[2]}
                      while((getline l<TX)>0){split(l,f,"\t"); if(f[11]=="MISLABELED") bad[f[1]]=1}
                      while((getline l<SU)>0){split(l,f,"\t"); if(f[1]=="accession")continue; unit[f[1]]=f[1]; sp[f[1]]=f[2]; grp[f[1]]=f[5]; smis[f[1]]=f[7]; isrep[f[1]]=f[8]
                          if(f[9]!="-"){n=split(f[9],mm,","); for(i=1;i<=n;i++) unit[mm[i]]=f[1]}}}
                { if($1 in bad) next
                  u=($1 in unit)?unit[$1]:$1; if(u in seenu) next; seenu[u]=1
                  qok=(!($1 in q)||q[$1]<=M); cok=($8<=C); sok=1
                  if((u in sp) && SM!="" && grp[u]>=2 && smis[u]>SM) sok=0
                  r[NR]=$0; pass[NR]=(qok&&cok&&sok); mis[NR]=$8; usp[NR]=(u in sp)?sp[u]:"acc:"$1; urep[NR]=(u in sp)?isrep[u]:0; usm[NR]=(u in sp)?smis[u]:0 }
                END{ k=0
                     # one candidate per species: its representative, else the unit closest to
                     # the species consensus, else the most similar
                     for(i=1;i<=NR;i++){ if(!(i in r)) continue; s=usp[i]; if(!(s in best)) best[s]=i
                         else { b=best[s]; if(urep[i]>urep[b] || (urep[i]==urep[b] && usm[i]<usm[b])) best[s]=i } }
                     for(i=1;i<=NR;i++){ if(!(i in r)) continue; if(best[usp[i]]!=i) used[i]=1 }
                     # passing candidates ranked by similarity minus W x within-species disagreement
                     # (features on which the record differs from the multi-group consensus of its species;
                     # 0 when the species has one group). W=0 keeps pure similarity order.
                     while(k<N){ bi=0; bs=-1e9; for(i=1;i<=NR;i++){ if(!((i in r)&&!used[i]&&pass[i])) continue; split(r[i],ff,"\t"); sc=ff[5]-W*((grp[unit[ff[1]]]>=2)?usm[i]:0); if(sc>bs){bs=sc; bi=i} }
                                 if(bi==0) break; k++; emit(bi,k) }
                     # backfill by fewest mismatches, then similarity order
                     while(k<N){ bi=0; for(i=1;i<=NR;i++) if((i in r)&&!pass[i]&&!used[i]&&(bi==0||mis[i]<mis[bi])) bi=i
                                 if(bi==0) break; k++; emit(bi,k) } }
                function emit(i,k,  f){ split(r[i],f,"\t"); used[i]=1
                     printf "%s\t%s\t%s\t%s\t%s\t%d\n", f[1],f[2],f[3],f[4],f[5],k
                     printf "%s\t%s\t%s\tjaccard=%s\tcompared=%s\tmismatches=%s\t%s\tspecies=%s\trep=%s\tspecies_mismatch=%s\t%s\n", f[1],f[3],f[4],f[5],f[7],f[8],(pass[i]?"pass":"backfill"),usp[i],urep[i],usm[i],f[10] > PICKS }' \
                PICKS="$WORK/guide_picks$SUF.tsv" "$WORK/ranking_scored$SUF.tsv" > "$WORK/ranking$SUF.tsv"
        else
            perl "$ANNOBTD_DIR/select_guides.pl" "$DATA/$target.fsa" "$DBSEL" -n "$NGUIDES" $GENUSOPT --exclude "$acc" ${EXGENUS:+--exclude-genus "$ANNOBTD_GENUS"} > "$WORK/ranking$SUF.tsv" 2>/dev/null
        fi
        }
        select_from_db "$DB" ""
        # Hybrid choice (--db2): the curated database gives better boundaries where it has
        # same-family guides (measured: 6 of 7 angiosperms), the full database where it does
        # not (Pinus, Marchantia). With at least ANNOBTD_HYBRID_MIN_FAMILY (3) same-family
        # guides passing from the primary, its list is kept and any empty slots are filled
        # from the secondary (new species only); otherwise the secondary's list is used.
        if [ -n "$DB2" ]; then
            select_from_db "$DB2" "_db2"
            fam_ok=$(awk -F'\t' -v F="$ANNOBTD_FAMILY" '$2==F && $7=="pass"' "$WORK/guide_picks.tsv" 2>/dev/null | wc -l | tr -d ' ')
            k1=$(wc -l < "$WORK/ranking.tsv" | tr -d ' ')
            if [ "$fam_ok" -ge "${ANNOBTD_HYBRID_MIN_FAMILY:-3}" ]; then
                choice="primary"; filled=0
                if [ "$k1" -lt "$NGUIDES" ]; then
                    awk -F'\t' -v N="$NGUIDES" -v k="$k1" 'NR==FNR{have[$1]=1; sp[$2"|"$3"|"$4]=1; next} !($1 in have) && k<N {k++; $6=k; print}' OFS='\t' "$WORK/ranking.tsv" "$WORK/ranking_db2.tsv" > "$WORK/ranking_fill.tsv"
                    filled=$(wc -l < "$WORK/ranking_fill.tsv" | tr -d ' '); cat "$WORK/ranking_fill.tsv" >> "$WORK/ranking.tsv"
                    awk -F'\t' 'NR==FNR{a[$1]=1; next} ($1 in a)' "$WORK/ranking_fill.tsv" "$WORK/guide_picks_db2.tsv" >> "$WORK/guide_picks.tsv"
                fi
            else
                choice="secondary"; filled=0; cp "$WORK/ranking_db2.tsv" "$WORK/ranking.tsv"; cp "$WORK/guide_picks_db2.tsv" "$WORK/guide_picks.tsv"
            fi
            printf "hybrid\t%s\tprimary_same_family_pass=%s\tprimary_picks=%s\tchoice=%s\tfilled_from_secondary=%s\n" "$target" "$fam_ok" "$k1" "$choice" "$filled" >> "$CONFIG"
        fi
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
