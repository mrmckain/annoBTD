#!/bin/bash
# Rebuild the derived reference tables from the GenBank flat-file cache.
#
# Everything the regression harness and the annotation pipeline read besides the
# scripts themselves is derived from two sources: the taxonomy table of every
# complete plastome in GenBank (../../taxonomy_all.tsv) and the flat-file cache
# (gb_cache/, one <accession.version>.gb per record). This script re-derives the
# rest, in dependency order, into a staging directory, and installs the staged
# set over the live one only when asked. The live files keep plain names; the
# previous set moves to archive/refs_<date>/ so a regression can be rerun
# against either.
#
# USAGE: ./refresh_references.sh [--offline] [--stage S[,S...]] [--force]
#                                [--workers N] [--install] [--stage-dir DIR]
#
#   default        every stage below, in order, skipping any whose staged output
#                  is newer than its inputs
#   --offline      skip the two stages that talk to NCBI (taxonomy, mine)
#   --stage S,...  run only the named stages (still in order)
#   --force        rebuild a stage even if its staged output looks current
#   --workers N    parallel sketch workers (default 8)
#   --install      after building, move the staged set into place; the current
#                  set goes to archive/refs_<date>/ and a summary of what changed
#                  is printed. Rerun the regression arms before trusting it
#                  (README, "Current reference numbers").
#
# STAGES and what each produces:
#   taxonomy   build_taxonomy_all.pl        new rows appended to ../../taxonomy_all.tsv
#                                           (existing rows, including assigned lineages, kept)
#   mine       mine_genbank_batch.sh        new flat files into gb_cache/ (resumable)
#   meta       extract_record_meta.pl       record_meta.tsv     submitter group + RefSeq source per record
#   lengths    extract_guide_lengths.pl     all_lengths.tsv     gene and exon lengths per record
#   sketches   build_reference_sketches.pl  ../reference_sketches_full.tsv, appended IN PLACE
#                                           for accessions not yet sketched (N workers)
#   taxcheck   check_taxonomy_consistency.pl taxonomy_consistency.tsv  MISLABELED / REVIEW / ASSIGNABLE
#   assign     assign_lineages.pl           lineage_assignments.tsv; ASSIGNABLE lineages written
#                                           IN PLACE into the taxonomy table and the sketch DB
#   consensus  build_species_consensus.pl   species_units.tsv species_lengths.tsv species_consensus.tsv
#                                           (one vote per submitter group, one form per species,
#                                           MISLABELED records excluded)
#   expect     build_expect_from_lengths.pl gene_expect.tsv gene_profile.tsv (species-voted,
#                                           genus / family / order / all strata)
#   counts     awk                          lineage_species_counts.tsv  species per family and
#                                           order (presence prior denominator)
#
# The in-place stages (sketches, assign) touch shared files deliberately: the
# sketch DB is append-only and keyed by accession, and assigned lineages are a
# correction to the taxonomy table itself, not a derived view of it. Both keep
# backups (assign: *.pre_assign, moved to the archive before a reapply). When an
# assignment lands, the sketch DB is newer than the staged taxonomy check, so the
# next run redoes taxcheck (about 7 minutes) with the corrected lineages; that
# second pass is what --install should follow.

set -uo pipefail
cd "$(dirname "$0")"
TEST="$(pwd)"
ANNOBTD_DIR="$(cd .. && pwd)"
TOP="$(cd ../.. && pwd)"
TAX="$TOP/taxonomy_all.tsv"
CACHE="$TEST/gb_cache"
SKETCH="$ANNOBTD_DIR/reference_sketches_full.tsv"
STAGE="$TEST/refresh.new"
ARCHIVE="$TEST/archive"
WORKERS=8
FORCE=""; INSTALL=""; OFFLINE=""; ONLY=""
ALL_STAGES="taxonomy mine meta lengths sketches taxcheck assign consensus expect counts"

while [ $# -gt 0 ]; do
    case "$1" in
        --offline)   OFFLINE=1; shift ;;
        --stage)     ONLY="$2"; shift 2 ;;
        --force)     FORCE=1; shift ;;
        --workers)   WORKERS="$2"; shift 2 ;;
        --install)   INSTALL=1; shift ;;
        --stage-dir) STAGE="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done
mkdir -p "$STAGE"
LOG="$STAGE/refresh.log"
say() { printf '%s  %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG" >&2; }
die() { say "ERROR: $*"; exit 1; }

# Does stage S run this time?
wanted() { case ",$ONLY," in ,,) [ -z "$OFFLINE" ] || { [ "$1" != taxonomy ] && [ "$1" != mine ]; } ;; *",$1,"*) return 0 ;; *) return 1 ;; esac; }
# Is OUT newer than every IN (and non-empty)? Then the stage is current.
current() { local out="$1"; shift; [ -n "$FORCE" ] && return 1; [ -s "$out" ] || return 1
    for i in "$@"; do [ -e "$i" ] || continue; [ "$out" -nt "$i" ] || return 1; done; return 0; }
rows() { [ -s "$1" ] && tail -n +2 "$1" | wc -l | tr -d ' ' || echo 0; }

[ -s "$TAX" ]  || die "no taxonomy table at $TAX (build_taxonomy_all.pl)"
[ -d "$CACHE" ] || die "no flat-file cache at $CACHE (mine_genbank_batch.sh)"
say "refresh: taxonomy $(rows "$TAX") rows, cache $(ls "$CACHE" | grep -c '\.gb$') records, sketches $(rows "$SKETCH"), staging $STAGE"

# ---- taxonomy: new accessions in GenBank --------------------------------------
if wanted taxonomy; then
    say "[taxonomy] build_taxonomy_all.pl (esearch + esummary + taxonomy efetch for new records)"
    perl "$ANNOBTD_DIR/build_taxonomy_all.pl" "$TAX" 2>&1 | tee -a "$LOG" >&2
    say "[taxonomy] $(rows "$TAX") rows"
fi

# ---- mine: flat files for anything not cached ---------------------------------
if wanted mine; then
    say "[mine] mine_genbank_batch.sh (skips cached records)"
    bash "$ANNOBTD_DIR/mine_genbank_batch.sh" "$TAX" "$CACHE" 2>&1 | tail -3 | tee -a "$LOG" >&2
    say "[mine] cache now $(ls "$CACHE" | grep -c '\.gb$') records"
fi

# Accessions the derived tables are built from: in the taxonomy table AND cached.
ACC="$STAGE/accessions.txt"
awk -F'\t' 'NR>1{sub(/\r$/,"",$1); print $1}' "$TAX" | while read -r a; do [ -s "$CACHE/$a.gb" ] && echo "$a"; done | sort > "$ACC"
say "accessions with a cached record: $(wc -l < "$ACC" | tr -d ' ') of $(rows "$TAX") in the taxonomy table"
[ -s "$ACC" ] || die "nothing to build from"

# ---- meta ---------------------------------------------------------------------
if wanted meta; then
    if current "$STAGE/record_meta.tsv" "$ACC" "$ANNOBTD_DIR/extract_record_meta.pl"; then say "[meta] current, skipped"; else
        say "[meta] extract_record_meta.pl over the cache"
        perl "$ANNOBTD_DIR/extract_record_meta.pl" "$CACHE" > "$STAGE/record_meta.tsv.part" && mv "$STAGE/record_meta.tsv.part" "$STAGE/record_meta.tsv" || die "meta failed"
        say "[meta] $(rows "$STAGE/record_meta.tsv") records"
    fi
fi

# ---- lengths ------------------------------------------------------------------
if wanted lengths; then
    if current "$STAGE/all_lengths.tsv" "$ACC" "$ANNOBTD_DIR/extract_guide_lengths.pl" "$ANNOBTD_DIR/gene_synonyms.tsv"; then say "[lengths] current, skipped"; else
        say "[lengths] extract_guide_lengths.pl over $(wc -l < "$ACC" | tr -d ' ') records"
        perl "$ANNOBTD_DIR/extract_guide_lengths.pl" "$CACHE" "$ACC" > "$STAGE/all_lengths.tsv.part" && mv "$STAGE/all_lengths.tsv.part" "$STAGE/all_lengths.tsv" || die "lengths failed"
        say "[lengths] $(rows "$STAGE/all_lengths.tsv") feature rows"
    fi
fi

# ---- sketches: append the missing accessions to the live sketch DB ------------
if wanted sketches; then
    [ -s "$SKETCH" ] || { printf '' > "$SKETCH"; }
    MISSING="$STAGE/sketch_missing.txt"
    comm -23 "$ACC" <(awk -F'\t' '{print $1}' "$SKETCH" | sort) > "$MISSING"
    nmiss=$(wc -l < "$MISSING" | tr -d ' ')
    if [ "$nmiss" = 0 ]; then say "[sketches] every cached accession is sketched"; else
        say "[sketches] $nmiss accessions to sketch with $WORKERS workers"
        P="$STAGE/sketch_parts"; rm -rf "$P"; mkdir -p "$P"
        per=$(( (nmiss + WORKERS - 1) / WORKERS ))
        split -l "$per" -a 2 -d "$MISSING" "$P/acc."
        for f in "$P"/acc.??; do
            perl "$ANNOBTD_DIR/build_reference_sketches.pl" "$TAX" "$f.db" --only-accessions "$f" --gb-dir "$CACHE" > "$f.log" 2>&1 &
        done
        wait
        n=0; for f in "$P"/acc.??.db; do [ -s "$f" ] || continue; cat "$f" >> "$SKETCH"; n=$((n + $(wc -l < "$f"))); done
        say "[sketches] appended $n sketches; $(grep -h SKIP "$P"/acc.??.log | wc -l | tr -d ' ') skipped (see $P/*.log); DB now $(rows "$SKETCH")"
    fi
fi

# ---- taxcheck -----------------------------------------------------------------
if wanted taxcheck; then
    if current "$STAGE/taxonomy_consistency.tsv" "$SKETCH" "$ANNOBTD_DIR/check_taxonomy_consistency.pl"; then say "[taxcheck] current, skipped"; else
        say "[taxcheck] check_taxonomy_consistency.pl over $(rows "$SKETCH") sketches"
        perl "$ANNOBTD_DIR/check_taxonomy_consistency.pl" "$SKETCH" "$STAGE/taxonomy_consistency.tsv.part" >> "$LOG" 2>&1 && mv "$STAGE/taxonomy_consistency.tsv.part" "$STAGE/taxonomy_consistency.tsv" || die "taxcheck failed"
        say "[taxcheck] verdicts: $(awk -F'\t' 'NR>1{v[$11]++} END{for(k in v) printf "%s=%d ",k,v[k]}' "$STAGE/taxonomy_consistency.tsv")"
    fi
fi

# ---- assign: proposed lineages for records that declare none -----------------
if wanted assign; then
    CHK="$STAGE/taxonomy_consistency.tsv"; [ -s "$CHK" ] || CHK="$TEST/taxonomy_consistency.tsv"
    # only records whose taxonomy row still says NA are pending (a rerun after an
    # install would otherwise re-apply the same assignments)
    nas=$(awk -F'\t' 'NR==FNR{ if (FNR>1 && $11=="ASSIGNABLE") a[$1]=1; next } FNR>1 && ($1 in a) && $6=="NA"' "$CHK" "$TAX" | wc -l | tr -d ' ')
    if [ "$nas" = 0 ]; then say "[assign] no pending ASSIGNABLE records in $CHK"; else
        say "[assign] $nas assignable records: writing lineages into the taxonomy table and the sketch DB"
        B="$ARCHIVE/pre_assign_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$B"
        for f in "$TAX.pre_assign" "$SKETCH.pre_assign"; do [ -e "$f" ] && mv "$f" "$B/"; done
        perl "$ANNOBTD_DIR/assign_lineages.pl" "$CHK" "$TAX" "$SKETCH" "$STAGE/lineage_assignments.new.tsv" 2>&1 | tee -a "$LOG" >&2 || die "assign failed"
        # the audit table is cumulative: earlier assignments stay on record
        { head -1 "$STAGE/lineage_assignments.new.tsv"
          { [ -s "$TEST/lineage_assignments.tsv" ] && tail -n +2 "$TEST/lineage_assignments.tsv"; tail -n +2 "$STAGE/lineage_assignments.new.tsv"; } | awk -F'\t' '!s[$1]++' | sort
        } > "$STAGE/lineage_assignments.tsv"
        say "[assign] audit table $(rows "$STAGE/lineage_assignments.tsv") rows (backups in $B)"
    fi
fi

# ---- consensus ----------------------------------------------------------------
if wanted consensus; then
    META="$STAGE/record_meta.tsv"; [ -s "$META" ] || META="$TEST/record_meta.tsv"
    LENS="$STAGE/all_lengths.tsv"; [ -s "$LENS" ] || LENS="$TEST/all_lengths.tsv"
    CHK="$STAGE/taxonomy_consistency.tsv"; [ -s "$CHK" ] || CHK="$TEST/taxonomy_consistency.tsv"
    if current "$STAGE/species_units.tsv" "$META" "$LENS" "$CHK" "$ANNOBTD_DIR/build_species_consensus.pl"; then say "[consensus] current, skipped"; else
        say "[consensus] build_species_consensus.pl ($META, $LENS, excluding MISLABELED from $CHK)"
        perl "$ANNOBTD_DIR/build_species_consensus.pl" "$META" "$LENS" "$STAGE/species" --exclude "$CHK" 2>&1 | tee -a "$LOG" >&2 || die "consensus failed"
        say "[consensus] $(rows "$STAGE/species_units.tsv") units, $(awk -F'\t' 'NR>1 && $8==1' "$STAGE/species_units.tsv" | wc -l | tr -d ' ') species representatives"
    fi
fi

# ---- expect -------------------------------------------------------------------
if wanted expect; then
    SL="$STAGE/species_lengths.tsv"; [ -s "$SL" ] || SL="$TEST/species_lengths.tsv"
    if current "$STAGE/gene_expect.tsv" "$SL" "$TAX" "$ANNOBTD_DIR/build_expect_from_lengths.pl"; then say "[expect] current, skipped"; else
        say "[expect] build_expect_from_lengths.pl on $SL"
        perl "$ANNOBTD_DIR/build_expect_from_lengths.pl" "$SL" "$TAX" "$STAGE/gene_expect.tsv" "$STAGE/gene_profile.tsv" 2>&1 | tee -a "$LOG" >&2 || die "expect failed"
        say "[expect] $(rows "$STAGE/gene_expect.tsv") cells: $(awk -F'\t' 'NR>1{t[$2]++} END{for(k in t) printf "%s=%d ",k,t[k]}' "$STAGE/gene_expect.tsv")"
    fi
fi

# ---- counts -------------------------------------------------------------------
if wanted counts; then
    SU="$STAGE/species_units.tsv"; [ -s "$SU" ] || SU="$TEST/species_units.tsv"
    if current "$STAGE/lineage_species_counts.tsv" "$SU" "$TAX"; then say "[counts] current, skipped"; else
        say "[counts] species per family and order from the representatives in $SU"
        awk -F'\t' 'NR==FNR{ if (FNR>1 && $8==1) r[$1]=1; next }
                    FNR>1 { sub(/\r$/,"",$7); if ($1 in r) { f[$6]++; o[$7]++ } }
                    END{ for(k in f) print "F\t"k"\t"f[k]; for(k in o) print "O\t"k"\t"o[k] }' "$SU" "$TAX" | sort > "$STAGE/lineage_species_counts.tsv"
        say "[counts] $(wc -l < "$STAGE/lineage_species_counts.tsv" | tr -d ' ') lineage rows"
    fi
fi

# ---- compare staged vs installed ----------------------------------------------
INSTALLABLE="record_meta.tsv all_lengths.tsv taxonomy_consistency.tsv lineage_assignments.tsv species_units.tsv species_lengths.tsv species_consensus.tsv gene_expect.tsv gene_profile.tsv lineage_species_counts.tsv"
say "staged vs installed:"
for f in $INSTALLABLE; do
    [ -s "$STAGE/$f" ] || continue
    if [ -s "$TEST/$f" ]; then
        if cmp -s "$STAGE/$f" "$TEST/$f"; then st="identical"; else st="differs ($(rows "$TEST/$f") -> $(rows "$STAGE/$f") rows)"; fi
    else st="new"; fi
    say "  $f: $st"
done
if [ -s "$STAGE/gene_expect.tsv" ] && [ -s "$TEST/gene_expect.tsv" ]; then
    say "  gene_expect settled cells whose median changed: $(awk -F'\t' 'NR==FNR{ if(FNR>1 && $6>=0.9) m[$1"\t"$2"\t"$3]=$5; next } FNR>1 && $6>=0.9 && ($1"\t"$2"\t"$3 in m) && m[$1"\t"$2"\t"$3]!=$5' "$TEST/gene_expect.tsv" "$STAGE/gene_expect.tsv" | wc -l | tr -d ' ')"
fi

# ---- install ------------------------------------------------------------------
if [ -n "$INSTALL" ]; then
    B="$ARCHIVE/refs_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$B"
    n=0
    for f in $INSTALLABLE; do
        [ -s "$STAGE/$f" ] || continue
        [ -e "$TEST/$f" ] && mv "$TEST/$f" "$B/$f"
        mv "$STAGE/$f" "$TEST/$f"; n=$((n+1))
    done
    say "installed $n files; previous set in $B. Rerun the regression arms (README: Current reference numbers) before relying on them."
fi
say "done"
