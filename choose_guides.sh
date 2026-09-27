#!/bin/bash
# Choose the guide genomes for a target plastome from a sketch database.
#
# USAGE: choose_guides.sh <target.fsa> <workdir> --db SKETCH_DB [--db2 SKETCH_DB]
#            -n N --family F --order O [--genus G] [--exclude ACC[,ACC..]]
#            [--exclude-genus] [--allow-congeners] [--refdb DIR]
#
# Writes <workdir>/ranking.tsv (the chosen guides, accession first) and
# <workdir>/guide_picks.tsv (why each was chosen), and prints one line
# "hybrid<TAB>..." when two databases were given. run_regression.sh and annobtd.sh
# both call this, so the benchmark and a real run choose guides the same way.
#
# Ranking is k-mer similarity (select_guides.pl), one candidate per genus unless
# --allow-congeners. Three filters then sit between the ranking and the list, each
# asked for ANNOBTD_CANDIDATES_MULT (6) x N candidates so the top N survive:
#   lineage consistency  score_guide_consistency.pl compares each candidate's CDS
#                        lengths (all_lengths.tsv) to the TARGET lineage's settled
#                        expectations (gene_expect.tsv); more than ANNOBTD_MAX_MISMATCH
#                        (2) genes off means the guide follows another lineage's
#                        annotation conventions and would transfer the wrong boundaries.
#   species collapse     species_units.tsv: one candidate per species (its representative
#                        record, else the unit agreeing best with the species consensus),
#                        so one well-sampled species cannot fill the list from one lab;
#                        RefSeq/INSDC copies of one genome count once.
#   taxonomy check       taxonomy_consistency.tsv: records whose sequence places them in
#                        another order than their declared family (MISLABELED) never guide.
# Passing candidates are taken in similarity order; empty slots are backfilled by
# fewest mismatches. With --db2 (hybrid): the primary list is kept when it holds at
# least ANNOBTD_HYBRID_MIN_FAMILY (3) passing same-family guides, empty slots filled
# from the secondary; otherwise the secondary's list is used.
#
# Reference tables come from --refdb (default ANNOBTD_REFDB, else <this dir>/refdb):
# gene_expect.tsv (or ANNOBTD_PROFILE), all_lengths.tsv (ANNOBTD_GUIDE_LENGTHS),
# species_units.tsv (ANNOBTD_SPECIES_UNITS), taxonomy_consistency.tsv
# (ANNOBTD_TAXONOMY_CHECK). ANNOBTD_NO_CONSISTENCY / _NO_SPECIES_COLLAPSE /
# _NO_TAXONOMY_CHECK switch a filter off; ANNOBTD_GUIDE_QUALITY / ANNOBTD_MAX_OUTLIERS
# (3) add the older outlier-count filter when a table is given.
set -uo pipefail
ANNOBTD_DIR="$(cd "$(dirname "$0")" && pwd)"
TARGET="${1:?usage: choose_guides.sh <target.fsa> <workdir> --db DB ...}"; WORK="${2:?workdir}"; shift 2
DB=""; DB2=""; NGUIDES=8; FAM="NA"; ORD="NA"; GEN="NA"; EXCL=""; EXGENUS=""; GENUSOPT="--one-per-genus"; REFDB="${ANNOBTD_REFDB:-$ANNOBTD_DIR/refdb}"
while [ $# -gt 0 ]; do
    case "$1" in
        --db) DB="$2"; shift 2 ;; --db2) DB2="$2"; shift 2 ;; -n) NGUIDES="$2"; shift 2 ;;
        --family) FAM="$2"; shift 2 ;; --order) ORD="$2"; shift 2 ;; --genus) GEN="$2"; shift 2 ;;
        --exclude) EXCL="$2"; shift 2 ;; --exclude-genus) EXGENUS=1; shift ;; --allow-congeners) GENUSOPT=""; shift ;;
        --refdb) REFDB="$2"; shift 2 ;;
        *) echo "choose_guides.sh: unknown option $1" >&2; exit 1 ;;
    esac
done
[ -n "$DB" ] || { echo "choose_guides.sh: --db is required" >&2; exit 1; }
mkdir -p "$WORK"
PROFILE="${ANNOBTD_PROFILE:-$REFDB/gene_expect.tsv}"; [ -s "$PROFILE" ] || PROFILE=""
SPECIES_UNITS="${ANNOBTD_SPECIES_UNITS:-$REFDB/species_units.tsv}"; [ -n "${ANNOBTD_NO_SPECIES_COLLAPSE:-}" ] && SPECIES_UNITS=""
TAX_CHECK="${ANNOBTD_TAXONOMY_CHECK:-$REFDB/taxonomy_consistency.tsv}"; [ -n "${ANNOBTD_NO_TAXONOMY_CHECK:-}" ] && TAX_CHECK=""
GUIDE_LENGTHS="${ANNOBTD_GUIDE_LENGTHS:-$REFDB/all_lengths.tsv}"; [ -n "${ANNOBTD_NO_CONSISTENCY:-}" ] && GUIDE_LENGTHS=""
# the excluded accessions also cover the other records of the same genome unit
# (a RefSeq copy and the INSDC record it came from), read from species_units.tsv
acc="$EXCL"
if [ -n "$acc" ] && [ -s "${SPECIES_UNITS:-/dev/null}" ]; then
    acc=$(echo "$acc" | tr ',' '\n' | while read -r b; do b0="${b%.*}"; awk -F'\t' -v b="$b0" 'NR>1{a=$1; sub(/\.[0-9]+$/,"",a); hit=(a==b); if(!hit && $9!="-"){n=split($9,m,","); for(i=1;i<=n;i++){x=m[i]; sub(/\.[0-9]+$/,"",x); if(x==b) hit=1}} if(hit){out=$1; if($9!="-") out=out","$9; print out; exit}}' "$SPECIES_UNITS" | grep . || echo "$b"; done | paste -sd, -)
fi
select_from_db() { local DBSEL="$1" SUF="$2"
    if { [ -n "${ANNOBTD_GUIDE_QUALITY:-}" ] && [ -s "$ANNOBTD_GUIDE_QUALITY" ]; } || { [ -s "$GUIDE_LENGTHS" ] && [ -n "$PROFILE" ]; }; then
        perl "$ANNOBTD_DIR/select_guides.pl" "$TARGET" "$DBSEL" -n $((NGUIDES*${ANNOBTD_CANDIDATES_MULT:-6})) $GENUSOPT ${acc:+--exclude "$acc"} ${EXGENUS:+--exclude-genus "$GEN"} > "$WORK/ranking_raw$SUF.tsv" 2>/dev/null
        if [ -s "$GUIDE_LENGTHS" ] && [ -n "$PROFILE" ]; then
            perl "$ANNOBTD_DIR/score_guide_consistency.pl" "$WORK/ranking_raw$SUF.tsv" "$FAM" "$ORD" \
                 "$PROFILE" "$GUIDE_LENGTHS" "${ANNOBTD_EXPECT_MIN_N:-10}" "${ANNOBTD_EXPECT_MIN_CONSENSUS:-0.90}" > "$WORK/ranking_scored$SUF.tsv"
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
                 for(i=1;i<=NR;i++){ if(!(i in r)) continue; s=usp[i]; if(!(s in best)) best[s]=i
                     else { b=best[s]; if(urep[i]>urep[b] || (urep[i]==urep[b] && usm[i]<usm[b])) best[s]=i } }
                 for(i=1;i<=NR;i++){ if(!(i in r)) continue; if(best[usp[i]]!=i) used[i]=1 }
                 while(k<N){ bi=0; bs=-1e9; for(i=1;i<=NR;i++){ if(!((i in r)&&!used[i]&&pass[i])) continue; split(r[i],ff,"\t"); sc=ff[5]-W*((grp[unit[ff[1]]]>=2)?usm[i]:0); if(sc>bs){bs=sc; bi=i} }
                             if(bi==0) break; k++; emit(bi,k) }
                 while(k<N){ bi=0; for(i=1;i<=NR;i++) if((i in r)&&!pass[i]&&!used[i]&&(bi==0||mis[i]<mis[bi])) bi=i
                             if(bi==0) break; k++; emit(bi,k) } }
            function emit(i,k,  f){ split(r[i],f,"\t"); used[i]=1
                 printf "%s\t%s\t%s\t%s\t%s\t%d\n", f[1],f[2],f[3],f[4],f[5],k
                 printf "%s\t%s\t%s\tjaccard=%s\tcompared=%s\tmismatches=%s\t%s\tspecies=%s\trep=%s\tspecies_mismatch=%s\t%s\n", f[1],f[3],f[4],f[5],f[7],f[8],(pass[i]?"pass":"backfill"),usp[i],urep[i],usm[i],f[10] > PICKS }' \
            PICKS="$WORK/guide_picks$SUF.tsv" "$WORK/ranking_scored$SUF.tsv" > "$WORK/ranking$SUF.tsv"
    else
        perl "$ANNOBTD_DIR/select_guides.pl" "$TARGET" "$DBSEL" -n "$NGUIDES" $GENUSOPT ${acc:+--exclude "$acc"} ${EXGENUS:+--exclude-genus "$GEN"} > "$WORK/ranking$SUF.tsv" 2>/dev/null
        awk -F'\t' '{print $1"\t"$3"\t"$4"\tjaccard="$5"\tcompared=-\tmismatches=-\tpass\tspecies=-\trep=-\tspecies_mismatch=-\t-"}' "$WORK/ranking$SUF.tsv" > "$WORK/guide_picks$SUF.tsv"
    fi
}
select_from_db "$DB" ""
if [ -n "$DB2" ]; then
    select_from_db "$DB2" "_db2"
    fam_ok=$(awk -F'\t' -v F="$FAM" '$2==F && $7=="pass"' "$WORK/guide_picks.tsv" 2>/dev/null | wc -l | tr -d ' ')
    k1=$(wc -l < "$WORK/ranking.tsv" | tr -d ' ')
    if [ "$fam_ok" -ge "${ANNOBTD_HYBRID_MIN_FAMILY:-3}" ]; then
        choice="primary"; filled=0
        if [ "$k1" -lt "$NGUIDES" ]; then
            awk -F'\t' -v N="$NGUIDES" -v k="$k1" 'NR==FNR{have[$1]=1; next} !($1 in have) && k<N {k++; $6=k; print}' OFS='\t' "$WORK/ranking.tsv" "$WORK/ranking_db2.tsv" > "$WORK/ranking_fill.tsv"
            filled=$(wc -l < "$WORK/ranking_fill.tsv" | tr -d ' '); cat "$WORK/ranking_fill.tsv" >> "$WORK/ranking.tsv"
            awk -F'\t' 'NR==FNR{a[$1]=1; next} ($1 in a)' "$WORK/ranking_fill.tsv" "$WORK/guide_picks_db2.tsv" >> "$WORK/guide_picks.tsv"
        fi
    else
        choice="secondary"; filled=0; cp "$WORK/ranking_db2.tsv" "$WORK/ranking.tsv"; cp "$WORK/guide_picks_db2.tsv" "$WORK/guide_picks.tsv"
    fi
    printf "hybrid\tprimary_same_family_pass=%s\tprimary_picks=%s\tchoice=%s\tfilled_from_secondary=%s\n" "$fam_ok" "$k1" "$choice" "$filled"
fi
[ -s "$WORK/ranking.tsv" ] || { echo "choose_guides.sh: no guides selected from $DB" >&2; exit 2; }
