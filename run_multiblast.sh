#!/bin/bash
# USAGE: run_multiblast.sh <genes.fsa> <trnas.fsa> <rrnas.fsa> <orf_prefix> <plastome.fsa>
#
# BLAST_BIN  - directory holding blastn/tblastx/makeblastdb. Defaults to the historical
#              install path, then falls back to whatever is on PATH.
# BLAST_THREADS - CPUs per search. Defaults to the machine's core count.
# BLAST_EVALUE  - e-value cutoff. Default 1e-1 preserves historical sensitivity;
#                 tighten it (1e-3, 1e-5) only once the regression harness confirms
#                 no short exons are lost.
# BLAST_MAXTARGET - -max_target_seqs cap (default 500). Must exceed the number of guide
#                 species, since one ORF legitimately hits the same gene in each.

set -euo pipefail

DEFAULT_BIN="/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin"
if [ -z "${BLAST_BIN:-}" ]; then
    if [ -x "$DEFAULT_BIN/tblastx" ]; then
        BLAST_BIN="$DEFAULT_BIN"
    elif command -v tblastx >/dev/null 2>&1; then
        BLAST_BIN="$(dirname "$(command -v tblastx)")"
    else
        echo "ERROR: cannot find BLAST+. Set BLAST_BIN to its bin directory." >&2
        exit 1
    fi
fi

if [ -z "${BLAST_THREADS:-}" ]; then
    BLAST_THREADS="$( { nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1; } )"
fi
BLAST_EVALUE="${BLAST_EVALUE:-1e-1}"
BLAST_MAXTARGET="${BLAST_MAXTARGET:-500}"

GENES="$1"; TRNAS="$2"; RRNAS="$3"; ORFPREFIX="$4"; PLASTOME="$5"

for f in "$GENES" "$TRNAS" "$RRNAS" "$ORFPREFIX" "$PLASTOME"; do
    [ -s "$f" ] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done

# Only rebuild a database when it is missing or older than its FASTA.
make_db_if_needed () {
    local fasta="$1"
    if [ ! -f "${fasta}.nin" ] || [ "$fasta" -nt "${fasta}.nin" ]; then
        "$BLAST_BIN/makeblastdb" -in "$fasta" -dbtype nucl > /dev/null
    fi
}

make_db_if_needed "$GENES"
make_db_if_needed "$TRNAS"
make_db_if_needed "$RRNAS"

echo "BLAST: bin=$BLAST_BIN threads=$BLAST_THREADS evalue=$BLAST_EVALUE max_target_seqs=$BLAST_MAXTARGET trna_task=${BLAST_TRNA_TASK:-blastn}" >&2

# -task blastn, not the default megablast. Megablast seeds on an exact 28 nt word,
# so a subject SHORTER than that can never be found - and spliced tRNA first exons
# are 23 nt. trnG_exon1 is 23 nt in every genome in the test set that has it, and
# megablast returned zero hits for it in all of them, so the exon was never
# annotated even when the reference was the genome's own sequence. With -task
# blastn (word size 11) it hits at 100% identity over 23/23 at the exact truth
# coordinates. Measured on Zea mays the switch is a strict superset: 36 references
# hit instead of 35, gaining trnG_exon1 and losing nothing.
#
# The rRNA search below is left on megablast: its subjects are 103 nt and up, well
# clear of the seed length, and megablast is the faster of the two.
BLAST_TRNA_TASK="${BLAST_TRNA_TASK:-blastn}"
"$BLAST_BIN/blastn" -query "$PLASTOME" -db "$TRNAS" -outfmt 6 -task "$BLAST_TRNA_TASK" \
    -evalue "$BLAST_EVALUE" -num_threads "$BLAST_THREADS" \
    -max_target_seqs "$BLAST_MAXTARGET" > "${ORFPREFIX}_trnas.blastn"

"$BLAST_BIN/blastn" -query "$PLASTOME" -db "$RRNAS" -outfmt 6 \
    -evalue "$BLAST_EVALUE" -num_threads "$BLAST_THREADS" \
    -max_target_seqs "$BLAST_MAXTARGET" > "${ORFPREFIX}_rrnas.blastn"

"$BLAST_BIN/tblastx" -query "$ORFPREFIX" -db "$GENES" -outfmt 6 \
    -evalue "$BLAST_EVALUE" -num_threads "$BLAST_THREADS" \
    -max_target_seqs "$BLAST_MAXTARGET" > "${ORFPREFIX}_genes.tblastx"
