#!/bin/bash
# Batched, resumable GenBank flat-file miner. Fetches BATCH accessions per efetch
# POST (gbwithparts), splits the response on record terminators, and writes one
# <accession.version>.gb per record. Skips anything already cached, so it can be
# stopped and restarted freely. HTTP/1.1 is forced: HTTP/2 against eutils was
# producing framing errors and multi-minute stalls.
#
# USAGE: mine_genbank_batch.sh <taxonomy.tsv> <cache_dir> [batch=50] [sleep=0.4]
set -u
TAX="$1"; CACHE="$2"; BATCH="${3:-50}"; SLEEP="${4:-0.4}"
mkdir -p "$CACHE"
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
KEY="${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}"
total=$(awk 'NR>1' "$TAX" | wc -l | tr -d ' ')
have=$(ls "$CACHE"/*.gb 2>/dev/null | wc -l | tr -d ' ')
echo "resuming: $have of $total already cached; batch=$BATCH" >&2
TMP="$CACHE/.batch.$$"; PEND="$TMP.pending"
trap 'rm -f "$TMP" "$PEND"' EXIT

# Split a multi-record flat file into per-accession files, keeping only records
# that reach ORIGIN. Prints the accessions written.
split_records() {
    perl -e '
        my ($in,$cache)=@ARGV; open my $fh,"<",$in or die; my ($acc,$buf,$ok);
        while (<$fh>) {
            if (/^LOCUS/) { $buf=""; $acc=undef; $ok=0 }
            $acc=$1 if /^VERSION\s+(\S+)/; $ok=1 if /^ORIGIN/; $buf.=$_;
            if (/^\/\/\s*$/) {
                if ($acc && $ok) { open my $o,">","$cache/$acc.gb" or die; print $o $buf, "\n"; close $o; print "$acc\n" }
                elsif ($acc) { print STDERR "incomplete: $acc\n" }
                $buf=""; $acc=undef; $ok=0;
            }
        }' "$1" "$2"
}

# Pending list goes to a file. Slicing a pipe with head over-reads its buffer
# and silently discards the rest of the queue.
awk -F'\t' 'NR>1{print $1}' "$TAX" | while read -r acc; do
    [ -s "$CACHE/$acc.gb" ] || echo "$acc"
done > "$PEND"
echo "pending: $(wc -l < "$PEND" | tr -d ' ')" >&2

n=0; done_n=0; start=1
while true; do
    ids=$(sed -n "${start},$((start+BATCH-1))p" "$PEND"); [ -z "$ids" ] && break
    start=$((start+BATCH))
    want=$(echo "$ids" | wc -l | tr -d ' ')
    idlist=$(echo "$ids" | paste -sd, -)
    got=0
    for attempt in 1 2 3; do
        if curl -sS --http1.1 -m 240 -X POST "$E" \
              --data "db=nuccore&rettype=gbwithparts&retmode=text&id=${idlist}${KEY}" -o "$TMP" \
           && grep -q '^ORIGIN' "$TMP"; then
            got=$(split_records "$TMP" "$CACHE" | wc -l | tr -d ' ')
            break
        fi
        echo "retry $attempt on batch starting $(echo "$ids" | head -1)" >&2; sleep 5
    done
    rm -f "$TMP"
    if [ "$got" -lt "$want" ]; then
        for a in $ids; do [ -s "$CACHE/$a.gb" ] || echo "missing: $a" >&2; done
    fi
    done_n=$((done_n+got)); n=$((n+1))
    echo "batch $n: $got/$want  (session total $done_n)" >&2
    sleep "$SLEEP"
done
echo "done" >&2
