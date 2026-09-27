#!/usr/bin/perl
use strict;
use warnings;

# Build a MinHash sketch database from a taxonomy table of chloroplast accessions.
#
# USAGE: build_reference_sketches.pl <taxonomy.tsv> <sketch_db.tsv> [--limit N] [--only-accessions FILE]
#
# The taxonomy table is: accession organism taxid tribe subfamily family order
# (a header line is skipped). One FASTA per accession is fetched from NCBI, reduced
# to a sketch, and appended to the sketch database.
#
# Why sketches: ranking a target against the full database by exact k-mer sets
# would mean holding ~150,000 k-mers for each of ~8,000 genomes. A MinHash sketch
# keeps the smallest N hash values instead, which estimates Jaccard similarity to
# within a few percent while making the whole database a few megabytes.
#
# The run is RESUMABLE: accessions already present in the sketch database are
# skipped, so it can be stopped and restarted. Expect roughly one accession per
# second against NCBI's rate limit.
#
# Sketch line format:
#   accession <TAB> organism <TAB> family <TAB> order <TAB> genome_len <TAB> h1,h2,...

my ($tax, $db, @rest) = @ARGV;
die "usage: $0 <taxonomy.tsv> <sketch_db.tsv> [--limit N] [--only-accessions FILE] [--gb-dir DIR]\n"
    unless defined $db;

my ($limit, $only_file, $gb_dir);
while (@rest) {
    my $a = shift @rest;
    $limit     = shift @rest if $a eq '--limit';
    $only_file = shift @rest if $a eq '--only-accessions';
    $gb_dir    = shift @rest if $a eq '--gb-dir';    # read sequences from cached flat files, no NCBI
}

my $K            = $ENV{ANNOBTD_KMER_K}      || 21;
my $SKETCH_SIZE  = $ENV{ANNOBTD_SKETCH_SIZE} || 2000;
my $EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';

my %only;
if ($only_file) {
    open my $fh, "<", $only_file or die "Cannot open $only_file: $!";
    while (<$fh>) { chomp; s/\s.*$//; $only{$_} = 1 if length }
    close $fh;
    printf STDERR "restricting to %d accessions from %s\n", scalar(keys %only), $only_file;
}

# Already-sketched accessions, so a stopped run can be resumed.
my %have;
if (-s $db) {
    open my $fh, "<", $db or die "Cannot read $db: $!";
    while (<$fh>) { $have{$1} = 1 if /^(\S+)\t/ }
    close $fh;
    printf STDERR "resuming: %d accessions already sketched\n", scalar(keys %have);
}

open my $tf, "<", $tax or die "Cannot open $tax: $!";
my @todo;
while (<$tf>) {
    chomp;
    next if /^accession\b/;
    my @f = split /\t/;
    next unless @f >= 7;
    next if $have{ $f[0] };
    next if %only && !$only{ $f[0] };
    push @todo, \@f;
}
close $tf;
@todo = @todo[0 .. $limit - 1] if $limit && @todo > $limit;
printf STDERR "%d accessions to fetch\n", scalar @todo;
exit 0 unless @todo;

open my $out, ">>", $db or die "Cannot append to $db: $!";
select((select($out), $| = 1)[0]);

my ($done, $failed) = (0, 0);
for my $r (@todo) {
    my ($acc, $org, undef, undef, undef, $fam, $ord) = @$r;

    my $seq = '';
    if ($gb_dir) {
        if (open my $gh, "<", "$gb_dir/$acc.gb") {
            my $in = 0;
            while (<$gh>) { if ($in) { last if /^\/\//; s/[^A-Za-z]//g; $seq .= $_ } elsif (/^ORIGIN/) { $in = 1 } }
            close $gh;
        }
    } else {
        my $url = "$EUTILS?db=nuccore&id=$acc&rettype=fasta&retmode=text";
        my $fa  = qx{curl -sS -m 120 --retry 2 --retry-delay 2 '$url' 2>/dev/null};
        if (defined $fa) {
            for my $line (split /\n/, $fa) {
                next if $line =~ /^>/;
                $line =~ s/\s+//g;
                $seq .= $line;
            }
        }
    }
    if (length($seq) < 20000) {          # a plastome is ~120-180 kb
        warn "  SKIP $acc ($org): got " . length($seq) . " nt\n";
        $failed++;
        sleep 1 unless $gb_dir;
        next;
    }

    my $sk = sketch(uc $seq, $K, $SKETCH_SIZE);
    print $out join("\t", $acc, $org, $fam, $ord, length($seq), join(",", @$sk)), "\n";
    $done++;
    printf STDERR "  [%d/%d] %s %s\n", $done, scalar(@todo), $acc, $org if $done % 25 == 0;
    sleep 1 unless $gb_dir;             # stay under the NCBI rate limit
}
close $out;
printf STDERR "sketched %d, skipped %d\n", $done, $failed;

# Bottom-N MinHash over canonical k-mers. Canonicalising (each k-mer stored as the
# lexicographically smaller of itself and its reverse complement) makes the sketch
# independent of which strand the assembly was deposited on.
sub sketch {
    my ($seq, $k, $n) = @_;
    my %seen;
    my @heap;
    my $len = length($seq) - $k;
    for my $i (0 .. $len) {
        my $kmer = substr($seq, $i, $k);
        next if $kmer =~ /[^ACGT]/;
        my $rc = reverse $kmer;
        $rc =~ tr/ACGT/TGCA/;
        my $canon = $kmer le $rc ? $kmer : $rc;
        my $h = hash32($canon);
        next if $seen{$h}++;
        push @heap, $h;
    }
    @heap = sort { $a <=> $b } @heap;
    @heap = @heap[0 .. $n - 1] if @heap > $n;
    return \@heap;
}

# FNV-1a, 32 bit. Deterministic across runs and platforms, which matters because
# sketches built on different days must remain comparable.
sub hash32 {
    my $s = shift;
    my $h = 2166136261;
    for my $c (unpack 'C*', $s) {
        $h ^= $c;
        $h = ($h * 16777619) & 0xFFFFFFFF;
    }
    return $h;
}
