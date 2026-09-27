#!/usr/bin/perl
use strict;
use warnings;

# Rank a sketch database against the genome being annotated and emit the closest
# guides.
#
# USAGE: select_guides.pl <target.fsa> <sketch_db.tsv> [-n 10] [--min-jaccard 0]
#                         [--exclude ACC,ACC] [--one-per-genus] [--exclude-genus GENUS]
#
# OUT: accession <TAB> organism <TAB> family <TAB> order <TAB> jaccard <TAB> rank
#
# Exon structure is lineage specific, so the guide a genome is annotated against
# should be a relative rather than whichever reference happens to be on hand.
# Ranking by k-mer similarity is a cheap proxy for relatedness that needs no
# taxonomy lookup and degrades gracefully for anything unplaced - though the
# organism and family columns are carried through so the choice can be sanity
# checked against taxonomy.
#
# --one-per-genus keeps the guide set from filling up with near-identical
# congeners, which adds reference bulk without adding information.

my ($target, $db, @rest) = @ARGV;
die "usage: $0 <target.fsa> <sketch_db.tsv> [-n N] [--min-jaccard X] " .
    "[--exclude ACC,...] [--one-per-genus] [--exclude-genus GENUS]\n" unless defined $db;

my $N          = 10;
my $MIN_J      = 0;
my $one_genus  = 0;
my $exclude_genus = '';   # benchmark: no guide from the target's own genus (a novel lineage)
my %exclude;
while (@rest) {
    my $a = shift @rest;
    if    ($a eq '-n')             { $N     = shift @rest }
    elsif ($a eq '--min-jaccard')  { $MIN_J = shift @rest }
    elsif ($a eq '--one-per-genus'){ $one_genus = 1 }
    elsif ($a eq '--exclude')      { $exclude{$_} = 1 for split /,/, (shift @rest // '') }
    elsif ($a eq '--exclude-genus'){ $exclude_genus = shift @rest // '' }
}

my $K           = $ENV{ANNOBTD_KMER_K}      || 21;
my $SKETCH_SIZE = $ENV{ANNOBTD_SKETCH_SIZE} || 2000;

open my $tf, "<", $target or die "Cannot open $target: $!";
my $seq = '';
while (<$tf>) { chomp; next if /^>/; s/\s+//g; $seq .= $_ }
close $tf;
die "No sequence in $target\n" unless length $seq;
my $tsketch = sketch(uc $seq, $K, $SKETCH_SIZE);
die "Target produced no usable ${K}-mers\n" unless @$tsketch;

open my $dh, "<", $db or die "Cannot open sketch db $db: $!";
my @hits;
while (<$dh>) {
    chomp;
    my ($acc, $org, $fam, $ord, $len, $hashes) = split /\t/;
    next unless defined $hashes;
    next if $exclude{$acc};
    if ($exclude_genus ne '') { my @w = grep { $_ ne 'x' } split /\s+/, ($org // ''); next if @w && $w[0] eq $exclude_genus }
    # Also drop an accession whose base (version stripped) was excluded.
    my $base = $acc; $base =~ s/\.\d+$//;
    next if $exclude{$base};
    my @s = split /,/, $hashes;
    next unless @s;
    my $j = sketch_jaccard($tsketch, \@s, $SKETCH_SIZE);
    next if $j < $MIN_J;
    push @hits, { acc => $acc, org => $org, fam => $fam, ord => $ord, j => $j };
}
close $dh;
die "Sketch database $db held no usable entries\n" unless @hits;

@hits = sort { $b->{j} <=> $a->{j} || $a->{acc} cmp $b->{acc} } @hits;

my (%genus_seen, @keep);
for my $h (@hits) {
    if ($one_genus) {
        my ($genus) = split /\s+/, ($h->{org} // '');
        next if defined $genus && length $genus && $genus_seen{$genus}++;
    }
    push @keep, $h;
    last if @keep >= $N;
}

my $rank = 0;
for my $h (@keep) {
    $rank++;
    printf "%s\t%s\t%s\t%s\t%.6f\t%d\n",
        $h->{acc}, $h->{org}, $h->{fam}, $h->{ord}, $h->{j}, $rank;
}
printf STDERR "%s: %d candidates, closest %s (%s, jaccard %.4f)\n",
    $target, scalar(@hits), $keep[0]{org}, $keep[0]{fam}, $keep[0]{j} if @keep;

# Estimated Jaccard from two bottom-N sketches: merge the two sorted lists, take
# the N smallest distinct hashes overall, and count how many are in both.
sub sketch_jaccard {
    my ($a, $b, $n) = @_;
    my %in_a = map { $_ => 1 } @$a;
    my %in_b = map { $_ => 1 } @$b;
    my ($i, $j, $shared, $seen) = (0, 0, 0, 0);
    my %done;
    while ($seen < $n && ($i < @$a || $j < @$b)) {
        my $v;
        if    ($j >= @$b) { $v = $a->[$i++] }
        elsif ($i >= @$a) { $v = $b->[$j++] }
        else { $v = $a->[$i] <= $b->[$j] ? $a->[$i++] : $b->[$j++] }
        next if $done{$v}++;
        $seen++;
        $shared++ if $in_a{$v} && $in_b{$v};
    }
    return $seen ? $shared / $seen : 0;
}

sub sketch {
    my ($s, $k, $n) = @_;
    my (%seen, @v);
    my $len = length($s) - $k;
    for my $i (0 .. $len) {
        my $kmer = substr($s, $i, $k);
        next if $kmer =~ /[^ACGT]/;
        my $rc = reverse $kmer;
        $rc =~ tr/ACGT/TGCA/;
        my $h = hash32($kmer le $rc ? $kmer : $rc);
        next if $seen{$h}++;
        push @v, $h;
    }
    @v = sort { $a <=> $b } @v;
    @v = @v[0 .. $n - 1] if @v > $n;
    return \@v;
}

sub hash32 {
    my $s = shift;
    my $h = 2166136261;
    for my $c (unpack 'C*', $s) {
        $h ^= $c;
        $h = ($h * 16777619) & 0xFFFFFFFF;
    }
    return $h;
}
