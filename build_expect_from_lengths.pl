#!/usr/bin/perl
# Per-gene expectation profile from a lengths table instead of the flat files, so
# the same cell definitions (median, consensus within 2%, alt mode) can be built
# from any set of votes: every record (reproduces build_gene_profiles.pl), or one
# vote per species (build_species_consensus.pl output), etc.
#
# USAGE: build_expect_from_lengths.pl <lengths.tsv> <taxonomy.tsv> <gene_expect.tsv> [gene_profile.tsv]
#   lengths.tsv: accession, name, exons, nt   (name = gene or gene_exonN)
#   gene_profile.tsv (optional): gene, family, order, n, exons_mode - the columns
#   score_guide_consistency.pl reads for its exon-count rule.
use strict; use warnings;
my ($lens, $tax, $exp_out, $prof_out) = @ARGV; die "usage: $0 <lengths.tsv> <taxonomy.tsv> <gene_expect.tsv> [gene_profile.tsv]\n" unless $exp_out;
# lineage per accession: family, order, and genus (first word of the organism name;
# hybrids written "x Genus" and unnamed "sp." entries still give a genus)
my %lin; open my $t, "<", $tax or die; while (<$t>) { chomp; s/\r//g; my @f = split /\t/; next if $. == 1; (my $a = $f[0]) =~ s/\.\d+$//;
    my $genus = 'NA'; if (defined $f[1]) { my @w = grep { $_ ne 'x' } split /\s+/, $f[1]; $genus = $w[0] if @w && $w[0] =~ /^[A-Z][a-z]+$/ }
    $lin{$a} = [$f[5] // 'NA', $f[6] // 'NA', $genus] } close $t;
my (%expect, %ex, %fo);
open my $l, "<", $lens or die;
while (<$l>) { chomp; my ($acc, $name, $exons, $nt) = split /\t/; next if $acc eq 'accession';
    (my $base = $acc) =~ s/\.\d+$//; my ($fam, $ord, $gen) = @{ $lin{$base} || ['NA', 'NA', 'NA'] };
    for my $lk ("G\t$gen", "F\t$fam", "O\t$ord", "ALL\tALL") { next if $lk =~ /\tNA$/;
        push @{ $expect{"$name\t$lk"} }, $nt }
    if ($name !~ /_exon\d+$/ && $fam ne 'NA') { $ex{"$name\t$fam\t$ord"}{$exons}++ }
}
close $l;
open my $E, ">", $exp_out or die;
print $E join("\t", qw(name lineage_type lineage n median_nt consensus alt_nt alt_frac)), "\n";
for my $k (sort keys %expect) { my @v = sort { $a <=> $b } @{ $expect{$k} }; my $m = $v[@v/2];
    my $in = grep { abs($_ - $m) / ($m || 1) <= 0.02 } @v; my %c; $c{$_}++ for grep { abs($_ - $m) / ($m || 1) > 0.02 } @v;
    my ($alt) = sort { $c{$b} <=> $c{$a} || $a <=> $b } keys %c;
    printf $E "%s\t%d\t%d\t%.3f\t%s\t%.3f\n", $k, scalar @v, $m, $in / @v, $alt // 0, $alt ? $c{$alt} / @v : 0 }
close $E;
if ($prof_out) {
    open my $P, ">", $prof_out or die; print $P join("\t", qw(gene family order n exons_mode)), "\n";
    for my $k (sort keys %ex) { my $h = $ex{$k}; my $n = 0; $n += $_ for values %$h;
        my ($mode) = sort { $h->{$b} <=> $h->{$a} || $a <=> $b } keys %$h; print $P join("\t", $k, $n, $mode), "\n" }
    close $P;
}
