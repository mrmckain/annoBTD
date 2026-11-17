#!/usr/bin/perl
# build_nterm_refs.pl - the N-terminal peptides of every protein gene in the hand-curated
# truth genomes, for the conserved-start rule (pipeline/ConservedStart.pm).
#
# USAGE: build_nterm_refs.pl <truth.tsv>... > ../tables/nterm_refs.tsv
#   each <x>.truth.tsv needs <x>.fsa (the genome) and <x>.gb (for the accession and order)
#   beside it: test/data/*.truth.tsv test/curated/*.truth.tsv
#
# OUT: gene <TAB> accession <TAB> order <TAB> organism <TAB> first 30 aa (the start codon read as M)
#
# One row per gene per genome (the first copy of an IR pair). A spliced gene is kept only
# when its first exon is at least 90 nt, so the 30 aa lie in that exon. rps12 (trans-
# spliced) is left out.
use strict; use warnings;
my %CODON; { my @b = qw(T C A G); my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; my $i = 0;
    for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1) } } } }
sub translate { my $s = shift; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr($s, $i, 3)} // 'X' } $p }
sub rc { my $x = reverse shift; $x =~ tr/ACGT/TGCA/; $x }
print join("\t", qw(gene accession order organism nterm)), "\n";
for my $t (@ARGV) { (my $base = $t) =~ s/\.truth\.tsv$//; my ($fa, $gb) = ("$base.fsa", "$base.gb"); next unless -s $fa && -s $gb;
    my $seq = ''; open my $f, "<", $fa or die; while (<$f>) { next if /^>/; s/\s//g; $seq .= uc } close $f;
    my ($acc, $org, $ord) = ('', '', '');
    open my $g, "<", $gb or die; my $in = 0; my $lin = '';
    while (<$g>) { $acc = $1 if /^VERSION\s+(\S+)/; if (/^  ORGANISM\s+(.+)/) { $org = $1; $in = 1; next } if ($in) { if (/^\S/) { $in = 0 } else { $lin .= $_ } } last if /^FEATURES/ } close $g;
    $lin =~ s/\s+//g; ($ord) = grep { /ales$/ } split /;/, $lin; $ord //= 'NA';
    my %rows; open my $h, "<", $t or die; while (<$h>) { chomp; my @c = split /\t/; next unless ($c[4] // '') eq 'CDS'; next if $c[0] =~ /_intron/;
        my ($gene, $ex) = $c[0] =~ /^(.+?)(?:_exon(\d+))?$/; next if $gene =~ /^rps12/ || $gene =~ /pseudo/i; push @{ $rows{$gene} }, [$ex // 1, @c[1 .. 3]] } close $h;
    for my $gene (sort keys %rows) { my @r = @{ $rows{$gene} };
        # one copy: the exons of one strand, nearest each other (IR copies sit far apart)
        my %by; push @{ $by{ $_->[0] } }, $_ for @r; my @ex;
        for my $n (sort { $a <=> $b } keys %by) { my ($pick) = @ex ? sort { abs($a->[1] - $ex[-1][1]) <=> abs($b->[1] - $ex[-1][1]) } grep { $_->[3] eq $ex[-1][3] } @{ $by{$n} } : @{ $by{$n} }; push @ex, $pick if $pick }
        next unless @ex; next if @ex > 1 && $ex[0][2] - $ex[0][1] + 1 < 90;
        my $x = join '', map { my $s = substr($seq, $_->[1] - 1, $_->[2] - $_->[1] + 1); $_->[3] eq '-' ? rc($s) : $s } @ex;
        my $p = translate(substr($x, 0, 93)); next if length $p < 30 || substr($p, 1, 29) =~ /\*/;
        print join("\t", $gene, $acc, $ord, $org, 'M' . substr($p, 1, 29)), "\n" } }
