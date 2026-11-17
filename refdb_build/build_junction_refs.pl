#!/usr/bin/perl
# build_junction_refs.pl - the peptides either side of every cis-splice junction of the protein genes in the
# hand-curated truth genomes, for the conserved-junction rule (pipeline/ConservedStart.pm).
#
# USAGE: build_junction_refs.pl <truth.tsv>... > ../tables/junction_refs.tsv
#   each <x>.truth.tsv needs <x>.fsa and <x>.gb beside it (test/data, test/curated, test/asp)
#
# OUT: gene  junction(1 = between exon 1 and exon 2)  accession  order  organism  phase  left  right
#   phase: nucleotides of the split codon that lie in the upstream exon (0, 1 or 2)
#   left:  the last 15 whole codons before the junction (fewer when the exons before it are short)
#   right: the first 15 whole codons after it
# One row per junction per genome (the first copy of an IR pair); rps12 (trans-spliced) is left out.
use strict; use warnings;
my %CODON; { my @b = qw(T C A G); my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; my $i = 0;
    for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1) } } } }
sub translate { my $s = shift; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr($s, $i, 3)} // 'X' } $p }
sub rc { my $x = reverse shift; $x =~ tr/ACGT/TGCA/; $x }
print join("\t", qw(gene junction accession order organism phase left right)), "\n";
for my $t (@ARGV) { (my $base = $t) =~ s/\.truth\.tsv$//; my ($fa, $gb) = ("$base.fsa", "$base.gb"); next unless -s $fa && -s $gb;
    my $seq = ''; open my $f, "<", $fa or die; while (<$f>) { next if /^>/; s/\s//g; $seq .= uc } close $f;
    my ($acc, $org) = ('', ''); open my $g, "<", $gb or die; my $in = 0; my $lin = '';
    while (<$g>) { $acc = $1 if /^VERSION\s+(\S+)/; if (/^  ORGANISM\s+(.+)/) { $org = $1; $in = 1; next } if ($in) { if (/^\S/) { $in = 0 } else { $lin .= $_ } } last if /^FEATURES/ } close $g;
    $lin =~ s/\s+//g; my ($ord) = grep { /ales$/ } split /;/, $lin; $ord //= 'NA';
    my %rows; open my $h, "<", $t or die; while (<$h>) { chomp; my @c = split /\t/; next unless ($c[4] // '') eq 'CDS'; next unless $c[0] =~ /^(.+)_exon(\d+)$/; my ($gene, $ex) = ($1, $2);
        next if $gene =~ /^rps12/ || $gene =~ /pseudo/i; push @{ $rows{$gene} }, [$ex, @c[1 .. 3]] } close $h;
    for my $gene (sort keys %rows) { my %by; push @{ $by{ $_->[0] } }, $_ for @{ $rows{$gene} }; my @ex;
        for my $n (sort { $a <=> $b } keys %by) { my ($pick) = @ex ? sort { abs($a->[1] - $ex[-1][1]) <=> abs($b->[1] - $ex[-1][1]) } grep { $_->[3] eq $ex[-1][3] } @{ $by{$n} } : @{ $by{$n} }; push @ex, $pick if $pick }
        next unless @ex >= 2; next if grep { abs($_->[1] - $ex[0][1]) > 6000 } @ex; next if grep { $ex[$_][0] != $_ + 1 } 0 .. $#ex;
        my @part = map { my $s = substr($seq, $_->[1] - 1, $_->[2] - $_->[1] + 1); $_->[3] eq '-' ? rc($s) : $s } @ex; my $cds = join '', @part;
        next if length($cds) % 3; my $pr = translate($cds); $pr =~ s/\*$//; next if $pr =~ /\*/;
        my $P = 0; for my $k (0 .. $#ex - 1) { $P += length $part[$k]; my $ph = $P % 3; my $le = $P - $ph; my $ls = $le - 45; $ls = 0 if $ls < 0; my $rs = $P + (3 - $ph) % 3;
            my ($L, $R) = (translate(substr($cds, $ls, $le - $ls)), translate(substr($cds, $rs, 45))); $R =~ s/\*.*//;
            print join("\t", $gene, $k + 1, $acc, $ord, $org, $ph, $L, $R), "\n" } } }
