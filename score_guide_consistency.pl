#!/usr/bin/perl
# Lineage-consistency of candidate guides: does each guide's annotation follow the
# conventions of the TARGET's lineage? For every gene and exon the guide annotates
# that has a settled expectation in the target's family (or order, when the family
# is thin), the guide's length is compared to the expected median. This is a
# different question from guide quality (is the record misannotated against its
# OWN family): a correctly annotated guide from a family with a different start
# convention transfers the wrong boundary.
#
# A mismatch is:
#   gene   length more than 2% off the settled median
#   exon   (only with ANNOBTD_CONSISTENCY_EXON_LENGTH=1) length more than 2% off,
#          OR off by a non-multiple of 3. Off by default: measured on the 9-genome
#          DB-guided set it cost recall (1090 -> 1081 both-exact, missing 92 -> 113)
#          because a cell "settled" at 2% tolerance hides one-nucleotide splits
#          (Solanaceae atpF exon 2: 410 in most records, median 411), so the rule
#          fires on whichever convention the median missed and evicts close guides.
#   exons  (only with ANNOBTD_CONSISTENCY_EXON_COUNT=1) the guide's exon count
#          differs from the lineage's modal count (a two-exon rps12 in a three-exon
#          family loses the 26-nt exon). Off by default: measured neutral overall
#          (1090 -> 1088 both-exact; Nicotiana +, Oryza -).
#
# USAGE: score_guide_consistency.pl <ranking.tsv> <family> <order> <gene_expect.tsv> <guide_lengths.tsv> [min_n=20] [min_cons=0.90]
#   ranking.tsv:       select_guides.pl output (acc, organism, family, order, jaccard, rank)
#   guide_lengths.tsv: extract_guide_lengths.pl output (accession, name, exons, nt)
#   gene_profile.tsv (exon-count modes) is looked for beside gene_expect.tsv
#   (gene_expectN.tsv -> gene_profileN.tsv); without it the exon-count rule is off.
# OUT: the ranking columns + compared, mismatches, mismatch_frac, mismatched_items
use strict; use warnings;
my ($rank, $fam, $ord, $exp, $lens, $min_n, $min_cons) = @ARGV;
die "usage: $0 <ranking.tsv> <family> <order> <gene_expect.tsv> <guide_lengths.tsv> [min_n] [min_cons]\n" unless $lens;
$min_n //= 10; $min_cons //= 0.90;   # a vote is a species since gene_expect6

# Expected length per gene/exon for the target lineage, with the same fallthrough
# rule as match_orfs expected_nt(): family if settled; order only if the family
# cell is thin (n < min_n); a populated but unsettled family cell blocks the item.
my %E; open my $e, "<", $exp or die "$exp: $!"; <$e>;
while (<$e>) { chomp; my ($g, $lt, $lin, $n, $med, $cons) = split /\t/; next if $lin eq 'NA';
    $E{$g}{"$lt:$lin"} = [$med, $n, $cons] }
close $e;
my %want;
for my $g (keys %E) {
    for my $lin ("F:$fam", "O:$ord") {
        my $c = $E{$g}{$lin}; next unless $c && $c->[1] >= $min_n;
        $want{$g} = $c->[0] if $c->[2] >= $min_cons;
        last;
    }
}
# Modal exon count per gene for the target lineage (family, else order).
my %exmode;
(my $prof = $exp) =~ s/gene_expect(\w*)\.tsv$/gene_profile$1.tsv/;
if ($prof ne $exp && -s $prof) {
    open my $p, "<", $prof or die; my $hdr = <$p>; chomp $hdr; my @h = split /\t/, $hdr;
    my %ix; @ix{@h} = 0 .. $#h;
    my (%byfam, %byord);
    while (<$p>) { chomp; my @f = split /\t/; next unless $f[$ix{n}] >= $min_n;
        $byfam{$f[0]} = $f[$ix{exons_mode}] if $f[1] eq $fam;
        $byord{$f[0]} = $f[$ix{exons_mode}] if $f[2] eq $ord && $f[1] ne 'NA' }
    close $p;
    # order-level mode only where the family is absent (thin families are not in
    # the profile at all below min_n), and only when the order's families agree
    my %oc; for my $g (keys %byord) { $oc{$g}{$byord{$g}}++ }
    for my $g (keys %byfam) { $exmode{$g} = $byfam{$g} }
    for my $g (keys %oc) { next if exists $exmode{$g}; my @m = keys %{ $oc{$g} }; $exmode{$g} = $m[0] if @m == 1 }
}

my %acc; open my $r, "<", $rank or die "$rank: $!"; my @rows;
while (<$r>) { chomp; next unless /\S/; my @f = split /\t/; push @rows, \@f; $acc{$f[0]} = 1 }
close $r;

my (%L, %X); open my $l, "<", $lens or die "$lens: $!";
while (<$l>) { chomp; my ($a, $g, $ex, $nt) = split /\t/; next if $a eq 'accession' || !$acc{$a};
    $L{$a}{$g} = $nt; $X{$a}{$g} = $ex if $g !~ /_exon\d+$/ }
close $l;

for my $f (@rows) {
    my ($a) = @$f; my ($n, $bad, @bad) = (0, 0);
    if ($L{$a}) {
        for my $g (sort keys %{ $L{$a} }) {
            my $m = $want{$g}; my $len = $L{$a}{$g};
            $m = undef if $g =~ /_exon\d+$/ && !$ENV{ANNOBTD_CONSISTENCY_EXON_LENGTH};
            if (defined $m) {
                $n++; my $d = $len - $m; my $off = abs($d) / $m > 0.02;
                $off ||= ($g =~ /_exon\d+$/ && abs($d) % 3);
                if ($off) { $bad++; push @bad, "$g:$len/$m" }
            }
            if ($g !~ /_exon\d+$/ && defined $exmode{$g} && $ENV{ANNOBTD_CONSISTENCY_EXON_COUNT}) {
                $n++; if ($X{$a}{$g} != $exmode{$g}) { $bad++; push @bad, "$g:${\$X{$a}{$g}}exons/$exmode{$g}" }
            }
        }
    }
    printf "%s\t%d\t%d\t%s\t%s\n", join("\t", @$f), $n, $bad, ($n ? sprintf("%.3f", $bad / $n) : "NA"),
        (@bad ? join(",", @bad) : "-");
}
