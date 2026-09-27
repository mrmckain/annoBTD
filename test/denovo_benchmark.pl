#!/usr/bin/perl
use strict;
use warnings;

# Benchmark denovo_exons.pl over every multi-exon protein-coding gene in the test
# set.
#
# USAGE: ./denovo_benchmark.pl [--pad N] [--gene NAME] [--genome NAME]
#
# For each gene in each genome:
#   * the locus window is the true span PADDED by --pad nucleotides on each side,
#     so the solver has to find the start and the stop itself rather than being
#     handed them
#   * the reference protein comes from the closest OTHER genome in the set, which
#     is what select_guides.pl would pick
#   * the exon count is not supplied - both 2 and 3 are tried and the better
#     scoring structure wins, so the decomposition is discovered rather than copied
#
# Results are classified rather than scored pass/fail, because a structure that
# differs from GenBank but encodes an identical protein is not an error in the
# same sense as one that does not.

my $PAD = 250;
my ($only_gene, $only_genome, $any_donor);
my $DIR = "data";
while (@ARGV) {
    my $a = shift @ARGV;
    $PAD         = shift @ARGV if $a eq '--pad';
    $only_gene   = shift @ARGV if $a eq '--gene';
    $only_genome = shift @ARGV if $a eq '--genome';
    $DIR         = shift @ARGV if $a eq '--dir';
    $any_donor   = 1            if $a eq '--any-donor';
}

my $TOOL = "../denovo_exons.pl";

my %CODON = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

# Closest other genome, by whole-plastome k-mer similarity as reported by
# rank_references_by_kmer.pl. This is the guide the pipeline would choose.
my %CLOSEST = (
    Nicotiana_tabacum    => 'Daucus_carota',
    Daucus_carota        => 'Nicotiana_tabacum',
    Arabidopsis_thaliana => 'Nicotiana_tabacum',
    Spinacia_oleracea    => 'Nicotiana_tabacum',
    Oryza_sativa         => 'Sorghum_bicolor',
    Sorghum_bicolor      => 'Zea_mays',
    Zea_mays             => 'Sorghum_bicolor',
    Pinus_thunbergii     => 'Nicotiana_tabacum',
    Marchantia_paleacea  => 'Pinus_thunbergii',
);

# With --dir the members and their closest relatives are discovered rather than
# hard coded: every genome is paired with whichever other one shares the most
# 12-mers, the same criterion select_guides.pl uses.
my @GEN;
if ($DIR ne 'data') {
    opendir(my $dh, $DIR) or die "Cannot read $DIR: $!";
    @GEN = sort map { s/\.truth\.tsv$//r } grep { /\.truth\.tsv$/ } readdir $dh;
    closedir $dh;
    %CLOSEST = ();
}
else { @GEN = sort keys %CLOSEST }

sub genome_seq {
    my $g = shift;
    my $s = '';
    open my $h, "<", "$DIR/$g.fsa" or return undef;
    while (<$h>) { chomp; next if /^>/; s/\s+//g; $s .= $_ }
    close $h;
    return uc $s;
}
# gene -> [ copies ], each copy an arrayref of exons sorted by exon number.
sub exons_of {
    my $g = shift;
    my %by;
    open my $t, "<", "$DIR/$g.truth.tsv" or return {};
    while (<$t>) {
        chomp;
        my @f = split /\t/;
        next unless ($f[4] // '') eq 'CDS';
        next unless $f[0] =~ /^(.+)_exon(\d+)$/;
        my ($gene, $n) = ($1, $2);
        next if $gene =~ /^rps12$/;          # trans-spliced: exons are not contiguous
        next if ($f[5] // '') =~ /trans_spliced/;
        push @{ $by{$gene} }, { n => $n, s => $f[1], e => $f[2], d => $f[3] };
    }
    my %out;
    for my $gene (keys %by) {
        # Split IR duplicates into separate copies on a coordinate gap.
        my @all = sort { $a->{s} <=> $b->{s} } @{ $by{$gene} };
        my @copies = ([ shift @all ]);
        for my $x (@all) {
            if ($x->{s} - $copies[-1][-1]{e} > 5000) { push @copies, [$x] }
            else { push @{ $copies[-1] }, $x }
        }
        $out{$gene} = \@copies;
    }
    return \%out;
}
sub cds_of {
    my ($seq, $copy) = @_;
    my @a = sort { $a->{s} <=> $b->{s} } @$copy;
    my $cds = '';
    $cds .= substr($seq, $_->{s} - 1, $_->{e} - $_->{s} + 1) for @a;
    if ($a[0]{d} eq '-') { $cds = reverse $cds; $cds =~ tr/ACGT/TGCA/ }
    return $cds;
}
sub prot_of {
    my $cds = shift;
    my $p = '';
    for (my $i = 0; $i + 2 < length $cds; $i += 3) {
        my $c = substr($cds, $i, 3);
        $p .= exists $CODON{$c} ? $CODON{$c} : 'X';
    }
    $p =~ s/\*$//;
    return $p;
}

my (%seq, %ex);
for my $g (@GEN) { $seq{$g} = genome_seq($g); $ex{$g} = exons_of($g) }

unless (%CLOSEST) {
    my %sk;
    for my $g (@GEN) {
        my %h;
        my $s = $seq{$g};
        for (my $i = 0; $i + 12 <= length $s; $i += 4) {
            my $k = substr($s, $i, 12);
            next if $k =~ /[^ACGT]/;
            $h{$k} = 1;
        }
        $sk{$g} = \%h;
    }
    for my $g (@GEN) {
        my ($best, $bs);
        for my $o (@GEN) {
            next if $o eq $g;
            my $sh = 0;
            for my $k (keys %{ $sk{$g} }) { $sh++ if $sk{$o}{$k} }
            if (!defined $bs || $sh > $bs) { ($best, $bs) = ($o, $sh) }
        }
        $CLOSEST{$g} = $best;
    }
}

my %tally;
my @rows;
for my $g (@GEN) {
    next if $only_genome && $g ne $only_genome;
    my $donor = $CLOSEST{$g};
    for my $gene (sort keys %{ $ex{$g} }) {
        next if $only_gene && $gene ne $only_gene;
        next unless $ex{$donor} && $ex{$donor}{$gene};

        # Reference protein: the donor's first copy of this gene, as one unit.
        my $refprot = prot_of(cds_of($seq{$donor}, $ex{$donor}{$gene}[0]));
        next unless length($refprot) > 20;

        my $copy = $ex{$g}{$gene}[0];
        my @a = sort { $a->{s} <=> $b->{s} } @$copy;
        my ($lo, $hi, $strand) = ($a[0]{s}, $a[-1]{e}, $a[0]{d});
        my $wlo = $lo - $PAD; $wlo = 1 if $wlo < 1;
        my $whi = $hi + $PAD; $whi = length($seq{$g}) if $whi > length($seq{$g});

        my $truth = join " ", map { "$_->{s}-$_->{e}" } @a;
        my $truth_prot = prot_of(cds_of($seq{$g}, $copy));

        # The reference's exon COUNT is used as a starting hint - not its
        # boundaries, which is the whole point of the exercise. A confident
        # solution at that count is accepted; otherwise the alternatives are
        # tried, so a genome that decomposes the gene differently is still found.
        my $hint = scalar @{ $ex{$donor}{$gene}[0] };
        my @try_counts = ($hint, grep { $_ != $hint } (2, 3));

        my ($best, $best_id);
        for my $n (@try_counts) {
            my @cmd = ("perl", $TOOL, "$DIR/$g.fsa", $wlo, $whi, $strand,
                       "-p", $refprot, "--exons", $n);
            push @cmd, "--any-donor" if $any_donor;
            open my $ph, "-|", @cmd or next;
            my $line = <$ph>;
            close $ph;
            next unless defined $line && $line =~ /^([\d.]+)\s+splice\s+([\d.]+)\s+(\d+) nt\s+(.*)$/;
            my ($id, $sp, $len, $co) = ($1, $2, $3, $4);
            if (!defined $best_id || $id > $best_id) { $best_id = $id; $best = { id=>$id, sp=>$sp, co=>$co, n=>$n } }
            last if $best_id >= 0.70;      # confident enough; do not pay for the rest
        }

        unless ($best) {
            push @rows, [ $g, $gene, scalar(@a), "no solution", "", "" ];
            $tally{"no solution"}++;
            next;
        }

        my @got = sort { (split /-/, $a)[0] <=> (split /-/, $b)[0] } split /\s+/, $best->{co};
        my $got = join " ", @got;
        my $verdict;
        if ($got eq $truth) { $verdict = "exact" }
        else {
            # Rebuild the found structure and compare the protein it encodes.
            my @cp = map { my ($s,$e)=split /-/; { s=>$s, e=>$e, d=>$strand } } @got;
            my $gp = prot_of(cds_of($seq{$g}, \@cp));
            $verdict = ($gp eq $truth_prot) ? "same protein" : "differs";
        }
        $tally{$verdict}++;
        push @rows, [ $g, $gene, scalar(@a), $verdict, $truth, $got ];
    }
}

printf "%-22s %-10s %-4s %-13s\n", "genome", "gene", "ex", "verdict";
for my $r (@rows) {
    printf "  %-20s %-10s %-4d %-13s\n", @$r[0,1,2,3];
    if ($r->[3] eq 'differs' || $r->[3] eq 'same protein') {
        printf "      truth %s\n      got   %s\n", $r->[4], $r->[5];
    }
}
my $tot = 0; $tot += $_ for values %tally;
print "\n  === $tot loci, pad ${PAD} nt, exon count discovered ===\n";
for my $k (qw(exact), 'same protein', 'differs', 'no solution') {
    printf "    %-14s %3d  (%.0f%%)\n", $k, $tally{$k} // 0, 100 * ($tally{$k} // 0) / ($tot || 1);
}
printf "    %-14s %3d  (%.0f%%)\n", "correct protein",
    ($tally{exact} // 0) + ($tally{'same protein'} // 0),
    100 * (($tally{exact} // 0) + ($tally{'same protein'} // 0)) / ($tot || 1);
