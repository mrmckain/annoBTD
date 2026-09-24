#!/usr/bin/perl
# Per-record CDS lengths for a set of accessions, gene and exon level, from cached
# GenBank flat files. Same parsing and exon numbering as build_gene_profiles.pl so
# the rows line up with gene_expect.tsv cells (rps12 gets a whole-gene row only:
# trans-spliced, its exons cannot be ordered by coordinate).
#
# USAGE: extract_guide_lengths.pl <gb_dir> <accessions.txt> > guide_lengths.tsv
# OUT:   accession <TAB> name <TAB> exons <TAB> nt      (name = gene or gene_exonN)
use strict; use warnings;
my ($dir, $list) = @ARGV; die "usage: $0 <gb_dir> <accessions.txt>\n" unless $list;
my %want; open my $l, "<", $list or die; while (<$l>) { chomp; s/\r//; $want{$_} = 1 if /\S/ } close $l;
print join("\t", qw(accession name exons nt)), "\n";
for my $acc (sort keys %want) {
    open my $h, "<", "$dir/$acc.gb" or next;
    my (@feat, $cur);
    while (<$h>) {
        chomp; last if /^ORIGIN/;
        if (/^     (\S+)\s+(\S.*)$/) { $cur = { key => $1, loc => $2, q => {} }; push @feat, $cur; next }
        next unless $cur;
        if (/^\s{21}\/(\w+)(?:=(.*))?$/) { my ($k, $v) = ($1, $2 // ''); $v =~ s/^"//; $cur->{q}{$k} = $v; $cur->{lastq} = $k }
        elsif (/^\s{21}(\S.*)$/) { if (defined $cur->{lastq}) { $cur->{q}{ $cur->{lastq} } .= $1 } else { $cur->{loc} .= $1 } }
    }
    close $h;
    for my $ft (@feat) {
        next unless $ft->{key} eq 'CDS';
        my $q = $ft->{q}; next if exists $q->{pseudo} || exists $q->{pseudogene};
        my $gene = $q->{gene} or next; $gene =~ s/"//g; $gene =~ s/\s+//g;
        $gene = canonical_gene($gene);
        my $prot = $q->{translation} or next; $prot =~ s/[^A-Z*]//g; next if length($prot) < 20;
        my $loc = $ft->{loc}; $loc =~ s/\s//g;
        my $strand = ($loc =~ /complement/) ? '-' : '+';
        my @r; while ($loc =~ /(\d+)\.\.(\d+)/g) { push @r, [$1, $2] } next unless @r;
        my @tr = $strand eq '+' ? sort { $a->[0] <=> $b->[0] } @r : sort { $b->[1] <=> $a->[1] } @r;
        my @ex = map { $_->[1] - $_->[0] + 1 } @tr; my $tot = 0; $tot += $_ for @ex;
        print join("\t", $acc, $gene, scalar(@r), $tot), "\n";
        if (@ex > 1 && $gene !~ /^rps12$/) { print join("\t", $acc, "${gene}_exon" . ($_ + 1), 1, $ex[$_]), "\n" for 0 .. $#ex }
    }
}

# Same symbol normalisation as test/gb_to_truth.pl, from annoBTD/gene_synonyms.tsv
# (synonym -> canonical, case-insensitive), then case-folding of known names.
{
    my (%syn, %canon, $loaded);
    sub canonical_gene {
        my $name = shift;
        unless ($loaded) { $loaded = 1;
            use FindBin;
            if (open my $fh, "<", "$FindBin::Bin/gene_synonyms.tsv") { while (<$fh>) { chomp; next if /^#/ || !/\S/; my ($a, $b) = split /\t/; $syn{lc $a} = $b if $b } close $fh }
            $canon{lc $_} = $_ for qw(accD atpA atpB atpE atpF atpH atpI ccsA cemA chlB chlL chlN clpP infA matK ndhA ndhB ndhC ndhD
                ndhE ndhF ndhG ndhH ndhI ndhJ ndhK petA petB petD petG petL petN psaA psaB psaC psaI psaJ psaM psbA psbB psbC psbD
                psbE psbF psbH psbI psbJ psbK psbL psbM psbN psbT psbZ rbcL rpl2 rpl14 rpl16 rpl20 rpl22 rpl23 rpl32 rpl33 rpl36
                rpoA rpoB rpoC1 rpoC2 rps2 rps3 rps4 rps7 rps8 rps11 rps12 rps14 rps15 rps16 rps18 rps19 ycf1 ycf2 ycf3 ycf4 ycf12
                ycf15 ycf68 ycf66);
            $canon{lc $_} = $_ for values %syn; }
        my $k = lc $name; return $syn{$k} // $canon{$k} // do {
            # "ycf10/cemA", "ycf3/pafI": both names of one gene joined with a slash
            my @p = $k =~ m{/} ? map { $syn{$_} // $canon{$_} } grep { exists $syn{$_} || exists $canon{$_} } split m{/}, $k : ();
            (@p && !grep { $_ ne $p[0] } @p) ? $p[0] : $name };
    }
}
