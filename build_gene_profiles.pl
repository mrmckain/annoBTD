#!/usr/bin/perl
use strict; use warnings;
# Mine every annotated CDS from a directory of GenBank flat files into a per-gene
# expectation: what a gene's protein usually looks like, stratified by family.
#
# USAGE: build_gene_profiles.pl <gb_dir> <taxonomy.tsv> <rows.tsv> <profile.tsv>
#
# rows.tsv    one line per CDS instance (accession, family, order, gene, exons,
#             strand, prot_len, nterm8, cterm8, start_codon, transl_except)
# profile.tsv one line per gene x family plus a gene x ALL line: n, length
#             median/min/max, modal N- and C-terminal 8-mers with their support,
#             ATG fraction, modal exon count.
#
# This is an ANNOTATION consensus, not a biological truth: nothing here has
# expression evidence behind it. Its value is statistical - a boundary error in
# one record is a minority against hundreds of records of the same gene, whereas
# transferred from a single guide it is simply copied. /pseudo features are
# excluded; /transl_except is recorded so a declared start is not read as an
# anomaly.
my ($dir, $tax, $rows_out, $prof_out, $expect_out) = @ARGV;
die "usage: $0 <gb_dir> <taxonomy.tsv> <rows.tsv> <profile.tsv> [expect.tsv]\n" unless $prof_out;
# expect.tsv is the compact table the annotator consumes (ANNOBTD_PROFILE):
#   name  lineage_type(F|O|ALL)  lineage  n  median_nt
# where name is a whole gene (nt = sum of exons, stop included) or gene_exonN in
# TRANSCRIPTION order, so a spliced gene's exons each get their own expectation.
my %expect;

my %lin;
open my $t, "<", $tax or die "$tax: $!";
while (<$t>) { chomp; my @f = split /\t/; next if $. == 1; (my $a = $f[0]) =~ s/\.\d+$//; $lin{$a} = [$f[5] // 'NA', $f[6] // 'NA'] }
close $t;

sub rc { my $x = reverse shift; $x =~ tr/ACGTacgt/TGCAtgca/; $x }

open my $R, ">", $rows_out or die;
print $R join("\t", qw(accession family order gene exons strand prot_len nterm cterm start transl_except)), "\n";
my (%prof, $nrec, $ncds);
opendir(my $dh, $dir) or die; my @files = sort grep { /\.gb$/ } readdir $dh; closedir $dh;
for my $f (@files) {
    (my $acc = $f) =~ s/\.gb$//; (my $base = $acc) =~ s/\.\d+$//;
    my ($fam, $ord) = @{ $lin{$base} || ['NA', 'NA'] };
    open my $h, "<", "$dir/$f" or next;
    my (@feat, $cur, $in_origin, $seq);
    while (<$h>) {
        chomp;
        if (/^ORIGIN/) { $in_origin = 1; next }
        if ($in_origin) { last if /^\/\//; s/[^A-Za-z]//g; $seq .= $_; next }
        if (/^     (\S+)\s+(\S.*)$/) {                 # new feature
            $cur = { key => $1, loc => $2, q => {} }; push @feat, $cur; next;
        }
        next unless $cur;
        if (/^\s{21}\/(\w+)(?:=(.*))?$/) {             # qualifier
            my ($k, $v) = ($1, $2 // ''); $v =~ s/^"//; $cur->{q}{$k} = $v; $cur->{lastq} = $k;
        } elsif (/^\s{21}(\S.*)$/) {                   # continuation
            if (defined $cur->{lastq}) { $cur->{q}{ $cur->{lastq} } .= $1 } else { $cur->{loc} .= $1 }
        }
    }
    close $h;
    $nrec++;
    $seq = uc($seq // '');
    for my $ft (@feat) {
        next unless $ft->{key} eq 'CDS';
        my $q = $ft->{q};
        next if exists $q->{pseudo} || exists $q->{pseudogene};
        my $gene = $q->{gene} or next; $gene =~ s/"//g; $gene =~ s/\s+//g;
        my $prot = $q->{translation} or next; $prot =~ s/[^A-Z*]//g;
        my $loc = $ft->{loc}; $loc =~ s/\s//g;
        my $strand = ($loc =~ /complement/) ? '-' : '+';
        my @r; while ($loc =~ /(\d+)\.\.(\d+)/g) { push @r, [$1, $2] }
        next unless @r;
        my $start = 'NA';
        if (length $seq) {
            if ($strand eq '+') { my ($s) = sort { $a->[0] <=> $b->[0] } @r; $start = substr($seq, $s->[0] - 1, 3) }
            else                { my ($e) = sort { $b->[1] <=> $a->[1] } @r; $start = rc(substr($seq, $e->[1] - 3, 3)) }
        }
        my $L = length $prot; next if $L < 20;
        # exon nt lengths in transcription order
        my @tr = $strand eq '+' ? sort { $a->[0] <=> $b->[0] } @r : sort { $b->[1] <=> $a->[1] } @r;
        my @ex_nt = map { $_->[1] - $_->[0] + 1 } @tr;
        my $tot_nt = 0; $tot_nt += $_ for @ex_nt;
        for my $lin ("F\t$fam", "O\t$ord", "ALL\tALL") { next if $lin =~ /\tNA$/;   # no NA stratum: records without taxonomy are a random mix
            push @{ $expect{"$gene\t$lin"} }, $tot_nt;
            # rps12 is trans-spliced: its 5' exon lies tens of kb from the others, so
            # sorting exons by coordinate does not give transcription order and the
            # per-exon rows came out scrambled (exon3 = 114 nt, which is exon1's length).
            # It keeps a whole-gene expectation only.
            if (@ex_nt > 1 && $gene !~ /^rps12$/) { push @{ $expect{"${gene}_exon" . ($_ + 1) . "\t$lin"} }, $ex_nt[$_] for 0 .. $#ex_nt }
        }
        my ($nt, $ct) = (substr($prot, 0, 8), substr($prot, -8));
        my $te = exists $q->{transl_except} ? 1 : 0;
        print $R join("\t", $acc, $fam, $ord, $gene, scalar(@r), $strand, $L, $nt, $ct, $start, $te), "\n";
        $ncds++;
        for my $key ("$gene\t$fam\t$ord", "$gene\tALL\tALL") { next if $key =~ /\tNA\t/;
            my $p = $prof{$key} ||= { len => [], nt => {}, ct => {}, st => {}, ex => {} };
            push @{ $p->{len} }, $L; $p->{nt}{$nt}++; $p->{ct}{$ct}++; $p->{st}{$start}++; $p->{ex}{ scalar @r }++;
        }
    }
}
close $R;
sub mode { my $h = shift; my ($k) = sort { $h->{$b} <=> $h->{$a} || $a cmp $b } keys %$h; return ($k, $h->{$k}) }
open my $P, ">", $prof_out or die;
print $P join("\t", qw(gene family order n len_median len_min len_max nterm nterm_frac cterm cterm_frac atg_frac exons_mode)), "\n";
for my $key (sort keys %prof) {
    my $p = $prof{$key}; my @l = sort { $a <=> $b } @{ $p->{len} }; my $n = @l;
    my ($nt, $ntn) = mode($p->{nt}); my ($ct, $ctn) = mode($p->{ct}); my ($ex) = mode($p->{ex});
    my $atg = ($p->{st}{ATG} // 0) / $n;
    printf $P "%s\t%d\t%d\t%d\t%d\t%s\t%.2f\t%s\t%.2f\t%.2f\t%s\n", $key, $n, $l[$n/2], $l[0], $l[-1], $nt, $ntn/$n, $ct, $ctn/$n, $atg, $ex;
}
close $P;
if ($expect_out) {
    open my $E, ">", $expect_out or die;
    # consensus = share of records within 2% of the median; alt = the commonest
    # length outside that band and its share. A cell where GenBank itself carries
    # two conventions (Poaceae atpF 57%/42%, Brassicaceae ndhB 51%/42%) has no
    # expectation worth enforcing, and the annotator refuses those below
    # ANNOBTD_EXPECT_MIN_CONSENSUS (default 0.90).
    print $E join("\t", qw(name lineage_type lineage n median_nt consensus alt_nt alt_frac)), "\n";
    for my $k (sort keys %expect) { my @v = sort { $a <=> $b } @{ $expect{$k} }; my $m = $v[@v/2];
        my $in = grep { abs($_ - $m) / ($m || 1) <= 0.02 } @v; my %c; $c{$_}++ for grep { abs($_ - $m) / ($m || 1) > 0.02 } @v;
        my ($alt) = sort { $c{$b} <=> $c{$a} || $a <=> $b } keys %c;
        printf $E "%s\t%d\t%d\t%.3f\t%s\t%.3f\n", $k, scalar @v, $m, $in / @v, $alt // 0, $alt ? $c{$alt} / @v : 0 }
    close $E;
}
printf STDERR "%d records, %d CDS instances, %d gene x lineage cells\n", $nrec, $ncds, scalar keys %prof;
