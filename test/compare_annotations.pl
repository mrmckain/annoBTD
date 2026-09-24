#!/usr/bin/perl
use strict;
use warnings;

# Boundary-accuracy diff between an annoBTD result and a GenBank truth table.
#
# USAGE: compare_annotations.pl <truth.tsv> <SPECIES_VERDANT_cleaned_annotation.txt> [--detail out.tsv] [--label NAME]
#
# truth.tsv comes from gb_to_truth.pl:  gene start end strand type flags
# The annoBTD file is:                  gene start end dir
#
# Reports, per run: recall, exact-boundary rate, separate 5' and 3' accuracy, and the
# signed offset distribution. The sign convention is relative to the direction of
# transcription, so a negative 5' offset means the predicted start is UPSTREAM of truth
# on both strands - that is what makes a systematic bias visible.

my ($truth_file, $pred_file, @rest) = @ARGV;
die "usage: $0 <truth.tsv> <annoBTD_annotation.txt> [--detail out.tsv] [--label NAME]\n"
    unless defined $pred_file;

my ($detail_file, $label);
while (@rest) {
    my $a = shift @rest;
    $detail_file = shift @rest if $a eq '--detail';
    $label       = shift @rest if $a eq '--label';
}
$label ||= $pred_file;

# Regions and intergenic spacers annoBTD emits that are not genes.
sub is_gene {
    my $n = shift;
    return 0 if $n =~ /^(LSC|SSC|IRA|IRB|FULL)$/;
    return 0 if $n =~ /~/;          # intergenic spacer, e.g. "psbA~matK"
    return 0 if $n =~ /_intron\d*$/;
    return 1;
}

# Normalise names so tobacco's "trnK" matches a reference set's "trnK-UUU".
sub norm {
    my $n = lc shift;
    $n =~ s/XXX\d+$//i;             # guide-species tag
    $n =~ s/^rrn45$/rrn4.5/;
    $n =~ s/\s+$//;
    return $n;
}
sub base_key {                       # anticodon- and exon-insensitive tRNA key
    my $n = norm(shift);
    return $n unless $n =~ /^trn/;   # tRNAs only; protein exons keep their number
    $n =~ s/-[a-z]{3}(?=(_exon\d+)?$)//;   # anticodon suffix
    $n =~ s/_exon\d+$//;                   # exon suffix
    return $n;
}
# Why the exon suffix goes too, for tRNAs. gb_to_truth.pl takes names from the
# source record, and the records disagree: Daucus and Sorghum name both exons of a
# spliced tRNA "trnG-GCC", while Arabidopsis, Nicotiana and the rest name them
# "trnG_exon1" / "trnG_exon2". Keeping the exon number meant Daucus's predicted
# "trnG_exon1 9027-9049" could never pair with its truth "trnG-GCC 9027-9049" even
# though the coordinates are identical - the comparator charged it once as a miss
# and again as a spurious call. base_key is only consulted when the exact name
# fails to match, and pairing is nearest-first by coordinate, so collapsing the
# exon number cannot mis-pair exons that really are distinct.

sub load {
    my ($file, $is_truth) = @_;
    open my $fh, "<", $file or die "Cannot open $file: $!";
    my @out;
    while (<$fh>) {
        chomp;
        next unless /\S/;
        my @f = split /\t/;
        @f = split /\s+/ unless @f >= 4;
        my ($name, $s, $e, $d) = @f[0..3];
        next unless defined $d && $s =~ /^\d+$/ && $e =~ /^\d+$/;
        next unless is_gene($name);
        ($s, $e) = ($e, $s) if $e < $s;
        push @out, {
            name => $name, norm => norm($name), base => base_key($name),
            start => $s, end => $e, strand => $d,
            type => $is_truth ? ($f[4] || 'CDS') : 'NA',
            flags => $is_truth ? ($f[5] || '-') : '-',
            matched => 0,
        };
    }
    close $fh;
    return \@out;
}

my $truth = load($truth_file, 1);
my $pred  = load($pred_file,  0);

# Index truth by exact name, then by anticodon-insensitive name.
my (%by_norm, %by_base);
push @{ $by_norm{ $_->{norm} } }, $_ for @$truth;
push @{ $by_base{ $_->{base} } }, $_ for @$truth;

# Match predictions to truth rows of the same gene, nearest pair first.
#
# Iterating predictions in coordinate order was wrong: where a gene is IR
# duplicated and the prediction set holds an extra copy, the first prediction
# encountered claimed the only truth row even though a later prediction sat
# exactly on it. The correct call was then scored as an error tens of thousands
# of nucleotides out, and the spurious one counted as the match. Building every
# candidate pair and assigning the closest first removes that ordering artifact.
my (@pairs, @spurious);
my @cand;
for my $p (@$pred) {
    my $cands = $by_norm{ $p->{norm} } || $by_base{ $p->{base} };
    unless ($cands) { $p->{nogene} = 1; next }
    my $pmid = ($p->{start} + $p->{end}) / 2;
    for my $t (@$cands) {
        push @cand, [ abs($pmid - ($t->{start} + $t->{end}) / 2), $p, $t ];
    }
}
for my $c (sort { $a->[0] <=> $b->[0]
               || $a->[1]{start} <=> $b->[1]{start}
               || $a->[2]{start} <=> $b->[2]{start} } @cand) {
    my ($d, $p, $t) = @$c;
    next if $p->{matched} || $t->{matched};
    $p->{matched} = 1;
    $t->{matched} = 1;
    push @pairs, [ $p, $t ];
}
@spurious = grep { !$_->{matched} } @$pred;
my @missing = grep { !$_->{matched} } @$truth;

# Offsets in transcription orientation: negative = predicted boundary lies upstream.
my (@d5, @d3, $exact_both, $exact5, $exact3, $strand_err, $inconsistent);
my @detail;
for my $pr (@pairs) {
    my ($p, $t) = @$pr;
    my ($o5, $o3);
    if ($t->{strand} eq q{-}) {
        $o5 = $t->{end}   - $p->{end};     # neg: predicted start extends upstream
        $o3 = $t->{start} - $p->{start};   # neg: predicted end falls short
    } else {
        $o5 = $p->{start} - $t->{start};
        $o3 = $p->{end}   - $t->{end};
    }
    push @d5, $o5; push @d3, $o3;
    $exact5++ if $o5 == 0;
    $exact3++ if $o3 == 0;
    $exact_both++ if $o5 == 0 && $o3 == 0;
    $strand_err++ if $p->{strand} ne $t->{strand};
    # A gene whose two IR copies are annotated at different lengths in the source
    # record cannot be satisfied on both: one reference sequence, two extents.
    $inconsistent++ if ($o5 || $o3) && $t->{flags} =~ /inconsistent_copies/;
    push @detail, [ $t->{name}, $t->{type}, $t->{strand}, $t->{start}, $t->{end},
                    $p->{start}, $p->{end}, $p->{strand}, $o5, $o3, $t->{flags} ];
}

sub median { my @s = sort { $a <=> $b } @_; return 0 unless @s; return $s[@s/2] }
sub mean   { return 0 unless @_; my $t = 0; $t += $_ for @_; return $t / @_ }
sub within { my ($n, @v) = @_; my $c = grep { abs($_) <= $n } @v; return @v ? 100*$c/@v : 0 }

my $nt = @$truth; my $np = @$pred; my $nm = @pairs;

printf "=== %s ===\n", $label;
printf "truth features      : %d\n", $nt;
printf "predicted features  : %d\n", $np;
printf "matched             : %d  (recall %.1f%%, precision %.1f%%)\n",
    $nm, $nt ? 100*$nm/$nt : 0, $np ? 100*$nm/$np : 0;
printf "missing             : %d\n", scalar @missing;
printf "spurious            : %d\n", scalar @spurious;
printf "strand errors       : %d\n", $strand_err || 0;
print  "\n";
printf "both boundaries exact : %d / %d  (%.1f%%)\n", $exact_both||0, $nm, $nm ? 100*($exact_both||0)/$nm : 0;
printf "5' exact              : %d / %d  (%.1f%%)   within 3nt %.1f%%   within 15nt %.1f%%\n",
    $exact5||0, $nm, $nm ? 100*($exact5||0)/$nm : 0, within(3,@d5), within(15,@d5);
printf "3' exact              : %d / %d  (%.1f%%)   within 3nt %.1f%%   within 15nt %.1f%%\n",
    $exact3||0, $nm, $nm ? 100*($exact3||0)/$nm : 0, within(3,@d3), within(15,@d3);
print  "\n";
printf "5' offset  mean %+.1f   median %+d   (negative = predicted start UPSTREAM of truth)\n", mean(@d5), median(@d5);
printf "3' offset  mean %+.1f   median %+d   (negative = predicted end SHORT of truth)\n",     mean(@d3), median(@d3);
if ($inconsistent) {
    printf "\n%d of the inexact features are flagged inconsistent_copies: the source record\n" .
           "annotates the two IR copies of that gene at different lengths, so no single\n" .
           "reference can satisfy both. Not an annotator error.\n", $inconsistent;
}

# A histogram makes a systematic directional bias obvious at a glance.
for my $set ([\@d5, "5' offset"], [\@d3, "3' offset"]) {
    my ($v, $t) = @$set;
    my %h;
    for my $x (@$v) {
        my $b = $x == 0        ? "     0"
              : abs($x) <= 3   ? ($x < 0 ? "  -1..-3" : "   1..3")
              : abs($x) <= 15  ? ($x < 0 ? " -4..-15" : "  4..15")
              : abs($x) <= 100 ? ($x < 0 ? "-16..-100" : " 16..100")
              :                  ($x < 0 ? "   < -100" : "    > 100");
        $h{$b}++;
    }
    print "\n$t distribution:\n";
    for my $b ("   < -100","-16..-100"," -4..-15","  -1..-3","     0","   1..3","  4..15"," 16..100","    > 100") {
        next unless $h{$b};
        printf "  %-10s %4d  %s\n", $b, $h{$b}, "#" x int(40 * $h{$b} / ($nm||1));
    }
}

if (@missing) {
    print "\nmissing (first 25): ", join(", ", map { "$_->{name}\[$_->{start}\]" } @missing[0 .. ($#missing > 24 ? 24 : $#missing)]), "\n";
}
if (@spurious) {
    print "\nspurious (first 25): ", join(", ", map { "$_->{name}\[$_->{start}\]" } @spurious[0 .. ($#spurious > 24 ? 24 : $#spurious)]), "\n";
}

if ($detail_file) {
    open my $dh, ">", $detail_file or die "Cannot write $detail_file: $!";
    print $dh join("\t", qw(gene type strand truth_start truth_end pred_start pred_end pred_strand off5 off3 flags)), "\n";
    print $dh join("\t", @$_), "\n" for sort { abs($b->[8]) + abs($b->[9]) <=> abs($a->[8]) + abs($a->[9]) } @detail;
    print $dh join("\t", $_->{name}, $_->{type}, $_->{strand}, $_->{start}, $_->{end}, "", "", "", "MISSING", "", $_->{flags}), "\n" for @missing;
    print $dh join("\t", $_->{name}, "NA", "", "", "", $_->{start}, $_->{end}, $_->{strand}, "SPURIOUS", "", ""), "\n" for @spurious;
    close $dh;
    print "\ndetail written to $detail_file (worst offsets first)\n";
}
