#!/usr/bin/perl
# Apply the lineages check_taxonomy_consistency.pl proposed (verdict ASSIGNABLE) to
# records that had none: family and order in the taxonomy table and in the sketch
# database. Pre-assignment copies are kept as .pre_assign; a side table records every assignment
# with its evidence so it can be audited or reverted.
#
# USAGE: assign_lineages.pl <taxonomy_consistency.tsv> <taxonomy_all.tsv> <sketch_db.tsv> <assignments_out.tsv>
use strict; use warnings;
my ($chk, $tax, $db, $out) = @ARGV; die "usage: $0 <taxonomy_consistency.tsv> <taxonomy_all.tsv> <sketch_db.tsv> <assignments_out.tsv>\n" unless $out;
my %as; open my $c, "<", $chk or die; <$c>;
while (<$c>) { chomp; my @f = split /\t/; next unless $f[10] eq 'ASSIGNABLE'; $as{$f[0]} = [$f[4], $f[5], $f[7], $f[9], $f[1]] }
close $c; printf STDERR "%d assignable records\n", scalar keys %as;
open my $o, ">", $out or die; print $o join("\t", qw(accession organism family order c_best j_best source)), "\n";
print $o join("\t", $_, $as{$_}[4], $as{$_}[0], $as{$_}[1], $as{$_}[2], $as{$_}[3], 'sketch'), "\n" for sort keys %as; close $o;
for my $spec ([$tax, 5, 6], [$db, 2, 3]) { my ($file, $fi, $oi) = @$spec;
    # keep a copy of the file as it was before this assignment (a .orig may already
    # exist from an earlier, unrelated repair and must not be read back)
    my $bak = "$file.pre_assign"; if (-e $bak) { die "$bak exists: assignments were already applied; remove it to reapply\n" }
    rename $file, $bak or die; open my $in, "<", $bak or die; open my $w, ">", $file or die; my $n = 0;
    while (<$in>) { chomp; my @f = split /\t/, $_, -1; if ($as{$f[0]} && $f[$fi] eq 'NA') { $f[$fi] = $as{$f[0]}[0]; $f[$oi] = $as{$f[0]}[1]; $n++ } print $w join("\t", @f), "\n" }
    close $in; close $w; printf STDERR "%s: %d rows assigned\n", $file, $n }
