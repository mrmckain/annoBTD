#!/usr/bin/perl
use strict; use warnings;
# Flag gene x family cells that do not follow the evolutionary pattern of the rest:
#   LENGTH  a family's median length far from other families' medians (robust z)
#   LOSS    a gene present in nearly all records overall but absent from most of
#           a family's records - the parasitic-lineage signature
# USAGE: flag_profile_outliers.pl <gene_expect.tsv> <cds_rows.tsv>
my ($exp, $rows) = @ARGV; die "usage\n" unless $rows;
my (%fam_med, %all_med, %fam_n, %recs, %gene_fam_recs, %gene_recs, $nrec);
open my $r, "<", $rows or die; my %seen;
while (<$r>) { chomp; my @f = split /\t/; next if $. == 1; my ($acc,$fam,$gene) = @f[0,1,3];
  $recs{$fam}{$acc} = 1; $gene_fam_recs{$gene}{$fam}{$acc} = 1; $gene_recs{$gene}{$acc} = 1; $seen{$acc} = 1 }
close $r; $nrec = keys %seen;
open my $e, "<", $exp or die;
while (<$e>) { chomp; my ($n,$lt,$lin,$cnt,$med) = split /\t/; next if $n eq 'name';
  if ($lt eq 'F' && $cnt >= 5) { $fam_med{$n}{$lin} = $med; $fam_n{$n}{$lin} = $cnt } elsif ($lt eq 'ALL') { $all_med{$n} = $med } }
close $e;
sub med { my @v = sort { $a <=> $b } @_; @v ? $v[@v/2] : undef }
print "# LENGTH outliers: family median vs other families (robust z >= 3 and >= 10% off)\n";
print join("\t", qw(gene family n family_median other_families_median ratio z)), "\n";
my @len;
for my $g (sort keys %fam_med) { my @fams = keys %{ $fam_med{$g} }; next if @fams < 6;
  for my $f (@fams) { my @others = map { $fam_med{$g}{$_} } grep { $_ ne $f } @fams; my $m = med(@others);
    my $mad = med(map { abs($_ - $m) } @others) || 1; my $z = abs($fam_med{$g}{$f} - $m) / (1.4826 * $mad);
    my $ratio = $fam_med{$g}{$f} / ($m || 1);
    push @len, [ $g, $f, $fam_n{$g}{$f}, $fam_med{$g}{$f}, $m, $ratio, $z ] if $z >= 3 && abs($ratio - 1) >= 0.10 } }
printf "%s\t%s\t%d\t%d\t%d\t%.2f\t%.1f\n", @$_ for sort { $b->[2] <=> $a->[2] } @len;
print "\n# LOSS candidates: gene in >=90% of all records but <=50% of a family's (family >=10 records)\n";
print join("\t", qw(gene family family_records with_gene fraction overall_fraction)), "\n";
my @loss;
for my $g (sort keys %gene_recs) { my $of = keys(%{ $gene_recs{$g} }) / $nrec; next if $of < 0.90;
  for my $f (sort keys %recs) { my $nf = keys %{ $recs{$f} }; next if $nf < 10; my $w = keys %{ $gene_fam_recs{$g}{$f} || {} };
    push @loss, [ $g, $f, $nf, $w, $w / $nf, $of ] if $w / $nf <= 0.50 } }
printf "%s\t%s\t%d\t%d\t%.2f\t%.2f\n", @$_ for sort { $a->[4] <=> $b->[4] || $b->[2] <=> $a->[2] } @loss;

# RECORD-level: a CDS instance that disagrees with its own family's consensus.
# This is the unit that "does not follow the evolutionary pattern": a whole family
# agreeing on an unusual length (infA 107 aa across 470 Poaceae records) IS the
# pattern and must be kept; one record 30% off its family's median is what a
# guide-quality filter should discount.
print "\n# RECORD outliers: instance length >=15% off its family median (family n>=10)\n";
print join("\t", qw(accession gene family prot_len family_median ratio)), "\n";
my (%fl, %fn);
open $r, "<", $rows or die; my %inst;
while (<$r>) { chomp; my @f = split /\t/; next if $. == 1; push @{ $inst{"$f[3]\t$f[1]"} }, [ $f[0], $f[6] ] }
close $r;
my ($nout, $nall, %bad_acc);
for my $k (sort keys %inst) { my @i = @{ $inst{$k} }; next if @i < 10; my $m = med(map { $_->[1] } @i);
  # Only a settled cell can convict: if GenBank itself is split on this gene in
  # this family, a record on the minority side is following a convention, not
  # misannotated. Require >=90% of the family's records within 2% of the median.
  my $settled = (grep { abs($_->[1] - $m) / ($m || 1) <= 0.02 } @i) / @i; next if $settled < 0.90;
  for my $x (@i) { $nall++; my $ratio = $x->[1] / ($m || 1); next unless abs($ratio - 1) >= 0.15;
    my ($g, $f) = split /\t/, $k; printf "%s\t%s\t%s\t%d\t%d\t%.2f\n", $x->[0], $g, $f, $x->[1], $m, $ratio; $nout++; $bad_acc{ $x->[0] }++ } }
printf STDERR "record outliers: %d of %d instances (%.1f%%), in %d records\n", $nout, $nall, 100 * $nout / ($nall || 1), scalar keys %bad_acc;
open my $q, ">", "guide_quality.tsv" or die; print $q "accession\toutlier_cds\n"; printf $q "%s\t%d\n", $_, $bad_acc{$_} for sort { $bad_acc{$b} <=> $bad_acc{$a} } keys %bad_acc; close $q;
