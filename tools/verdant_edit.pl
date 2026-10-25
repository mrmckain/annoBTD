#!/usr/bin/perl
# Apply hand curation to an annoBTD annotation table and regenerate its spacer rows.
#
# USAGE: verdant_edit.pl <annotation.txt> <genome.fsa> <edits.tsv> [--out FILE]
#   edits.tsv rows (tab separated; a first column naming the genome is ignored):
#     rename  <name>  <start> <end>  <new_name>        e.g. rps16_exon1 4902 4941 rps16_pseudogene
#     add     <name>  <start> <end>  <strand>          a new feature row
#     delete  <name>  <start> <end>
#     move    <name>  <start> <end>  <new_start> <new_end>
#   A row is matched by name and both coordinates. Spacer rows (geneA~geneB) are
#   rebuilt from the surviving features; LSC/SSC/IRA/IRB/FULL rows are kept.
#   Pseudogenes have no special syntax in this format: mark them with the suffix
#   _pseudogene, which every annoBTD check skips (is_gene treats it as a gene name).
use strict; use warnings;
my ($ann, $fsa, $edits, @rest) = @ARGV; die "usage: $0 <annotation.txt> <genome.fsa> <edits.tsv> [--out FILE]\n" unless $edits;
my $out = $ann; while (@rest) { my $a = shift @rest; $out = shift @rest if $a eq '--out' }
my $len = 0; if (open my $f, "<", $fsa) { while (<$f>) { next if /^>/; s/\s//g; $len += length } close $f }
my (@rows, @regions); open my $h, "<", $ann or die "$ann: $!";
while (<$h>) { chomp; next unless /\S/; my @f = split /\t/; next if $f[0] =~ /~/; if ($f[0] =~ /^(LSC|SSC|IRA|IRB|FULL)$/) { push @regions, \@f } else { push @rows, \@f } } close $h;
my ($applied, $missed) = (0, 0); open my $e, "<", $edits or die "$edits: $!";
while (<$e>) { chomp; next if /^#/ || !/\S/; my @f = split /\t/; shift @f unless $f[0] =~ /^(rename|add|delete|move)$/;   # optional leading genome column
    my ($act, $name, $s, $e2, @x) = @f;
    if ($act eq 'add') { if (grep { $_->[0] eq $name && $_->[1] == $s && $_->[2] == $e2 && $_->[3] eq $x[0] } @rows) { $missed++; next }   # already there: an add is idempotent
        push @rows, [$name, $s, $e2, $x[0]]; $applied++; next }
    my ($i) = grep { $rows[$_][0] eq $name && $rows[$_][1] == $s && $rows[$_][2] == $e2 } 0 .. $#rows;
    unless (defined $i) { warn "no row $name $s-$e2 in $ann\n"; $missed++; next }
    if ($act eq 'rename') { $rows[$i][0] = $x[0] } elsif ($act eq 'delete') { splice @rows, $i, 1 } elsif ($act eq 'move') { ($rows[$i][1], $rows[$i][2]) = @x[0, 1] }
    $applied++ } close $e;
my @sorted = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @rows; my @gaps; my ($pn, $pe) = ('Start', 0);
for my $r (@sorted) { push @gaps, ["$pn~$r->[0]", $pe + 1, $r->[1] - 1, '+'] if $r->[1] > $pe + 1; if ($r->[2] > $pe) { $pe = $r->[2]; $pn = $r->[0] } }
push @gaps, ["$pn~End", $pe + 1, $len, '+'] if $len && $len > $pe;
my @merged = sort { $a->[1] <=> $b->[1] || ($a->[0] =~ /~/ ? 0 : 1) <=> ($b->[0] =~ /~/ ? 0 : 1) } (@sorted, @gaps);
open my $w, ">", $out or die; print $w join("\t", @$_), "\n" for @merged, @regions; close $w;
printf STDERR "%s: %d edits applied, %d not matched -> %s\n", $ann, $applied, $missed, $out;
