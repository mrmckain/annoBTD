#!/usr/bin/perl
use strict;
use warnings;

# Pull the headline numbers out of one compare_annotations.pl report and append
# them as a row to the cross-genome summary table.
#
# USAGE: summarize.pl <compare.txt> <summary.tsv> <genome> <seconds>

my ($report, $summary, $genome, $seconds) = @ARGV;
die "usage: $0 <compare.txt> <summary.tsv> <genome> <seconds>\n" unless defined $seconds;

open my $fh, "<", $report or die "Cannot open $report: $!";
my $t = do { local $/; <$fh> };
close $fh;

my %v;
$v{truth}  = $1 if $t =~ /truth features\s*:\s*(\d+)/;
$v{pred}   = $1 if $t =~ /predicted features\s*:\s*(\d+)/;
($v{match}, $v{recall}) = ($1, $2) if $t =~ /matched\s*:\s*(\d+)\s*\(recall\s*([\d.]+)/;
$v{both}   = $1 if $t =~ /both boundaries exact\s*:.*?\(([\d.]+)%\)/;
$v{p5}     = $1 if $t =~ /5. exact\s*:.*?\(([\d.]+)%\)/;
$v{p3}     = $1 if $t =~ /3. exact\s*:.*?\(([\d.]+)%\)/;
$v{m5}     = $1 if $t =~ /5. offset\s+mean\s+\S+\s+median\s+([+-]?\d+)/;
$v{m3}     = $1 if $t =~ /3. offset\s+mean\s+\S+\s+median\s+([+-]?\d+)/;
$v{miss}   = $1 if $t =~ /^missing\s*:\s*(\d+)/m;
$v{spur}   = $1 if $t =~ /^spurious\s*:\s*(\d+)/m;

open my $out, ">>", $summary or die "Cannot append to $summary: $!";
print $out join("\t", $genome,
    map { defined $v{$_} ? $v{$_} : "" }
        qw(truth pred match recall both p5 p3 m5 m3 miss spur)),
    "\t$seconds\n";
close $out;
