#!/usr/bin/env perl
# check_homopolymers.pl: count each homopolymer run that annoBTD reported in
# the sequencing reads, and say whether the reads back the assembly. The report
# lists the runs where a gene's reading frame shifts and every run of 10 or more
# identical bases inside a gene that is in frame as assembled.
#
#   check_homopolymers.pl <assembly.fasta> <prefix>_homopolymers.tsv \
#        --reads R1.fastq.gz [R2.fastq.gz ...] [--out FILE] [--fix FASTA]
#
# For every run in the report it finds the reads that cross the run with
# --flank (default 15) exactly matching bases on both sides, on either strand,
# and tabulates the run length each of them shows.
#
# Verdicts
#   confirmed       the assembly's length is the commonest in the reads, by
#                   at least 1.5 times the next length
#   reads differ    another length is at least 1.5 times as common as the
#                   assembly's; that length is suggested, whether or not it
#                   keeps the gene in frame (the note says which)
#   unresolved      two or more lengths are about equally common. If one of
#                   them keeps the gene in frame it is suggested, because an
#                   intact gene is the better guess when the reads cannot
#                   decide
#   too few reads   fewer than --min-reads (default 10) reads cross the run
#
# Reads undercount long runs more often than they overcount them (PCR and
# sequencing slippage), so the shorter neighbour of the true length is always
# well represented. Above about 18 bases the commonest length is itself
# uncertain.
#
# --fix FASTA writes the assembly with every suggested length applied; when no
# length is suggested the file is not written. The coordinates of everything
# downstream of a changed run shift, so annotate the fixed sequence again.
use strict;
use warnings;
use Getopt::Long;

my ($out, $fix, @reads);
my ($flank, $minreads, $maxrun) = (15, 10, 80);
GetOptions('reads=s{1,}' => \@reads, 'out=s' => \$out, 'fix=s' => \$fix,
           'flank=i' => \$flank, 'min-reads=i' => \$minreads)
    or die "check_homopolymers.pl: bad option\n";
my ($fa, $report) = @ARGV;
die "usage: check_homopolymers.pl <assembly.fasta> <prefix>_homopolymers.tsv --reads R1.fastq.gz [R2.fastq.gz ...] [--out FILE] [--fix FASTA]\n"
    unless $fa && $report && @reads;
for ($fa, $report, @reads) { -s $_ or die "check_homopolymers.pl: $_ is missing or empty\n" }

open my $F, '<', $fa or die "$fa: $!\n";
my ($head, $g) = ('', '');
while (<$F>) {
    if (/^>/) { die "check_homopolymers.pl: $fa has more than one sequence\n" if $head; $head = $_; next }
    s/\s//g; $g .= uc;
}
close $F;

sub rc { my $s = reverse $_[0]; $s =~ tr/ACGT/TGCA/; $s }

# The report's runs, each widened to the whole run in the assembly.
my @S; my $fcol;
open my $R, '<', $report or die "$report: $!\n";
while (<$R>) {
    s/\r?\n\z//;
    my @c = split /\t/;
    if (@c && $c[0] eq 'gene') { ($fcol) = grep { $c[$_] eq 'gene_frame' } 0 .. $#c; next }
    next if !@c || @c < 9;
    my ($s, $e) = @c[4, 5];
    my $b = substr($g, $s - 1, 1);
    $s-- while $s > 1 && substr($g, $s - 2, 1) eq $b;
    $e++ while $e < length($g) && substr($g, $e, 1) eq $b;
    next if $s - 1 < $flank || $e + $flank > length $g;
    my $L = substr($g, $s - 1 - $flank, $flank);
    my $Rf = substr($g, $e, $flank);
    my ($rL, $rR, $rb) = (rc($Rf), rc($L), rc($b));
    push @S, { gene => $c[0], s => $s, e => $e, b => $b, rb => $rb, len => $e - $s + 1,
               inframe => { map { $_ => 1 } $c[8] =~ /(\d+)/g }, inframe_text => $c[8], lf => $L, rlf => $rL,
               frame => (defined $fcol && defined $c[$fcol] && length $c[$fcol] ? $c[$fcol] : 'shifted at this run'),
               f => qr/$L([ACGTN]{1,$maxrun}?)$Rf/, r => qr/$rL([ACGTN]{1,$maxrun}?)$rR/,
               n => {}, other => 0 };
}
close $R;

my $nreads = 0;
# one pass decides whether a read can cross any run at all: the left flank of some run, on either strand
my $any; if (@S) { my $alt = join '|', map { ($_->{lf}, $_->{rlf}) } @S; $any = qr/$alt/ }
if (@S) {
    for my $file (@reads) {
        my $fh;
        if ($file =~ /\.gz$/) { open $fh, '-|', 'gzip', '-dc', $file or die "$file: $!\n" }
        else                  { open $fh, '<', $file or die "$file: $!\n" }
        while (defined(my $h = <$fh>)) {
            my $seq = <$fh>; my $plus = <$fh>; my $qual = <$fh>;
            last unless defined $qual;
            $nreads++;
            next unless $seq =~ $any;
            for my $x (@S) {
                if ($seq =~ $x->{f}) {
                    my $m = $1;
                    if ($m =~ /^$x->{b}+$/) { $x->{n}{ length $m }++ } else { $x->{other}++ }
                }
                elsif ($seq =~ $x->{r}) {
                    my $m = $1;
                    if ($m =~ /^$x->{rb}+$/) { $x->{n}{ length $m }++ } else { $x->{other}++ }
                }
            }
        }
        close $fh;
    }
}

my $O;
if ($out) { open $O, '>', $out or die "$out: $!\n" } else { $O = \*STDOUT }
print $O join("\t", qw(gene run_start run_end base_forward_strand assembly_length in_frame_if_run_were
                       reads_crossing reads_with_another_base lengths_in_reads verdict suggested_length note gene_frame)), "\n";
my @edits;
for my $x (@S) {
    my %n = %{ $x->{n} };
    my $tot = 0; $tot += $_ for values %n;
    my @len = sort { $n{$b} <=> $n{$a} || $a <=> $b } keys %n;
    my $dist = join(' ', map { "$_:$n{$_}" } sort { $a <=> $b } keys %n) || 'none';
    my $asm = $n{ $x->{len} } || 0;
    my ($verdict, $suggest, $note) = ('', '', '');
    if ($tot < $minreads) {
        $verdict = 'too few reads';
        $note = "only $tot reads cross the run";
    }
    else {
        my $top = $n{ $len[0] };
        # Lengths the reads cannot tell apart from the commonest one.
        my @close = sort { $a <=> $b } grep { $n{$_} * 1.5 > $top } @len;
        if (@close == 1 && $close[0] == $x->{len}) {
            $verdict = 'confirmed';
            $note = sprintf('%d of %d reads (%.0f%%) show %d', $asm, $tot, 100 * $asm / $tot, $x->{len});
        }
        elsif (@close == 1) {
            $verdict = 'reads differ';
            $suggest = $close[0];
            $note = sprintf('%d of %d reads (%.0f%%) show %d; %d show the assembly\'s %d; %d %s',
                            $top, $tot, 100 * $top / $tot, $close[0], $asm, $x->{len}, $close[0],
                            $x->{inframe}{ $close[0] } ? 'keeps the gene in frame' : 'puts the gene out of frame');
        }
        else {
            $verdict = 'unresolved';
            my @in = grep { $x->{inframe}{$_} } @close;
            $note = 'reads split between ' . join(' and ', map { "$_ ($n{$_})" } @close);
            if (@in) {
                ($suggest) = sort { $n{$b} <=> $n{$a} } @in;
                $note .= "; $suggest keeps the gene in frame";
            }
            else { $note .= '; none of them keeps the gene in frame' }
        }
    }
    $suggest = '' if $suggest ne '' && $suggest == $x->{len};
    push @edits, [$x->{s}, $x->{e}, $x->{b}, $suggest] if $suggest ne '';
    print $O join("\t", $x->{gene}, $x->{s}, $x->{e}, $x->{b}, $x->{len}, $x->{inframe_text},
                  $tot, $x->{other}, $dist, $verdict, $suggest, $note, $x->{frame}), "\n";
}
close $O if $out;

if ($fix && !@edits) {
    unlink $fix if -e $fix;
    print STDERR "check_homopolymers.pl: no run needs another length; $fix not written\n";
}
elsif ($fix) {
    my $new = $g;
    my %done;
    # Right to left, so earlier coordinates stay valid; each run once.
    for my $e (sort { $b->[0] <=> $a->[0] } @edits) {
        next if $done{ $e->[0] }++;
        substr($new, $e->[0] - 1, $e->[1] - $e->[0] + 1) = $e->[2] x $e->[3];
    }
    open my $W, '>', $fix or die "$fix: $!\n";
    print $W $head;
    print $W "$1\n" while $new =~ /(.{1,60})/g;
    close $W;
    printf STDERR "check_homopolymers.pl: %d run(s) changed, %d -> %d nt, written to %s\n",
        scalar(keys %done), length $g, length $new, $fix;
}
printf STDERR "check_homopolymers.pl: %d run(s) checked in %d reads\n", scalar @S, $nreads;
