#!/usr/bin/perl
use strict;
use warnings;

# Splice-junction check on a finished annotation.
#
# USAGE: check_splice.pl <plastome.fsa> <annotation.txt> [--tsv out.tsv] [--fix out.txt]
#
# check_annotation.pl translates what was called and catches frames, starts, stops
# and lengths. It cannot see the commonest real error in curated plastomes: an
# intron drawn a few nucleotides away from where it actually is. When the displaced
# junction happens to encode the same protein - which it usually does, because the
# shift is chosen to keep the reading frame - every protein-level check passes and
# the wrong boundary is published, then propagated to everything annotated from it.
#
# The signal that survives is the splice site itself. Plastid introns are group II
# and open with GT - not the spliceosomal GT..AG. Measured over 44 plastomes, the
# genes whose junction is unambiguous agree: ndhB 43/43 GT..AC, petD 44/44 GT..AT,
# rpl16 44/44 GT..AC. So a non-GT donor is worth a look.
#
# Only the donor is used as the test. The acceptor was measured too and is not
# invariant enough to gate on: AC and AT dominate, but rpl2 reads GT..AA in every
# record annotated with a GT donor - tobacco, five Asparagaceae, and the grass
# records from three independent labs - and clpP intron 2 reads GT..TC. Requiring
# AY flagged all of those as errors, which they are not. The acceptor is reported
# so it can be read, and it breaks ties between candidates, but it never convicts.
#
# It is not automatically an error, and this tool does not claim it is. Three
# outcomes are separated:
#
#   RELOCATE  a different boundary is canonical AND encodes the same protein.
#             Nothing downstream changes except the intron; this is a safe fix.
#   REVIEW    a canonical boundary exists but encodes a different protein. Only a
#             human (or transcript data) can say which is right.
#   edit?     the donor is GC. C->U editing at that position yields GU, restoring
#             the canonical site in the transcript though not in the genome. This
#             is a biological finding, not a mistake, and must not be "fixed".
#   ok        non-canonical with no better alternative - accepted as real.

my ($plastome_file, $ann_file, @rest) = @ARGV;
die "usage: $0 <plastome.fsa> <annotation.txt> [--tsv out.tsv] [--fix out.txt]\n"
    unless defined $ann_file;

my ($tsv_file, $fix_file, $refprot_file, %review_ok);
while (@rest) {
    my $a = shift @rest;
    $tsv_file = shift @rest if $a eq '--tsv';
    $fix_file = shift @rest if $a eq '--fix';
    # --fix-review names the genes whose PROTEIN-CHANGING canonical alternative
    # may also be written. It is deliberately not a blanket switch: applying one
    # of these rewrites the protein, so the gene has to be named by someone who
    # has looked at it.
    if ($a eq '--fix-review') { $review_ok{$_} = 1 for split /,/, (shift(@rest) // '') }
    # --ref-prot supplies a trusted protein per gene (FASTA, header = gene name).
    # Without it a REVIEW gene can only be corrected by SLIDING the junction, which
    # moves both boundaries by the same amount and therefore cannot change the
    # coding length. That is the wrong shape of correction for a gene whose true
    # intron is a different size: ycf3 exon1 132->124 goes with exon2 228->230, a
    # net loss of two codons, so no slide can reach it. With a reference the two
    # boundaries are searched independently and the structure whose protein best
    # matches the reference wins.
    $refprot_file = shift @rest if $a eq '--ref-prot';
}

my %CODON = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

my $SLIDE = 15;      # how far a junction may be slid while looking for GT..AY
my $ROAM  = 1500;    # how far a short exon may be relocated
my $SHORT = 15;      # an exon this small can sit almost anywhere by chance

my %refprot;
if ($refprot_file) {
    my $g;
    open my $rh, "<", $refprot_file or die "Cannot open $refprot_file: $!";
    # One gene may appear several times, one record per trusted source. They are
    # kept apart: concatenating them would build a chimera that matches nothing.
    while (<$rh>) {
        chomp;
        if (/^>\s*(\S+)/) { ($g = $1) =~ s/_exon\d+$//; push @{ $refprot{$g} }, '' }
        elsif (defined $g) { $refprot{$g}[-1] .= $_ }
    }
    close $rh;
}
# Alignment-free similarity, so a genuine indel between reference and target does
# not throw the comparison off the way a positional match would.
#
# SYMMETRIC (Dice over 4-mer sets), which matters here. One-directional
# containment scores a protein that is a strict subset of the reference as a
# perfect match, so a candidate missing a residue at the junction beats the
# candidate that has it - which is exactly the wrong way round when the question
# being asked is where the junction goes.
sub kmer_sim {
    my ($x, $y) = @_;
    return 0 unless length($x) >= 4 && length($y) >= 4;
    my (%X, %Y);
    $X{ substr($x, $_, 4) } = 1 for 0 .. length($x) - 4;
    $Y{ substr($y, $_, 4) } = 1 for 0 .. length($y) - 4;
    my $shared = grep { $Y{$_} } keys %X;
    my $tot = keys(%X) + keys(%Y);
    return $tot ? 2 * $shared / $tot : 0;
}

my $seq = '';
open my $pf, "<", $plastome_file or die "Cannot open $plastome_file: $!";
while (<$pf>) { chomp; next if /^>/; s/\s+//g; $seq .= $_ }
close $pf;
die "No sequence in $plastome_file\n" unless length $seq;
$seq = uc $seq;

sub rc { my $x = reverse shift; $x =~ tr/ACGT/TGCA/; return $x }
sub translate {
    my $s = shift;
    my $p = '';
    for (my $i = 0; $i + 2 < length $s; $i += 3) { $p .= $CODON{ substr($s, $i, 3) } // 'X' }
    $p =~ s/\*$//;
    return $p;
}
# Exons are concatenated in ASCENDING coordinate order and the whole strand
# reverse complemented once. Reversing the exon list as well would undo itself.
sub cds_of {
    my ($ex, $dir) = @_;
    my @a = sort { $a->{s} <=> $b->{s} } @$ex;
    my $c = '';
    $c .= substr($seq, $_->{s} - 1, $_->{e} - $_->{s} + 1) for @a;
    return $dir eq '-' ? rc($c) : $c;
}

my (%parts, @order);
open my $af, "<", $ann_file or die "Cannot open $ann_file: $!";
while (<$af>) {
    chomp;
    next unless /\S/;
    my ($name, $s, $e, $d) = split /\t/;
    next unless defined $d && defined $s && $s =~ /^\d+$/ && $e =~ /^\d+$/;
    next unless $d eq '+' || $d eq '-';
    next if $name =~ /^(LSC|SSC|IRA|IRB|FULL)$/;
    next if $name =~ /~/ || $name =~ /_intron\d*$/;
    next if $name =~ /^trn/ || $name =~ /^rrn/;
    next if $name =~ /^rps12(_|$)/;          # trans-spliced; its exons are not one locus
    next unless $name =~ /^(.+)_exon(\d+)$/;  # only spliced genes have junctions
    my ($gene, $n) = ($1, $2);
    ($s, $e) = ($e, $s) if $e < $s;
    push @{ $parts{"$gene\t$d"} }, { s => $s, e => $e, n => $n };
    push @order, "$gene\t$d" unless grep { $_ eq "$gene\t$d" } @order;
}
close $af;

my (@report, %fix);
my %tally;
for my $key (@order) {
    my ($gene, $dir) = split /\t/, $key;
    my @ex = sort { $a->{s} <=> $b->{s} } @{ $parts{$key} };

    # IR copies of the same gene are annotated twice; split on a coordinate gap.
    my @copies = ([ shift @ex ]);
    for my $x (@ex) {
        if ($x->{s} - $copies[-1][-1]{e} > 5000) { push @copies, [$x] }
        else { push @{ $copies[-1] }, $x }
    }

    for my $copy (@copies) {
        next if @$copy < 2;
        # Transcription order: descending coordinates on the minus strand.
        my @o = sort { $a->{s} <=> $b->{s} } @$copy;
        @o = reverse @o if $dir eq '-';
        my $base = translate(cds_of(\@o, $dir));

        # Donor and acceptor of intron $i, with the junction slid by $k nt.
        my $sig = sub {
            my ($i, $k) = @_;
            my ($a, $b) = ($o[$i], $o[$i+1]);
            return ($dir eq '+')
                ? ( substr($seq, $a->{e} + $k, 2), substr($seq, $b->{s} - 3 + $k, 2) )
                : ( rc(substr($seq, $a->{s} - 3 - $k, 2)), rc(substr($seq, $b->{e} - $k, 2)) );
        };

        for my $i (0 .. $#o - 1) {
            my ($don, $acc) = $sig->($i, 0);
            my $intron = $dir eq '+' ? $o[$i+1]{s} - $o[$i]{e} - 1
                                     : $o[$i]{s} - $o[$i+1]{e} - 1;
            $tally{"$don..$acc"}++;
            next if $don eq 'GT';

            my (@same, @diff);

            # (1) Slide the junction. Both sides move together so the coding
            # length is preserved, but the sequence is not - it is checked.
            for my $k (-$SLIDE .. $SLIDE) {
                next unless $k;
                my ($d2, $a2) = $sig->($i, $k);
                next unless $d2 eq 'GT';
                my @n = map { { %$_ } } @o;
                if ($dir eq '+') { $n[$i]{e} += $k; $n[$i+1]{s} += $k }
                # Minus strand: transcription runs toward LOWER coordinates, so
                # sliding the junction downstream moves BOTH boundaries down. The
                # acceptor must travel the same way as the donor or the intron is
                # shortened from both ends instead of slid, and the two IR copies
                # of one gene then disagree with each other.
                else             { $n[$i]{s} -= $k; $n[$i+1]{e} -= $k }
                next if $n[$i]{e} < $n[$i]{s} || $n[$i+1]{e} < $n[$i+1]{s};
                my $p = translate(cds_of(\@n, $dir));
                my $rec = { k => $k, don => $d2, acc => $a2, ex => \@n,
                            what => sprintf("slide %+d nt", $k) };
                push @{ $p eq $base ? \@same : \@diff }, $rec;
            }

            # (2) Relocate a short exon. petB and petD open with 6-8 nt, which
            # recurs by chance every few hundred bases, so the wrong copy can be
            # picked without changing the protein at all - only the intron.
            my $elen = $o[$i]{e} - $o[$i]{s} + 1;
            if ($i == 0 && $elen <= $SHORT) {
                for my $off (-$ROAM .. $ROAM) {
                    next unless $off;
                    my @n = map { { %$_ } } @o;
                    $n[0]{s} += $off; $n[0]{e} += $off;
                    next if $n[0]{s} < 1 || $n[0]{e} > length $seq;
                    # must still lie outside exon 2, on the correct side
                    next if $dir eq '+' ? $n[0]{e} >= $o[1]{s} : $n[0]{s} <= $o[1]{e};
                    my ($d2, $a2) = ($dir eq '+')
                        ? ( substr($seq, $n[0]{e}, 2), $acc )
                        : ( rc(substr($seq, $n[0]{s} - 3, 2)), $acc );
                    next unless $d2 eq 'GT';
                    my $p = translate(cds_of(\@n, $dir));
                    next unless $p eq $base;
                    push @same, { k => $off, don => $d2, acc => $a2, ex => \@n,
                                  what => sprintf("move exon1 %+d nt", $off) };
                }
            }

            # Search donor and acceptor independently against a trusted protein.
            # Reachable from any non-canonical case, not only from the @diff one:
            # a GC donor is normally left alone, but once the gene is named in
            # --fix-review the reference is the better authority.
            my $ref_search = sub {
                my $refs = $refprot{$gene} or return undef;
                my $best;
                for my $kd (-$SLIDE .. $SLIDE) {
                    for my $ka (-$SLIDE .. $SLIDE) {
                        my @n = map { { %$_ } } @o;
                        if ($dir eq '+') { $n[$i]{e} += $kd; $n[$i+1]{s} += $ka }
                        else             { $n[$i]{s} -= $kd; $n[$i+1]{e} -= $ka }
                        next if $n[$i]{e} < $n[$i]{s} || $n[$i+1]{e} < $n[$i+1]{s};
                        my $cds = cds_of(\@n, $dir);
                        next if length($cds) % 3;
                        my ($d2, $a2) = ($dir eq '+')
                            ? ( substr($seq, $n[$i]{e}, 2), substr($seq, $n[$i+1]{s} - 3, 2) )
                            : ( rc(substr($seq, $n[$i]{s} - 3, 2)), rc(substr($seq, $n[$i+1]{e}, 2)) );
                        next unless $d2 eq 'GT';
                        my $p = translate($cds);
                        next if $p =~ /\*./;            # internal stop
                        # Best agreement with any one reference, never the average:
                        # a single divergent source must not veto.
                        my $sc = 0;
                        for my $r (@$refs) { my $x = kmer_sim($r, $p); $sc = $x if $x > $sc }
                        next if $best && $sc <= $best->{sc};
                        $best = { sc => $sc, ex => \@n, don => $d2, acc => $a2,
                                  kd => $kd, ka => $ka, len => length $p };
                    }
                }
                return $best;
            };
            my $applied = sub {
                my $c = shift;
                $fix{"$gene\t$copy->[0]{s}"} = $c->{ex};
                return sprintf("; against the reference: donor %+d acceptor %+d gives %s..%s, %d aa, %.0f%% match [applied]",
                               $c->{kd}, $c->{ka}, $c->{don}, $c->{acc}, $c->{len}, 100 * $c->{sc});
            };

            my ($verdict, $note);
            if (@same) {
                # Prefer the alternative whose intron length is most ordinary.
                @same = sort { ($b->{acc} =~ /^A/) <=> ($a->{acc} =~ /^A/)
                            || abs($a->{k}) <=> abs($b->{k}) } @same;
                my $c = $same[0];
                my $nint = $dir eq '+' ? $c->{ex}[$i+1]{s} - $c->{ex}[$i]{e} - 1
                                       : $c->{ex}[$i]{s} - $c->{ex}[$i+1]{e} - 1;
                $verdict = 'RELOCATE';
                $note = sprintf("%s gives %s..%s, intron %d->%d nt, protein unchanged",
                                $c->{what}, $c->{don}, $c->{acc}, $intron, $nint);
                $fix{"$gene\t$copy->[0]{s}"} = $c->{ex};
            }
            elsif ($don eq 'GC') {
                $verdict = 'edit?';
                $note = "GC donor: C->U editing yields GU, canonical in the transcript"
                      . (@diff ? sprintf("; a GT donor %+d nt away would change the protein", $diff[0]{k}) : "");
                if ($review_ok{$gene}) {
                    my $cand = $ref_search->();
                    if ($cand) { $verdict = 'REVIEW'; $note .= $applied->($cand) }
                }
            }
            elsif (@diff) {
                @diff = sort { ($b->{acc} =~ /^A/) <=> ($a->{acc} =~ /^A/)
                            || abs($a->{k}) <=> abs($b->{k}) } @diff;
                $verdict = 'REVIEW';
                $note = sprintf("%s gives %s..%s but changes the protein",
                                $diff[0]{what}, $diff[0]{don}, $diff[0]{acc});
                if ($review_ok{$gene}) {
                    if ($refprot{$gene}) {
                        my $cand = $ref_search->();
                        $note .= $cand ? $applied->($cand)
                                       : "; no GT structure matched the reference, left alone";
                    }
                    else {
                        $fix{"$gene\t$copy->[0]{s}"} = $diff[0]{ex};
                        $note .= " [slide applied on request]";
                    }
                }
            }
            else {
                $verdict = 'ok';
                $note = "no GT donor reachable within $SLIDE nt; accepted as real";
            }
            push @report, [ $verdict, $gene, $i + 1, $dir, "$don..$acc", $intron, $note ];
        }
    }
}

my %rank = (RELOCATE => 0, REVIEW => 1, 'edit?' => 2, ok => 3);
my $act = grep { $_->[0] eq 'RELOCATE' || $_->[0] eq 'REVIEW' } @report;
printf "=== %s: %d junctions, %d non-GT donor, %d actionable ===\n",
    $ann_file, scalar(grep { $_ } map { $tally{$_} } keys %tally) && do { my $t=0; $t+=$_ for values %tally; $t },
    scalar(@report), $act;
for my $r (sort { $rank{$a->[0]} <=> $rank{$b->[0]} || $a->[1] cmp $b->[1] } @report) {
    printf "  %-9s %-10s intron%d %s  %-6s %5d nt  %s\n", @$r[0,1,2,3,4,5,6];
}

if ($tsv_file) {
    open my $th, ">", $tsv_file or die "Cannot write $tsv_file: $!";
    print $th join("\t", qw(verdict gene intron strand signal intron_len note)), "\n";
    print $th join("\t", @$_), "\n" for @report;
    close $th;
}
if ($fix_file) {
    open my $fh, ">", $fix_file or die "Cannot write $fix_file: $!";
    for my $k (sort keys %fix) {
        my ($g) = split /\t/, $k;
        my @e = sort { $a->{s} <=> $b->{s} } @{ $fix{$k} };
        printf $fh "%s_exon%d\t%d\t%d\n", $g, $_->{n}, $_->{s}, $_->{e} for @e;
    }
    close $fh;
    printf STDERR "wrote %s: %d gene copies with canonical boundaries\n",
        $fix_file, scalar keys %fix;
}
exit($act > 255 ? 255 : $act);
