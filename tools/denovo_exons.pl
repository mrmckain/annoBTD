#!/usr/bin/perl
use strict;
use warnings;

# EXPERIMENT: solve a multi-exon gene's internal structure de novo inside an
# anchored locus, using the gene's protein as a single unit.
#
# USAGE: denovo_exons.pl <plastome.fsa> <locus_start> <locus_end> <strand> \
#                        <reference_protein.fasta|-p PROTEIN> [--exons N]
#
# The production path matches each exon against its own per-exon reference, so a
# guide that decomposes the gene differently transfers the wrong decomposition and
# every boundary inherits whatever that guide said. This asks a different question:
# given a window that contains the whole gene and the protein it should encode,
# where must the exon boundaries be?
#
# The structure is constrained rather than copied:
#   * the joined exons must translate to something close to the reference protein
#   * introns carry a group II donor (GT) and an A-ending acceptor
#   * introns are of plausible plastid length
#   * the reading frame is continuous across the joins - the total is a multiple
#     of three, there is a start codon, a terminal stop and no internal stops
#
# Nothing about the guide's own exon boundaries is used, only its protein.

my ($fasta, $lstart, $lend, $strand, $ref, @rest) = @ARGV;
die "usage: $0 <plastome.fsa> <start> <end> <strand> <ref_protein.fa|-p SEQ> [--exons N]\n"
    unless defined $ref;

# With -p the protein is the next argument; take it before the option loop, which
# would otherwise consume it.
my $inline_prot;
$inline_prot = shift @rest if $ref eq '-p';

my $NEXON = 3;
my $VERBOSE = 0;
my $ANY_DONOR = 0;
while (@rest) {
    my $a = shift @rest;
    $NEXON     = shift @rest if $a eq '--exons';
    $VERBOSE   = 1            if $a eq '--verbose';
    $ANY_DONOR = 1            if $a eq '--any-donor';
}

my %CODON = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
# Plastid CDS do not all begin ATG. Measured over 2,248 CDS in 30 curated grass
# plastomes: ATG 96.9%, ACG 1.6% (psbC, rpl2), GTG 1.5% (psbC, rps19), GCG once.
# ACG is the RNA-editing case - C->U editing at position 2 makes the transcript
# AUG even though the genome reads ACG - so a genomic non-ATG start is normal and
# must not be treated as an error. ATA/ATT/TTG are documented plastid initiators
# that did not occur in that set but are accepted rather than risk a false flag.
# The set below is the initiator list of NCBI translation table 11, the bacterial /
# archaeal / plant plastid code these records declare: TTG CTG ATT ATC ATA ATG GTG.
# Those are initiators by definition, not anomalies. ATC was missing and cost real
# genes: Nicotiana psbI reads ATC TAT TCT and Spinacia ndhD reads ATC ACG AAT, and
# both GenBank records translate that first codon as M.
#
# ACG and GCG are kept on top of table 11 as the RNA-editing cases - C->U editing
# makes a genomic ACG read as AUG, and GCG as GUG. Anemarrhena asphodeloides
# (NC_032698.1) ndhK is a CTG example, reading CTG where Hosta and Camassia read
# ATG with the next twelve nucleotides identical.
#
# Even this list is not exhaustive, because a start can be declared rather than
# inferred: Sorghum bicolor rpl23 begins TAC and the record carries
# /transl_except=(pos:59411..59413,aa:Met). Sequence alone cannot see that, which
# is why an unrecognised start is reported but never treated as fatal.
my %START_OK = map { $_ => 1 } qw(ATG ACG GTG GCG ATA ATT ATC TTG CTG);

# petB, petD and rpl16 open with exons of 6, 8 and 9 nt. A minimum exon length of
# 20 made those structures unreachable and the solver returned nothing at all for
# every copy of those three genes.
my $MIN_EXON   = 3;
my $MIN_INTRON = 150;
my $MAX_INTRON = 3000;

my $genome = '';
open my $fh, "<", $fasta or die "Cannot open $fasta: $!";
while (<$fh>) { chomp; next if /^>/; s/\s+//g; $genome .= $_ }
close $fh;
$genome = uc $genome;

# The locus in transcription orientation, so every offset below is 5'->3'.
my $locus = substr($genome, $lstart - 1, $lend - $lstart + 1);
if ($strand eq '-') { $locus = reverse $locus; $locus =~ tr/ACGT/TGCA/ }
my $L = length $locus;

my $refprot;
if ($ref eq '-p') { $refprot = $inline_prot }
else {
    open my $rf, "<", $ref or die "Cannot open $ref: $!";
    while (<$rf>) { chomp; next if /^>/; s/\s+//g; $refprot .= $_ }
    close $rf;
}
$refprot = uc $refprot;
$refprot =~ s/\*+$//;
die "No reference protein\n" unless length $refprot;
my $target_nt = (length($refprot) + 1) * 3;      # protein plus its stop codon

sub tr_ {
    my $s = shift;
    my $p = '';
    for (my $i = 0; $i + 2 < length $s; $i += 3) {
        my $c = substr($s, $i, 3);
        $p .= exists $CODON{$c} ? $CODON{$c} : 'X';
    }
    return $p;
}
# Indel-tolerant similarity: the fraction of the reference's 4-mers present in the
# candidate. Positional identity was tried first and is wrong here - Spinacia's
# ycf3 carries three fewer N-terminal residues than tobacco's, and comparing
# position by position scored the correct structure at 0.05.
sub sim {
    my ($a, $b) = @_;
    return 0 unless length($a) >= 4 && length($b) >= 4;
    my %A;
    $A{ substr($a, $_, 4) } = 1 for 0 .. length($a) - 4;
    my ($hit, $tot) = (0, 0);
    for my $i (0 .. length($b) - 4) { $tot++; $hit++ if $A{ substr($b, $i, 4) } }
    return $tot ? $hit / $tot : 0;
}
# Positional identity, kept for the cheap prefix test where the candidate and the
# reference prefix start at the same residue by construction.
sub ident {
    my ($a, $b) = @_;
    my $n = length($a) < length($b) ? length($a) : length($b);
    return 0 unless $n;
    my $m = 0;
    for my $i (0 .. $n - 1) { $m++ if substr($a, $i, 1) eq substr($b, $i, 1) }
    my $denom = length($a) > length($b) ? length($a) : length($b);
    return $m / ($denom || 1);
}

# Candidate splice positions, as offsets into the locus.
# A donor is the first base of an intron and must read GT.
# An acceptor is the last base of an intron; plastid group II introns end in A
# followed by a pyrimidine, and the exon resumes at acceptor+1.
# GT is the group II donor and holds for most plastid introns, but not all: clpP
# intron 2 reads TG in four of five genomes and TT in the fifth, so requiring GT
# makes the true structure unreachable. --any-donor drops the requirement and
# leaves donor quality to the score, at the cost of a much larger search.
my (@donor, @acceptor);
for my $i (0 .. $L - 2) {
    push @donor,    $i if $ANY_DONOR || substr($locus, $i, 2) eq 'GT';
    push @acceptor, $i + 1 if substr($locus, $i, 1) =~ /^[AG]$/;
}
printf STDERR "locus %d nt, %d donor sites, %d acceptor sites, reference protein %d aa\n",
    $L, scalar(@donor), scalar(@acceptor), length($refprot) if $VERBOSE;

# Start codons in the first part of the locus. The gene is anchored here, so the
# true start sits near the beginning by construction.
my @starts;
for my $i (0 .. ($L > 400 ? 400 : $L - 3)) {
    push @starts, $i if $START_OK{ substr($locus, $i, 3) };
}

my @hits;

# Branch and bound. A partial structure is only extended while the protein it has
# produced so far still matches the reference prefix - without that the search
# over donor x acceptor pairs for every start position does not terminate in
# useful time.
my $MIN_PREFIX_ID = 0.45;

# A search for the wrong number of exons has no good structure to prune against
# and will grind indefinitely, so it is given a budget. Reaching it means "no
# confident answer" rather than "no answer exists", which is the right outcome
# when the exon count being tried is not the gene's.
my $MAX_NODES = $ENV{ANNOBTD_DENOVO_NODES} || 120000;
my $nodes = 0;

for my $s (@starts) {
    my @stack = ([ $s, 0, [] ]);                 # exon start, introns placed, segments
    while (@stack) {
        my ($estart, $placed, $segs) = @{ pop @stack };
        if (++$nodes > $MAX_NODES) { @stack = (); last }

        if ($placed == $NEXON - 1) {
            my $carried = 0;
            $carried += $_->[1] - $_->[0] + 1 for @$segs;
            my $phase = $carried % 3;
            my $need  = $phase ? 3 - $phase : 0;
            for (my $p = $estart + $need; $p + 2 < $L; $p += 3) {
                my $c = substr($locus, $p, 3);
                next unless exists $CODON{$c} && $CODON{$c} eq '*';
                my @full = (@$segs, [ $estart, $p + 2 ]);
                my $cds = join '', map { substr($locus, $_->[0], $_->[1] - $_->[0] + 1) } @full;
                next if length($cds) % 3;
                next if abs(length($cds) - $target_nt) > 0.15 * $target_nt;
                my $prot = tr_($cds);
                next if $prot =~ /\*.+/;
                $prot =~ s/\*$//;
                my $id = sim($prot, $refprot);
                next unless $id > 0.5;
                # Tiebreak on splice quality: a canonical A-ending acceptor is
                # preferred, but only between structures the protein cannot
                # separate, so it never overrides sequence evidence.
                my $sp = 0;
                for my $k (0 .. $#full - 1) {
                    my $istart = $full[$k][1] + 1;
                    my $iend   = $full[$k + 1][0] - 1;
                    $sp += 0.5  if substr($locus, $istart, 2) eq 'GT';
                    $sp += 0.25 if substr($locus, $istart, 2) =~ /^(GC|AT)$/;
                    $sp += 0.5  if substr($locus, $iend - 1, 2) =~ /^A[CT]$/;
                    $sp += 0.25 if substr($locus, $iend - 1, 2) =~ /^A[AG]$/;
                }
                push @hits, { id => $id, sp => $sp, segs => \@full, cds => length($cds) };
                last;
            }
            next;
        }

        for my $d (@donor) {
            next if $d < $estart + $MIN_EXON;
            last if $d > $L - $MIN_INTRON;

            # Score what this exon would contribute before exploring any acceptor.
            my @try = (@$segs, [ $estart, $d - 1 ]);
            my $part = join '', map { substr($locus, $_->[0], $_->[1] - $_->[0] + 1) } @try;
            my $whole = int(length($part) / 3) * 3;
            my $pp = $whole ? tr_(substr($part, 0, $whole)) : '';
            next if $pp =~ /\*/;                 # no stop inside an internal exon
            next if length($pp) > length($refprot);
            # Too short to judge: a two-residue prefix carries no signal, so let
            # the branch through rather than rejecting a genuinely tiny exon.
            if ($whole >= 15) {
            # Compare against a generous window of the reference prefix, with the
            # same indel-tolerant measure - a positional test here rejected the
            # correct Spinacia structure outright.
                my $win = length($pp) + 12;
                $win = length($refprot) if $win > length($refprot);
                next if sim(substr($refprot, 0, $win), $pp) < $MIN_PREFIX_ID;
            }

            my $carried = length($part) % 3;
            my $need    = $carried ? 3 - $carried : 0;
            for my $a (@acceptor) {
                next if $a < $d + $MIN_INTRON;
                last if $a > $d + $MAX_INTRON;
                next if $a + $MIN_EXON >= $L;
                # The exon opening at $a+1 must extend the protein, not merely
                # exist. Check its first few codons before taking the branch.
                my $probe = substr($locus, $a + 1 + $need, 24);
                next if length($probe) < 12;
                my $pr = tr_(substr($probe, 0, int(length($probe) / 3) * 3));
                next if $pr =~ /\*/;
                next unless sim($refprot, $pr) > 0.25;
                push @stack, [ $a + 1, $placed + 1, [ @try ] ];
            }
        }
    }
}

printf STDERR "search visited %d nodes%s\n", $nodes,
    ($nodes > $MAX_NODES ? " (budget reached)" : "") if $VERBOSE;
unless (@hits) {
    print "no structure found satisfying the constraints\n";
    exit 1;
}
@hits = sort { $b->{id} <=> $a->{id} || $b->{sp} <=> $a->{sp} } @hits;

# Report the best few, converted back to genome coordinates.
my $shown = 0;
for my $h (@hits) {
    last if $shown++ >= 3;
    my @coords;
    for my $sg (@{ $h->{segs} }) {
        my ($gs, $ge);
        if ($strand eq '-') { $gs = $lend - $sg->[1]; $ge = $lend - $sg->[0] }
        else                { $gs = $lstart + $sg->[0]; $ge = $lstart + $sg->[1] }
        push @coords, "$gs-$ge";
    }
    @coords = reverse @coords if $strand eq '-';
    printf "%.3f  splice %.2f  %d nt  %s\n", $h->{id}, $h->{sp}, $h->{cds}, join("  ", @coords);
}
printf STDERR "%d candidate structures scored\n", scalar(@hits) if $VERBOSE;
