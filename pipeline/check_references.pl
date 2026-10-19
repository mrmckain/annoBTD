#!/usr/bin/perl
use strict;
use warnings;

# Screen a guide reference set for bad annotations before they are used.
#
# USAGE: check_references.pl <genes.fsa> [--tsv out.tsv] [--drop out.fsa]
#
# <genes.fsa> is the pipeline's *_annotated_regions_fromverdant_genes.fsa: one
# record per gene per guide, named <gene>[_exonN]XXX<guide>.
#
# GenBank holds real annotation errors, and annoBTD transfers whatever the winning
# reference says - so a reference whose own coding sequence does not translate is
# a direct route to a bad annotation. Two independent signals are used:
#
#   intrinsic  - the reference's own protein is malformed: internal stop codons,
#                no recognised start, length not a multiple of three. This needs
#                no comparison and catches an outright bad record. ON by default.
#   comparative- the reference disagrees with the other guides for the same gene,
#                in length or in sequence identity. OFF by default: --comparative.
#
# The comparative checks are off because they were measured and they cost more
# than they returned. Their premise - that a guide far from the peer median is
# probably mis-annotated - does not hold for plastid genes, whose length varies
# genuinely between lineages. Over a 9-genome leave-one-out they dropped 39-47
# reference copies per genome and the roll-call is the genes best known for real
# length variation: rps18 (492 vs median 306), cemA (1305 vs 693), accD (951 vs
# 1476), rpl22 (600 vs 468), infA (324 vs 237 - the grass copies).
#
# The cost was 10 features lost across the set (277 -> 287 missing) for +0.2 points
# of boundary exactness, which is noise. infA was the clearest casualty: screening
# removed the only close references, leaving the scorer to compare a grass ORF
# against eudicot guides at 72% identity, and infA then went uncalled in all seven
# genomes that have it.
#
# So a length outlier is reported when asked for, but it is not evidence of a bad
# annotation on its own. If the check is ever re-enabled by default it should
# compare each guide against its k-mer-nearest peers rather than a global median,
# so a grass is judged against grasses.
#
# --drop writes a copy of the FASTA with the flagged references removed, ready to
# feed back into the pipeline.

my ($genes_file, @rest) = @ARGV;
die "usage: $0 <genes.fsa> [--tsv out.tsv] [--drop out.fsa]\n" unless defined $genes_file;

my ($tsv_file, $drop_file, $comparative);
while (@rest) {
    my $a = shift @rest;
    $tsv_file  = shift @rest if $a eq '--tsv';
    $drop_file = shift @rest if $a eq '--drop';
    $comparative = 1         if $a eq '--comparative';
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

# Read the reference FASTA and reassemble each guide's copy of each gene by
# joining its exons in transcription order (exon1 first).
my ($id, %seq, @ids);
open my $fh, "<", $genes_file or die "Cannot open $genes_file: $!";
while (<$fh>) { chomp; if (/^>(.+)/) { $id = $1; push @ids, $id } else { $seq{$id} .= $_ } }
close $fh;
die "No records in $genes_file\n" unless @ids;

my %gene;      # gene -> guide -> { seq, parts }
for my $k (@ids) {
    next unless $k =~ /^(.+?)XXX(.+)$/;
    my ($name, $guide) = ($1, $2);
    my $exon = 0;
    if ($name =~ /^(.+)_exon(\d+)$/) { ($name, $exon) = ($1, $2) }
    next if $name =~ /^(trn|rrn)/;          # not translated
    next if $name =~ /^rps12(_|$)/;         # trans-spliced; exons are not contiguous
    push @{ $gene{$name}{$guide} }, { n => $exon, s => uc($seq{$k} // ''), id => $k };
}

sub translate {
    my $s = shift;
    my $p = '';
    for (my $i = 0; $i + 2 < length $s; $i += 3) {
        my $c = substr($s, $i, 3);
        $p .= exists $CODON{$c} ? $CODON{$c} : 'X';
    }
    return $p;
}
sub median { my @v = sort { $a <=> $b } @_; return @v ? $v[@v/2] : 0 }

# Fraction of positions shared, compared as a bag of 5-mers. Alignment free, so an
# indel does not throw the whole comparison off the way a positional match would.
sub kmer_sim {
    my ($a, $b) = @_;
    return 0 unless length($a) >= 5 && length($b) >= 5;
    my %A;
    $A{ substr($a, $_, 5) } = 1 for 0 .. length($a) - 5;
    my ($hit, $tot) = (0, 0);
    for my $i (0 .. length($b) - 5) { $tot++; $hit++ if $A{ substr($b, $i, 5) } }
    return $tot ? $hit / $tot : 0;
}

my (@flag, %drop);
for my $g (sort keys %gene) {
    my (%prot, %len);
    for my $guide (keys %{ $gene{$g} }) {
        my @ex = sort { $a->{n} <=> $b->{n} } @{ $gene{$g}{$guide} };
        my $cds = join '', map { $_->{s} } @ex;
        next unless length $cds;
        $prot{$guide} = translate($cds);
        $len{$guide}  = length $cds;
    }
    next unless keys %len;
    my $med = median(values %len);

    # Mean similarity of each copy to the others. Comparing against a fixed
    # threshold made one bad copy drag its innocent peers below it too, so each
    # copy is judged against the typical similarity in this gene instead.
    my %sim;
    if (keys(%prot) >= 3) {
        for my $a (keys %prot) {
            my @s = map { kmer_sim($prot{$_}, $prot{$a}) } grep { $_ ne $a } keys %prot;
            my $m = 0; $m += $_ for @s; $sim{$a} = @s ? $m / @s : 1;
        }
    }
    my $typical = %sim ? median(values %sim) : 1;

    for my $guide (sort keys %prot) {
        my $p = $prot{$guide};
        my @why;

        # Intrinsic checks on the reference's own coding sequence.
        push @why, "length $len{$guide} not a multiple of 3" if $len{$guide} % 3;
        my $body = $p; $body =~ s/\*$//;
        my $internal = ($body =~ tr/\*//);
        push @why, "$internal internal stop" . ($internal > 1 ? "s" : "") if $internal;

        my @ex = sort { $a->{n} <=> $b->{n} } @{ $gene{$g}{$guide} };
        my $start = substr(join('', map { $_->{s} } @ex), 0, 3);
        # Reported, but NOT fatal - see the note on %START_OK. A start can be
        # declared by /transl_except rather than inferable from sequence, so an
        # unrecognised one is a curiosity, not proof the reference is unusable.
        # Dropping on it removed Nicotiana psbI and Spinacia ndhD outright: in a
        # self-annotation there is only one reference, so dropping it loses the
        # gene. Only an internal stop or a broken frame makes a reference unusable.
        my $odd_start = $START_OK{$start} ? 0 : 1;
        push @why, "start codon $start (reported only)" if $odd_start;

        # Comparative checks, only meaningful with peers to compare against - and
        # only when explicitly asked for. See the note at the top of the file.
        if ($comparative && keys(%len) >= 3 && $med) {
            my $ratio = $len{$guide} / $med;
            push @why, sprintf("length %d vs peer median %d (%.2fx)", $len{$guide}, $med, $ratio)
                if $ratio < 0.85 || $ratio > 1.18;

            if (exists $sim{$guide} && $typical > 0.3
                && $sim{$guide} < 0.5 * $typical && $sim{$guide} < 0.5) {
                push @why, sprintf("protein shares %.0f%% of 5-mers with its peers, against %.0f%% typical for this gene",
                    100 * $sim{$guide}, 100 * $typical);
            }
        }

        next unless @why;
        push @flag, [ $g, $guide, $len{$guide}, join("; ", @why) ];
        # Drop only when something other than the start codon is wrong.
        next if $odd_start && @why == 1;
        $drop{$_->{id}} = 1 for @{ $gene{$g}{$guide} };
    }
}

printf "=== %s: %d genes, %d reference copies flagged (%s) ===\n",
    $genes_file, scalar(keys %gene), scalar(@flag),
    $comparative ? "intrinsic + comparative" : "intrinsic only";
for my $f (sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @flag) {
    printf "  %-16s guide %-6s %5d nt  %s\n", @$f;
}

if ($tsv_file) {
    open my $th, ">", $tsv_file or die "Cannot write $tsv_file: $!";
    print $th join("\t", qw(gene guide cds_len problems)), "\n";
    print $th join("\t", @$_), "\n" for @flag;
    close $th;
}
if ($drop_file) {
    open my $dh, ">", $drop_file or die "Cannot write $drop_file: $!";
    my $kept = 0;
    for my $k (@ids) {
        next if $drop{$k};
        print $dh ">$k\n$seq{$k}\n";
        $kept++;
    }
    close $dh;
    printf STDERR "wrote %s: %d of %d records kept\n", $drop_file, $kept, scalar(@ids);
}
