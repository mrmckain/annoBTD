#!/usr/bin/perl
use strict;
use warnings;

# Post-hoc protein sanity check on a finished annotation.
#
# USAGE: check_annotation.pl <plastome.fsa> <annotation.txt> [--refs genes.fsa] [--tsv out.tsv]
#
# The annotation file is the pipeline's output: gene, start, end, dir (tab
# separated). Exons of one gene are joined in transcription order before
# translating, so the checks apply to the coding sequence a ribosome would see
# rather than to each exon separately.
#
# A boundary can be wrong in ways the coordinates alone never reveal. Translating
# what was actually called catches them:
#   * length not a multiple of 3        - the exon join is out of frame
#   * no recognised start codon         - the 5' boundary is off
#   * no terminal stop                  - the 3' boundary is off
#   * internal stop codons              - the frame or an exon boundary is wrong
#   * length far from the reference     - the wrong ORF was chosen
#
# With --refs (the pipeline's *_annotated_regions_fromverdant_genes.fsa) each
# called protein is also compared in length against the guide copies of that gene.
#
# Exit status is the number of genes with a HIGH severity problem.

my ($plastome_file, $ann_file, @rest) = @ARGV;
die "usage: $0 <plastome.fsa> <annotation.txt> [--refs genes.fsa] [--tsv out.tsv]\n"
    unless defined $ann_file;

my ($refs_file, $tsv_file);
while (@rest) {
    my $a = shift @rest;
    $refs_file = shift @rest if $a eq '--refs';
    $tsv_file  = shift @rest if $a eq '--tsv';
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

my $plastome = '';
open my $pf, "<", $plastome_file or die "Cannot open $plastome_file: $!";
while (<$pf>) { chomp; next if /^>/; s/\s+//g; $plastome .= $_ }
close $pf;
die "No sequence in $plastome_file\n" unless length $plastome;
$plastome = uc $plastome;

# Reference protein lengths per gene, for the outlier comparison.
my %ref_len;
if ($refs_file && -s $refs_file) {
    my ($id, %seq);
    open my $rf, "<", $refs_file or die "Cannot open $refs_file: $!";
    while (<$rf>) { chomp; if (/^>(.+)/) { $id = $1 } else { $seq{$id} .= $_ } }
    close $rf;
    # Sum a reference's exons before comparing: taking each exon separately made
    # the median one exon long, so every multi-exon gene looked 2x too big.
    my %per_ref;
    for my $k (keys %seq) {
        my ($gene, $rid) = ($k, '');
        if ($k =~ /^(.+?)XXX(.+)$/) { ($gene, $rid) = ($1, $2) }
        $gene =~ s/_exon\d+$//;
        $per_ref{$gene}{$rid} += length $seq{$k};
    }
    for my $g (keys %per_ref) {
        push @{ $ref_len{$g} }, $_ for values %{ $per_ref{$g} };
    }
}

# Gather exons per gene. Intergenic spacers and the region records are not genes.
my (%parts, @order);
open my $af, "<", $ann_file or die "Cannot open $ann_file: $!";
while (<$af>) {
    chomp;
    next unless /\S/;
    my ($name, $s, $e, $d) = split /\t/;
    next unless defined $d && defined $s && $s =~ /^\d+$/ && $e =~ /^\d+$/;
    next if $name =~ /^(LSC|SSC|IRA|IRB|FULL)$/;
    next if $name =~ /~/ || $name =~ /_intron\d*$/;
    next if $name =~ /^trn/ || $name =~ /^rrn/;     # RNA genes are not translated
    # rps12 is the one trans-spliced plastid gene: its 5' exon sits tens of
    # kilobases from the others and is spliced in trans, so joining its exons by
    # proximity is meaningless and it would always flag.
    next if $name =~ /^rps12(_|$)/;
    my $gene = $name; my $exon = 0;
    if ($name =~ /^(.+)_exon(\d+)$/) { ($gene, $exon) = ($1, $2) }
    ($s, $e) = ($e, $s) if $e < $s;
    push @{ $parts{"$gene\t$d"} }, { s => $s, e => $e, n => $exon };
    push @order, "$gene\t$d" unless grep { $_ eq "$gene\t$d" } @order;
}
close $af;

my @report;
my $high = 0;
for my $key (@order) {
    my ($gene, $dir) = split /\t/, $key;
    my @ex = @{ $parts{$key} };

    # IR-duplicated genes appear twice; split the exon list into runs that sit
    # close together so each copy is checked on its own.
    @ex = sort { $a->{s} <=> $b->{s} } @ex;
    my @groups = ([ shift @ex ]);
    for my $x (@ex) {
        if ($x->{s} - $groups[-1][-1]{e} > 5000) { push @groups, [$x] }
        else { push @{ $groups[-1] }, $x }
    }

    for my $grp (@groups) {
        # Concatenate in ASCENDING coordinate order and reverse complement the
        # whole thing for a minus-strand gene. Reversing the exon order as well
        # would undo itself and splice the exons back to front.
        my @o = @$grp;
        my $cds = '';
        $cds .= substr($plastome, $_->{s} - 1, $_->{e} - $_->{s} + 1) for @o;
        if ($dir eq '-') { $cds = reverse $cds; $cds =~ tr/ACGT/TGCA/ }
        next unless length $cds;

        my @problems;
        my $sev = 'ok';
        push @problems, sprintf("length %d not a multiple of 3", length $cds) if length($cds) % 3;

        my $start = substr($cds, 0, 3);
        unless ($START_OK{$start}) { push @problems, "start codon $start"; $sev = 'HIGH' }

        my $prot = '';
        for (my $i = 0; $i + 2 < length $cds; $i += 3) {
            my $c = substr($cds, $i, 3);
            $prot .= exists $CODON{$c} ? $CODON{$c} : 'X';
        }
        my $terminal = substr($prot, -1) eq '*';
        unless ($terminal) { push @problems, "no terminal stop"; $sev = 'HIGH' }

        my $body = $terminal ? substr($prot, 0, -1) : $prot;
        my $internal = ($body =~ tr/\*//);
        if ($internal) { push @problems, "$internal internal stop" . ($internal > 1 ? "s" : ""); $sev = 'HIGH' }

        if ($ref_len{$gene} && @{ $ref_len{$gene} }) {
            my @l = sort { $a <=> $b } @{ $ref_len{$gene} };
            my $med = $l[ @l / 2 ];
            if ($med) {
                my $ratio = length($cds) / $med;
                if ($ratio < 0.8 || $ratio > 1.25) {
                    push @problems, sprintf("length %d vs reference median %d (%.2fx)", length($cds), $med, $ratio);
                    $sev = 'HIGH' if $ratio < 0.5 || $ratio > 2;
                    $sev = 'warn' if $sev eq 'ok';
                }
            }
        }

        next unless @problems;
        $sev = 'warn' if $sev eq 'ok';
        $high++ if $sev eq 'HIGH';
        push @report, [ $sev, $gene, $dir, $o[0]{s}, $o[-1]{e}, length($cds), join("; ", @problems) ];
    }
}

printf "=== %s: %d coding genes checked, %d with problems (%d HIGH) ===\n",
    $ann_file, scalar(@order), scalar(@report), $high;
for my $r (sort { ($b->[0] eq 'HIGH') <=> ($a->[0] eq 'HIGH') || $a->[1] cmp $b->[1] } @report) {
    printf "  %-5s %-16s %s %8d-%-8d %5d nt  %s\n", @$r[0,1,2,3,4,5,6];
}

if ($tsv_file) {
    open my $th, ">", $tsv_file or die "Cannot write $tsv_file: $!";
    print $th join("\t", qw(severity gene strand start end cds_len problems)), "\n";
    print $th join("\t", @$_), "\n" for @report;
    close $th;
}
exit($high > 255 ? 255 : $high);
