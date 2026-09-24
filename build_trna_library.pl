#!/usr/bin/perl
# Build the curated plastid tRNA sequence library from truth tables.
#
# post_filter_annotation.pl names each predicted tRNA by the identity of its best
# BLAST hit in this library (--trna-library). Names transferred from guides are
# only as good as GenBank, which confuses the three CAU tRNAs (trnI-CAU, trnM-CAU,
# trnfM-CAU), swaps trnG-GCC and trnG-UCC, and hands a spliced tRNA's exons to the
# wrong letter when no close guide exists; the sequence settles all of these.
#
# USAGE: build_trna_library.pl <out.fasta> <truth.tsv> [truth.tsv ...]
#   each truth table needs its genome FASTA beside it (<name>.fsa) and
#   plastid_trna_names.txt beside this script.
#
# One entry per identity per genome (IR copies collapse), so a genome can be left
# out by name for an honest leave-one-out (post filter --trna-library-exclude).
# A short row (<= 60 nt) of a spliceable letter is an exon even without _exonN.
# A bare truth name (Arabidopsis "trnS") takes its anticodon from the overlapping
# ARAGORN gene; the six intron-containing identities are fixed by letter, since
# ARAGORN misnames those.
#   >identity|accession|genome
use strict; use warnings; use FindBin;
my ($out, @truths) = @ARGV; die "usage: $0 <out.fasta> <truth.tsv>...\n" unless @truths;
my $aragorn = "$FindBin::Bin/bin/aragorn";
my %ok; if (open my $nh, "<", "$FindBin::Bin/plastid_trna_names.txt") { while (<$nh>) { next if /^#/; $ok{$_} = 1 for /(trn\S+)/g } close $nh }
my %spliced = (trnI => 'trnI-GAU', trnA => 'trnA-UGC', trnK => 'trnK-UUU', trnG => 'trnG-UCC', trnL => 'trnL-UAA', trnV => 'trnV-UAC');
open my $o, ">", $out or die "$out: $!"; my ($nent, %per_id);
for my $truth (@truths) {
    (my $fsa = $truth) =~ s/\.truth\.tsv$/.fsa/; (my $genome = $truth) =~ s{.*/}{}; $genome =~ s/\.truth\.tsv$//;
    open my $fh, "<", $fsa or do { warn "no FASTA for $truth\n"; next };
    my ($acc, $seq) = ('', ''); while (<$fh>) { if (/^>(\S+)/) { $acc = $1; next } s/\s//g; $seq .= uc $_ } close $fh;
    my @ar; if (-x $aragorn && open my $ah, "-|", "$aragorn -t -i -w -gcbact $fsa 2>/dev/null") {
        while (<$ah>) { next unless /^\d+\s+tRNA-\w+\s+(c?)\[(\d+),(\d+)\]\s+\d+\s+\((\w+)\)/; (my $ac = uc $4) =~ tr/T/U/; push @ar, [$2, $3, $1 ? '-' : '+', $ac] } close $ah }
    # group exon rows of one base name and strand within 3.5 kb
    my (@groups, %seen_id); open my $th, "<", $truth or die;
    while (<$th>) { chomp; my @f = split /\t/; next unless ($f[4] // '') eq 'tRNA' && $f[0] =~ /^trn/;
        # exon rows carry _exonN, or (Daucus, Sorghum) the plain gene name on each short piece
        (my $base = $f[0]) =~ s/_exon\d+$//; my $ex = ($f[0] =~ /_exon\d+$/ || ($spliced{ $base =~ /^(trn\w)/ ? $1 : '' } && $f[2] - $f[1] + 1 <= 60)) ? 1 : 0; my $g;
        if ($ex) { for my $c (@groups) { next unless $c->{ex} && $c->{base} eq $base && $c->{d} eq $f[3] && @{ $c->{rows} } < 2; if (grep { abs($_->[0] - $f[1]) < 3500 } @{ $c->{rows} }) { $g = $c; last } } }
        unless ($g) { $g = { base => $base, d => $f[3], ex => $ex, rows => [] }; push @groups, $g } push @{ $g->{rows} }, [$f[1], $f[2]] }
    close $th;
    for my $g (@groups) { next if $g->{ex} && @{ $g->{rows} } != 2;
        my @r = sort { $a->[0] <=> $b->[0] } @{ $g->{rows} }; my ($s0, $e0) = ($r[0][0], $r[-1][1]);
        my $t = join '', map { substr($seq, $_->[0] - 1, $_->[1] - $_->[0] + 1) } @r; if ($g->{d} eq '-') { $t = reverse $t; $t =~ tr/ACGT/TGCA/ }
        my $id; if ($g->{base} =~ /^trn[A-Za-z]{1,2}-[ACGU]{3}$/) { $id = $g->{base} }
        elsif ($g->{base} =~ /^(trn[A-Za-z]{1,2})$/) { my $letter = $1;
            if ($g->{ex}) { $id = $spliced{$letter} }
            else { my ($bo, $bx) = (0, undef); for my $x (@ar) { next unless $x->[2] eq $g->{d}; my $lo = $s0 > $x->[0] ? $s0 : $x->[0]; my $hi = $e0 < $x->[1] ? $e0 : $x->[1]; my $ov = $hi - $lo + 1; ($bo, $bx) = ($ov, $x) if $ov > $bo }
                $id = "$letter-$bx->[3]" if $bx && $bo >= 0.5 * ($e0 - $s0 + 1) } }
        unless ($id && $ok{$id}) { warn sprintf("skip %s %s %d-%d %s: identity %s\n", $genome, $g->{base}, $s0, $e0, $g->{d}, $id // 'unresolved'); next }
        next if $seen_id{$id}++;   # one entry per identity per genome (IR copies)
        print $o ">$id|$acc|$genome\n$t\n"; $nent++; $per_id{$id}++ }
}
close $o; printf STDERR "%d entries, %d identities -> %s\n", $nent, scalar keys %per_id, $out;
printf STDERR "  %-10s %d\n", $_, $per_id{$_} for sort keys %per_id;
