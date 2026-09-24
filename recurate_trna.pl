#!/usr/bin/perl
# Recurate tRNA boundaries in a truth table against ARAGORN.
#
# Truth tRNA windows inherited from GenBank carry submission-era conventions; the
# grass tables place every inverted-repeat tRNA 4 nt from where the gene folds.
# ARAGORN (Laslett & Canback 2004) finds the gene; but its reported window is one
# base longer at each end than the curated convention for many tRNAs, so a per-
# identity convention table (truth minus ARAGORN offsets, derived from tables
# that were curated by hand) turns its call into the house convention. A window
# is replaced only when it disagrees with that expectation.
#
# USAGE: recurate_trna.pl <truth.tsv> <aragorn_batch_output.txt> <convention.tsv> <out.tsv> [--log FILE]
#   convention.tsv: identity (trnX-YYY) <TAB> off5/off3   (truth minus ARAGORN, transcript orientation)
use strict; use warnings;
my ($truth, $arout, $conv, $out, @rest) = @ARGV; die "usage: $0 <truth.tsv> <aragorn.txt> <convention.tsv> <out.tsv> [--log FILE]\n" unless $out;
my %o; while (@rest) { my $a = shift @rest; $o{$1} = shift @rest if $a =~ /^--(.+)$/ }
my %aa3 = (Ala=>'A',Arg=>'R',Asn=>'N',Asp=>'D',Cys=>'C',Gln=>'Q',Glu=>'E',Gly=>'G',His=>'H',Ile=>'I',Leu=>'L',Lys=>'K',Met=>'M',Phe=>'F',Pro=>'P',Ser=>'S',Thr=>'T',Trp=>'W',Tyr=>'Y',Val=>'V',fMet=>'fM');
# curated exon lengths per spliced identity, beside the convention file
my %exlen; { (my $ef = $conv) =~ s{[^/]*$}{trna_exon_lengths.tsv}; if (open my $eh, "<", $ef) { while (<$eh>) { next if /^#/; my ($id, $l) = split; $exlen{$id} = [$1, $2] if $l && $l =~ m{^(\d+)/(\d+)$} } close $eh } }
my %conv; open my $c, "<", $conv or die; while (<$c>) { my ($id, $off) = split; next unless $off && $off =~ m{^(-?\d+)/(-?\d+)$}; $conv{$id} = [$1, $2] } close $c;
my @ar; open my $afh, "<", $arout or die;
while (<$afh>) { next unless /^\d+\s+tRNA-(\w+)\s+(c?)\[(\d+),(\d+)\]\s+\d+\s+\((\w+)\)(.*)/; my ($aa, $cc, $s, $e, $ac, $rest) = ($1, $2, $3, $4, $5, $6);
    my ($ipos, $ilen) = $rest =~ /i\((\d+),(\d+)\)/ ? ($1, $2) : (0, 0); (my $acu = uc $ac) =~ tr/T/U/;
    push @ar, { id => "trn" . ($aa3{$aa} // $aa) . "-$acu", s => $s, e => $e, d => ($cc ? '-' : '+'), ipos => $ipos, ilen => $ilen } }
close $afh;
(my $fsa = $truth) =~ s/\.truth\.tsv.*$/.fsa/; my $genome = ''; if (open my $fh, "<", $fsa) { <$fh>; $genome = do { local $/; <$fh> }; $genome =~ s/\s//g; $genome = uc $genome; close $fh } else { die "no genome fasta beside $truth ($fsa)\n" }
my @rows; open my $h, "<", $truth or die; while (<$h>) { chomp; push @rows, [split /\t/] } close $h;
# group truth tRNA rows into genes
my @genes; for my $i (0 .. $#rows) { my $f = $rows[$i]; next unless ($f->[4] // '') eq 'tRNA'; my $ex = ($f->[0] =~ /_exon(\d+)$/) ? $1 : 0; (my $n = $f->[0]) =~ s/_exon\d+$//; my $grp;
    my $short = ($f->[2] - $f->[1] + 1) <= 60;   # an exon of a spliced tRNA, not a whole one
    for my $g (@genes) { next unless $g->{n} eq $n && $g->{d} eq $f->[3]; next unless $ex ? $g->{ex} : ($short && $g->{short});
        if (grep { abs($rows[$_][1] - $f->[1]) < 3500 } @{ $g->{rows} }) { $grp = $g; last } }
    unless ($grp) { $grp = { n => $n, d => $f->[3], rows => [], ex => $ex ? 1 : 0, short => $short }; push @genes, $grp } push @{ $grp->{rows} }, $i; $grp->{ex} = 1 if @{ $grp->{rows} } == 2 && !$ex }
my $log = $o{log} ? do { open my $l, ">>", $o{log} or die; $l } : \*STDERR; my ($nchg, $nkeep, $nnone) = (0, 0, 0);
(my $gname = $truth) =~ s{.*/}{}; $gname =~ s/\.truth\.tsv$//;
my %spliced_id = (trnI => "trnI-GAU", trnA => "trnA-UGC", trnK => "trnK-UUU", trnG => "trnG-UCC", trnL => "trnL-UAA", trnV => "trnV-UAC");
for my $g (@genes) { my @r = map { $rows[$_] } @{ $g->{rows} }; my ($s, $e) = (1e12, 0); for my $f (@r) { $s = $f->[1] if $f->[1] < $s; $e = $f->[2] if $f->[2] > $e }
    my ($best, $bo) = (undef, 0); for my $x (@ar) { next unless $x->{d} eq $g->{d}; next unless ($x->{ipos} ? 1 : 0) == ($g->{ex} ? 1 : 0); my $lo = $s > $x->{s} ? $s : $x->{s}; my $hi = $e < $x->{e} ? $e : $x->{e}; my $ov = $hi - $lo + 1; if ($ov > $bo) { ($best, $bo) = ($x, $ov) } }
    if (!$best || $bo < 0.5 * ($e - $s + 1)) {
        if ($g->{ex} && @r == 2) { $best = { s => $s, e => $e, id => ($spliced_id{ $g->{n} } // $g->{n}), ipos => 1, ilen => 0, noaragorn => 1 } }
        else { $nnone++; printf $log "%s\t%s\t%d-%d\t%s\tNO_ARAGORN_GENE\n", $gname, $g->{n}, $s, $e, $g->{d}; next } }
    my $id = $g->{n} =~ /^trn\w+-[ACGU]{3}$/ ? $g->{n} : $g->{ex} ? ($spliced_id{ $g->{n} } // $best->{id}) : $best->{id};
    # expected window: ARAGORN gene ends placed structurally, plus the identity's residual
    my ($ns, $ne) = $best->{noaragorn} ? ($s, $e) : trna_place(\$genome, $best->{s}, $best->{e}, $g->{d}, $id, \%conv);
    # a spliced gene whose ends already agree still gets its exon lengths normalised
    my $keep_ends = ($ns == $s && $ne == $e);
    if ($keep_ends && !($g->{ex} && @r == 2 && $best->{ipos})) { $nkeep++; next }
    # rewrite: exon rows from the intron (gene-relative, transcript orientation), else the single row
    my @new;
    if ($best->{ipos} && $g->{ex} && @r != 2) { printf $log "%s\t%s\t%d-%d\t%s\tSKIPPED_ODD_EXONS\n", $gname, $g->{n}, $s, $e, $g->{d}; next }
    if ($best->{ipos} && $g->{ex}) {
        # ARAGORN's intron split is off by a base or two; keep the curated exon lengths and
        # anchor them to the corrected gene ends (exon1 at the 5' end, exon2 at the 3' end)
        my ($r1, $r2) = $g->{d} eq '+' ? (sort { $a->[1] <=> $b->[1] } @r) : (sort { $b->[1] <=> $a->[1] } @r);   # exon1, exon2 in transcript order
        my $eid = $g->{n} =~ /-[ACGU]{3}$/ ? $g->{n} : $id; my ($l1, $l2) = $exlen{$eid} ? @{ $exlen{$eid} } : ($r1->[2] - $r1->[1] + 1, $r2->[2] - $r2->[1] + 1);
        if ($keep_ends && $l1 == $r1->[2] - $r1->[1] + 1 && $l2 == $r2->[2] - $r2->[1] + 1) { $nkeep++; next }
        if ($g->{d} eq '+') { @new = ([$ns, $ns + $l1 - 1], [$ne - $l2 + 1, $ne]) } else { @new = ([$ne - $l1 + 1, $ne], [$ns, $ns + $l2 - 1]) } }
    else { @new = ([$ns, $ne]) }
    if ($g->{ex} && @new == 2 && @r == 2) { my ($o1, $o2) = $g->{d} eq '+' ? (sort { $a->[1] <=> $b->[1] } @r) : (sort { $b->[1] <=> $a->[1] } @r);
        if (1) {
            printf $log "%s\t%s\t%d-%d,%d-%d\t%s\tREPLACED\t%d-%d,%d-%d\n", $gname, $g->{n}, $o1->[1], $o1->[2], $o2->[1], $o2->[2], $g->{d}, $new[0][0], $new[0][1], $new[1][0], $new[1][1];
            ($o1->[1], $o1->[2]) = @{ $new[0] }; ($o2->[1], $o2->[2]) = @{ $new[1] }; $nchg++ }
        else { printf $log "%s\t%s\t%d-%d\t%s\tSKIPPED_ODD_EXONS\n", $gname, $g->{n}, $s, $e, $g->{d} } }
    elsif (!$g->{ex} && @new == 1 && @r == 1) { printf $log "%s\t%s\t%d-%d\t%s\tREPLACED\t%d-%d\n", $gname, $g->{n}, $s, $e, $g->{d}, $ns, $ne; ($r[0][1], $r[0][2]) = ($ns, $ne); $nchg++ }
    else { printf $log "%s\t%s\t%d-%d\t%s\tSKIPPED_SHAPE_MISMATCH\n", $gname, $g->{n}, $s, $e, $g->{d} } }
open my $w, ">", $out or die; print $w join("\t", @$_), "\n" for @rows; close $w;
printf $log "%s\tSUMMARY\treplaced %d, unchanged %d, no ARAGORN gene %d\n", $gname, $nchg, $nkeep, $nnone;

# Structural placement of a tRNA gene inside an ARAGORN window. ARAGORN pads its
# window with a flanking base at either end whenever that base can pair, so the
# window is trimmed (0-2 bases per end) to the frame in which positions 1-7 pair
# with the seven bases before the last (discriminator) base - most pairs, then U8,
# then the smallest trim; His pairs its first base with its last (G-1:C73). A few
# identities then carry a fixed residual (trna_window_residuals.tsv: trnT, trnH,
# trnM, trnP), consistent across every hand-curated source.
sub trna_trim { my ($t, $his) = @_; my %pair = map { $_ => 1 } qw(AT TA GC CG GT TG); my ($bs, $bt, $bsc) = (0, 0, -1);
    for my $a (0 .. 2) { for my $b (0 .. 2) { my $w = substr($t, $a, length($t) - $a - $b); my $L = length $w; next if $L < 60;
        my $n = 0; for my $i (0 .. 6) { $n++ if $pair{ substr($w, $i, 1) . substr($w, $L - ($his ? 1 : 2) - $i, 1) } }
        my $sc = $n * 10 + (substr($w, 7, 1) eq 'T' ? 3 : 0) + (2 - ($a > $b ? $a : $b));
        if ($sc > $bsc) { ($bs, $bt, $bsc) = ($a, $b, $sc) } } } return ($bs, $bt) }
# gene ends from an ARAGORN window [as,ae] on strand d, sequence ref \$seq, identity id
sub trna_place { my ($seq, $as, $ae, $d, $id, $resid) = @_; my $t = substr($$seq, $as - 1, $ae - $as + 1); if ($d eq '-') { $t = reverse $t; $t =~ tr/ACGT/TGCA/ }
    my ($a, $b) = trna_trim($t, $id =~ /^trnH/ ? 1 : 0); my $off = $resid->{$id} // [0, 0];
    return $d eq '+' ? ($as + $a + $off->[0], $ae - $b + $off->[1]) : ($as + $b - $off->[1], $ae - $a - $off->[0]) }
1;
