package ConservedStart;
# The start of a protein gene from the N-terminus that distant lineages conserve.
#
# GenBank start codons follow each lab's habit, and a family can share a wrong one: about
# 740 of 950 Liliaceae records start ndhK 37 codons upstream of the MVL... N-terminus that
# tobacco, the grasses and Marchantia conserve (the extension is the first ATG, read through
# ndhC in another frame), and the Tulipa/Amana records start rbcL 7 codons before the
# MSPQTE... that every land plant conserves. Neither the guides' majority nor the first
# in-frame ATG finds those starts; conservation across orders does.
#
# tables/nterm_refs.tsv (refdb_build/build_nterm_refs.pl) holds the first 30 aa of every
# protein gene of the hand-curated truth genomes, which span seven orders from Marchantia
# to the grasses. Each reference peptide is aligned (Smith-Waterman) to the target's
# translation in the call's frame; where it aligns (score >= 20) it says where its first
# residue falls (its first three residues must be in the alignment). A position counts when references from at least two ORDERS put their first
# residue within 2 codons of it, so one well-sampled order (30 curated grasses) is one vote.
# The references are evidence in two directions only; small differences where lineages genuinely annotate a start
# differently (ndhK at MVL... in tobacco and the curated grasses, at MN... in Arabidopsis and the RefSeq grasses; a
# 7-codon rpoC1 extension in Apiales and Solanales) are left to the guides:
#   - a guide majority may not extend a call whose start two orders already share and none puts upstream
#     (veto_extension: Liliaceae ndhK stays at MVL... instead of the 894 nt extension through ndhC)
#   - a call is trimmed when no order supports its start and two orders agree downstream (trim: rbcL from the
#     Tulipa records' upstream ATG back to MSPQTE...)
use strict; use warnings;

my %CODON; { my @b = qw(T C A G); my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; my $i = 0;
    for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1) } } } }
sub translate { my $s = shift; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr($s, $i, 3)} // 'X' } $p }

# load(file, exclude list) -> { gene => [[order, peptide, accession], ...] }: at most 3 distinct peptides per
# order per gene; references whose accession (with or without version) or organism is excluded are skipped
sub load { my ($file, $excl) = @_; my %ex; for (split /,/, ($excl // '')) { next unless length; $ex{$_} = 1; (my $b = $_) =~ s/\.\d+$//; $ex{$b} = 1 }
    my (%refs, %n); open my $h, "<", $file or return {};
    while (<$h>) { chomp; my ($g, $acc, $ord, $org, $p) = split /\t/; next if $g eq 'gene' || !$p; (my $b = $acc) =~ s/\.\d+$//; (my $o = $org // '') =~ s/ /_/g;
        next if $ex{$acc} || $ex{$b} || $ex{$o}; $n{$g}{$ord}{$p}++ }
    close $h;
    for my $g (keys %n) { for my $ord (sort keys %{ $n{$g} }) { my $c = $n{$g}{$ord}; my @p = sort { $c->{$b} <=> $c->{$a} || $a cmp $b } keys %$c; splice(@p, 3) if @p > 3; push @{ $refs{$g} }, map { [$ord, $_, ''] } @p } }
    \%refs }

# Smith-Waterman (match 2, mismatch -1, gap -2) of a peptide against a translated window: score, and the window
# index where the peptide's first residue falls (alignment start in the window minus the peptide residues skipped)
sub anchor { my @r = sw(@_); return @r[0 .. 2] }
# sw(peptide, window) -> (score, window index of the peptide's first residue, peptide residues skipped before the alignment,
#   peptide residues aligned through (its end index), window index just past the alignment)
sub sw { my ($g, $t) = @_; my ($best, $bi, $bj) = (0, 0, 0); my (@H, @P); my ($lg, $lt) = (length $g, length $t);
    $H[$_][0] = 0 for 0 .. $lg; $H[0][$_] = 0 for 0 .. $lt;
    for my $i (1 .. $lg) { my $gc = substr($g, $i - 1, 1); for my $j (1 .. $lt) { my $tc = substr($t, $j - 1, 1); my $m = ($gc eq $tc && $tc ne '*') ? 2 : -1;
        my ($d, $u, $l) = ($H[$i - 1][$j - 1] + $m, $H[$i - 1][$j] - 2, $H[$i][$j - 1] - 2); my ($v, $k) = (0, 0); ($v, $k) = ($d, 1) if $d > $v; ($v, $k) = ($u, 2) if $u > $v; ($v, $k) = ($l, 3) if $l > $v;
        $H[$i][$j] = $v; $P[$i][$j] = $k; ($best, $bi, $bj) = ($v, $i, $j) if $v > $best } }
    my ($i, $j) = ($bi, $bj); while ($i > 0 && $j > 0 && $H[$i][$j] > 0) { my $k = $P[$i][$j]; last unless $k; if ($k == 1) { $i--; $j-- } elsif ($k == 2) { $i-- } else { $j-- } }
    return ($best, $j - $i, $i, $bi, $bj) }

# anchors(refs, oriented x, number of codons of x to translate) -> [[codon position, order], ...] for the references
# whose own first residues align (score >= 20, at most 2 leading residues skipped)
sub anchors { my ($refs, $x, $nwin) = @_; $nwin = int(length($x) / 3) if $nwin > int(length($x) / 3); my $win = translate(substr($x, 0, 3 * $nwin)); my @a;
    for my $r (@{ $refs || [] }) { my ($sc, $at, $skip) = anchor($r->[1], $win); push @a, [$at, $r->[0]] if $sc >= 20 && $at >= 0 && $skip <= 2 } \@a }
# orders whose references put their first residue in [lo, hi] (codon positions)
sub orders_in { my ($an, $lo, $hi) = @_; my %o = map { $_->[1] => 1 } grep { $_->[0] >= $lo && $_->[0] <= $hi } @$an; my @r = sort keys %o; @r }   # a list; in scalar context the count (sort itself is undefined there)

# The two calls the evidence supports, both relative to the current start at codon k0 of x:
#   veto_extension(anchors, k0): true when at least two orders put the N-terminus within 3 codons of the call's start and
#     none puts it more than 5 codons upstream: the call already starts where distant lineages start, and a guide
#     majority wanting it longer is a shared convention (Liliaceae ndhK).
#   trim(anchors, x, k0, maxdown, acg, whole): when no order supports a start within 5 codons of the call's start or
#     upstream of it, and two or more orders agree on a position downstream, the in-frame ATG (else GTG; ACG for the
#     edited genes) within 3 codons of the most upstream such position, the frame open from there. -> {pos, codon, orders} or undef
sub veto_extension { my ($an, $k0) = @_; my @at = orders_in($an, $k0 - 3, $k0 + 3); my @up = orders_in($an, 0, $k0 - 6); return (@at >= 2 && !@up) ? join(',', @at) : '' }
sub trim { my ($an, $x, $k0, $maxdown, $acg, $whole) = @_;
    return undef if orders_in($an, 0, $k0 + 5);
    my %pos; $pos{ $_->[0] } = 1 for grep { $_->[0] > $k0 + 5 && $_->[0] <= $k0 + $maxdown + 3 } @$an;
    for my $p (sort { $a <=> $b } keys %pos) { my @o = orders_in($an, $p - 2, $p + 2); next unless @o >= 2;
        my @c = sort { abs($a - $p) <=> abs($b - $p) || $a <=> $b } grep { $_ > $k0 && $_ <= $k0 + $maxdown } ($p - 3 .. $p + 3);
        for my $want ('ATG', 'GTG', ($acg ? 'ACG' : ())) { for my $c (@c) { next unless substr($x, 3 * $c, 3) eq $want;
                my $pr = translate(substr($x, 3 * $c)); if ($whole) { $pr =~ s/\*$//; next if $pr =~ /\*/ } else { next if $pr =~ /\*/ }
                return { pos => $c, codon => $want, orders => join(',', @o) } } } }
    undef }

# near_start(anchors, x, k0, acg, whole): for a call on a rare initiator (TTG, ATT, ...), the ATG (ACG for the edited genes) within
# 3 codons of a position within 5 codons of the call's start where two or more orders start the protein, the frame open from
# there (a Lilium ndhK on an ATT two codons above the MVL... ATG). -> {pos, codon, orders} or undef
sub near_start { my ($an, $x, $k0, $acg, $whole) = @_; my %pos; $pos{ $_->[0] } = 1 for grep { abs($_->[0] - $k0) <= 5 } @$an;
    for my $p (sort { abs($a - $k0) <=> abs($b - $k0) || $a <=> $b } keys %pos) { my @o = orders_in($an, $p - 2, $p + 2); next unless @o >= 2;
        for my $want ('ATG', ($acg ? 'ACG' : ())) { for my $c (sort { abs($a - $p) <=> abs($b - $p) || $a <=> $b } grep { $_ >= 0 && $_ != $k0 } ($p - 3 .. $p + 3)) { next unless substr($x, 3 * $c, 3) eq $want;
                my $pr = translate(substr($x, 3 * $c)); if ($whole) { $pr =~ s/\*$//; next if $pr =~ /\*/ } else { next if $pr =~ /\*/ }
                return { pos => $c, codon => $want, orders => join(',', @o) } } } }
    undef }

# extend(anchors, x, k0, acg): the mirror of trim, for a call that starts inside the protein: no order starts the protein
# within 5 codons of the call's start, two or more orders start it at one position upstream; the ATG (else GTG; ACG for the
# edited genes) within 3 codons of the nearest such position, the frame open from there to the call (a Lilium ycf1 called at
# 2844 of ~5550 nt, where only one guide carries ycf1 because most Liliaceae RefSeq records keep it as a pseudogene).
# -> {pos, codon, orders} or undef
sub extend { my ($an, $x, $k0, $acg) = @_;
    return undef if orders_in($an, $k0 - 5, $k0 + 5);
    my %pos; $pos{ $_->[0] } = 1 for grep { $_->[0] < $k0 - 5 } @$an;
    for my $p (sort { $b <=> $a } keys %pos) { my @o = orders_in($an, $p - 2, $p + 2); next unless @o >= 2;
        my @c = sort { abs($a - $p) <=> abs($b - $p) || $a <=> $b } grep { $_ >= 0 && $_ < $k0 } ($p - 3 .. $p + 3);
        for my $want ('ATG', 'GTG', ($acg ? 'ACG' : ())) { for my $c (@c) { next unless substr($x, 3 * $c, 3) eq $want;
                next if translate(substr($x, 3 * $c, 3 * ($k0 - $c))) =~ /\*/;
                return { pos => $c, codon => $want, orders => join(',', @o) } } } }
    undef }


# ---- splice junctions ---------------------------------------------------------------------------------------------
# A junction a few nucleotides off changes the protein AT the junction (petB exon 2 annotated 642, 648 or 657 nt long
# in the Liliaceae adds or drops residues after MS|), so the proteins of other orders can place it. tables/junction_refs.tsv
# (refdb_build/build_junction_refs.pl) holds, for every cis-splice junction of the curated genomes, the 15 whole codons
# before it, the 15 after it and the phase (nucleotides of the split codon in the upstream exon).
#
# load_junctions(file, exclude list) -> { gene => { junction number => [[order, phase, left, right], ...] } }, the three
# commonest (phase, left, right) per order
sub load_junctions { my ($file, $excl) = @_; my %ex; for (split /,/, ($excl // '')) { next unless length; $ex{$_} = 1; (my $b = $_) =~ s/\.\d+$//; $ex{$b} = 1 }
    my (%refs, %n); open my $h, "<", $file or return {};
    while (<$h>) { chomp; my ($g, $j, $acc, $ord, $org, $ph, $L, $R) = split /\t/; next if $g eq 'gene' || !defined $R; (my $b = $acc) =~ s/\.\d+$//; (my $o = $org // '') =~ s/ /_/g;
        next if $ex{$acc} || $ex{$b} || $ex{$o}; $n{$g}{$j}{$ord}{"$ph\t$L\t$R"}++ }
    close $h;
    for my $g (keys %n) { for my $j (keys %{ $n{$g} }) { for my $ord (sort keys %{ $n{$g}{$j} }) { my $c = $n{$g}{$j}{$ord}; my @k = sort { $c->{$b} <=> $c->{$a} || $a cmp $b } keys %$c; splice(@k, 3) if @k > 3;
        push @{ $refs{$g}{$j} }, map { [$ord, split /\t/, $_, -1] } @k } } }
    \%refs }
sub oseq { my ($seq, $lo, $hi, $d) = @_; return '' if $lo < 1 || $hi > length($$seq) || $hi < $lo; my $x = substr($$seq, $lo - 1, $hi - $lo + 1); if ($d eq '-') { $x = reverse $x; $x =~ tr/ACGT/TGCA/ } $x }

# fix_junctions(junction refs of the gene, \genome, strand, exons in transcript order [[s, e], ...])
#   -> (new exon list, [[junction number, old donor, old acceptor, new donor, new acceptor, orders], ...])
# For each junction, every reference says where it falls in the target: its left peptide is aligned to the upstream exon
# in that exon's frame (fixed by the exons before it) and must align through its LAST residue; its right peptide is aligned
# to the downstream exon in the frame the exons after it fix and must align from its FIRST residue; the phase then gives
# the donor and acceptor nucleotides. A reference whose junction-adjacent residue does not match abstains, so an indel at
# the junction in the target's lineage is never extrapolated across. A first exon too short to align (petB, 6 nt) keeps its
# donor when its phase matches. The junction moves when two or more orders agree on one (donor, acceptor) pair, more orders
# than any other pair, the intron stays >= 50 nt and the spliced CDS stays intact.
sub fix_junctions { my ($jr, $seq, $d, $exons) = @_; my @ex = map { [@$_] } @$exons; my @chg; return (\@ex, \@chg) unless $jr && @ex >= 2;
    my $slen = length $$seq;
    for my $k (0 .. $#ex - 1) { my $refs = $jr->{ $k + 1 } or next; my ($u, $w) = ($ex[$k], $ex[$k + 1]);
        my $Ok = 0; $Ok += $ex[$_][1] - $ex[$_][0] + 1 for 0 .. $k - 1; my $Lk = $u->[1] - $u->[0] + 1; my $fk = (3 - $Ok % 3) % 3;
        my $Taft = 0; $Taft += $ex[$_][1] - $ex[$_][0] + 1 for $k + 2 .. $#ex; my $r = (3 - $Taft % 3) % 3;
        # upstream exon and 45 nt of the intron; downstream exon with 45 nt of the intron before it
        my ($glo, $ghi) = $d eq '+' ? ($u->[0], $u->[1] + 45) : ($u->[0] - 45, $u->[1]); my $G = oseq($seq, $glo, $ghi, $d); next unless $G;
        my ($hlo, $hhi) = $d eq '+' ? ($w->[0] - 45, $w->[1]) : ($w->[0], $w->[1] + 45); my $H = oseq($seq, $hlo, $hhi, $d); next unless $H;
        my $c0 = int(($Lk - $fk) / 3) - 30; $c0 = 0 if $c0 < 0; my $wd = translate(substr($G, $fk + 3 * $c0));
        my $M = length $H; my $h0 = ($M - $r) % 3; my $wa = translate(substr($H, $h0, 135));
        my $cur_don = $d eq '+' ? $u->[1] : $u->[0]; my $cur_acc = $d eq '+' ? $w->[0] : $w->[1]; my (%vote, %seen);
        for my $rf (@$refs) { my ($ord, $ph, $L, $R) = @$rf; next unless length($R) >= 8; my ($don, $acc);
            if (length($L) >= 8) { my ($sc, undef, undef, $bi, $bj) = sw($L, $wd); next unless $sc >= 16 && $bi == length $L;
                my $o = $fk + 3 * ($c0 + $bj) + $ph - 1; $don = $d eq '+' ? $glo + $o : $ghi - $o }
            # a first exon too short to align (petB 6 nt, petD 8, rpl16 9): when the target's exon begins with the reference's
            # own codons, the reference's exon length gives the donor (a petB exon 1 called 4 nt long becomes 6); otherwise
            # the call's donor stands if its phase matches
            elsif ($k == 0 && length($L) && translate(substr($G, 0, 3 * length $L)) eq $L) { my $o = 3 * length($L) + $ph - 1; $don = $d eq '+' ? $glo + $o : $ghi - $o }
            else { next unless ($Ok + $Lk) % 3 == $ph; $don = $cur_don }
            my ($sc, $at, $skip) = sw($R, $wa); next unless $sc >= 16 && $skip == 0 && $at >= 0;
            my $o = $h0 + 3 * $at - (3 - $ph) % 3; next if $o < 0; $acc = $d eq '+' ? $hlo + $o : $hhi - $o;
            next if $seen{"$ord\t$don|$acc"}++; $vote{"$don|$acc"}{$ord} = 1 }
        my @keys = sort { keys %{ $vote{$b} } <=> keys %{ $vote{$a} } || $a cmp $b } keys %vote; next unless @keys;
        my $n1 = keys %{ $vote{ $keys[0] } }; my $n2 = @keys > 1 ? keys %{ $vote{ $keys[1] } } : 0; next unless $n1 >= 2 && $n1 > $n2;
        my ($don, $acc) = split /\|/, $keys[0]; next if $don == $cur_don && $acc == $cur_acc;
        my @new = map { [@$_] } @ex; if ($d eq '+') { $new[$k][1] = $don; $new[$k + 1][0] = $acc } else { $new[$k][0] = $don; $new[$k + 1][1] = $acc }
        next if $new[$k][1] < $new[$k][0] || $new[$k + 1][1] < $new[$k + 1][0]; next if ($d eq '+' ? $acc - $don - 1 : $don - $acc - 1) < 50;
        my $cds = join '', map { oseq($seq, $_->[0], $_->[1], $d) } @new; next if length($cds) % 3; my $p = translate($cds); next unless substr($p, -1) eq '*'; $p =~ s/\*$//; next if $p =~ /\*/;
        push @chg, [$k + 1, $cur_don, $cur_acc, $don, $acc, join(',', sort keys %{ $vote{ $keys[0] } })]; @ex = @new }
    return (\@ex, \@chg) }


# ---- frameshifts at homopolymer runs ------------------------------------------------------------------------------
# homopolymer_site(cds) -> (offset of the run in the CDS, its length, its base, "the lengths that would restore the frame")
# For a span that runs from a start codon to a stop but is not one open frame: the run of >= 7 identical bases whose
# length, changed by one or two, makes the whole span an intact CDS (a multiple of 3, one stop, at the end). That run is
# where the frame is lost; everything before it and everything after it is in frame. The longest such run; () when no
# run does it (the frameshift is not at a homopolymer) or the span is already a multiple of 3.
sub homopolymer_site { my $x = shift; return () unless length($x) % 3; my @best;
    while ($x =~ /(A{7,}|C{7,}|G{7,}|T{7,})/g) { my ($o, $l, $b) = ($-[0], length $1, substr($1, 0, 1)); my @ok;
        for my $dl (-2, -1, 1, 2) { next if (length($x) + $dl) % 3; my $y = $x; substr($y, $o, $l) = $b x ($l + $dl); my $p = translate($y); next unless substr($p, -1) eq '*'; chop $p; next if index($p, '*') >= 0; push @ok, $l + $dl }
        @best = ($o, $l, $b, join(' or ', @ok)) if @ok && (!@best || $l > $best[1]) }
    @best }
1;
