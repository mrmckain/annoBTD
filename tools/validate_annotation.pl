#!/usr/bin/perl
# validate_annotation.pl - one verdict per annotated feature.
#
# USAGE: validate_annotation.pl <annotation.txt> <genome.fsa> --guides-dir DIR --guide-seqs DIR
#            [--refdb DIR] [--family F] [--order O] [--genus G] [--val table2asn.val]
#            [--out report.tsv] [--blast-bin DIR] [--window N]
#
# Every CDS, tRNA and rRNA row of the annotation is tested against what GenBank
# leads us to expect: the lineage profile (refdb/gene_expect.tsv), the guide
# genomes' copies of the gene (exon count, length, protein), the translation
# itself, the neighbouring genes and the IR boundaries. A feature that passes is
# PASS. A feature that fails is re-derived from the guide protein (tblastn) and
# the verdict says what kind of failure it is:
#
#   ARTIFACT    the sequence carries a canonical form the annotator missed
#               (wrong start, boundary one codon off, a remnant copy at the IR
#               boundary); the proposed edit is given in verdant_edit.pl format
#   REAL        no canonical form exists: the guide protein aligns in two frames
#               or stops early, or part of the gene is gone; a pseudogene call is
#               proposed with the evidence
#   UNRESOLVED  ambiguity codes in the span, an IR boundary inside the gene, or a
#               spliced gene whose exons cannot be placed; look at the assembly
#
# Canonical genes the guides carry but the annotation lacks are tested the same
# way. With --val, table2asn's validator messages are attached to the rows they
# concern, so GenBank's view and ours sit side by side.
#
# OUT (tab-separated): feature  start  end  strand  status  flags  evidence  edit  table2asn
#   edit: one or more verdant_edit.pl rows (without the genome column, fields separated by spaces), joined by " ; "
use strict; use warnings; use FindBin;

my ($ann, $fsa, @rest) = @ARGV;
die "usage: $0 <annotation.txt> <genome.fsa> --guides-dir DIR --guide-seqs DIR [--refdb DIR] [--family F --order O --genus G] [--val FILE] [--out FILE]\n" unless $fsa;
my %o = (window => 3000);
while (@rest) { my $k = shift @rest; $k =~ s/^--//; $o{$k} = shift @rest }
my $ANNOBTD = $ENV{ANNOBTD_DIR} || "$FindBin::Bin/..";
$o{refdb} //= $ENV{ANNOBTD_REFDB} || "$ANNOBTD/refdb";
$o{family} //= $ENV{ANNOBTD_FAMILY}; $o{order} //= $ENV{ANNOBTD_ORDER}; $o{genus} //= $ENV{ANNOBTD_GENUS};
my $out; if ($o{out}) { open $out, ">", $o{out} or die "$o{out}: $!" } else { $out = \*STDOUT }

# ---- genome ---------------------------------------------------------------------
my $seq = ''; open my $fh, "<", $fsa or die "$fsa: $!"; while (<$fh>) { next if /^>/; s/\s//g; $seq .= uc $_ } close $fh;
my $len = length $seq; die "no sequence in $fsa\n" unless $len;
my %CODON; { my @b = qw(T C A G); my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; my $i = 0; for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1) } } } }
# Start codons. ATG and GTG are accepted everywhere; ACG only for the genes whose
# C-to-U editing of the start is documented (ndhD, psbL, rpl2, psbC). The table-11
# initiators (TTG CTG ATT ATC ATA) and ACG elsewhere are flagged: an open frame that
# begins at one of them two codons from an out-of-frame ATG is how a frameshifted
# gene passes as intact (Schoenolirion rps19, Hosta infA).
my %INIT = map { $_ => 1 } qw(ATG GTG); my %EDITABLE = map { $_ => 1 } qw(ndhD psbL rpl2 psbC);
sub init_ok { my ($c, $g) = @_; $INIT{$c} || ($c eq 'ACG' && $g && $EDITABLE{$g}) }
sub translate { my $s = shift; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr($s, $i, 3)} // 'X' } $p }
sub rc { my $x = reverse shift; $x =~ tr/ACGTRYKMSWN/TGCAYRMKSWN/; $x }
sub spliced { my ($d, @ex) = @_; my $x = ''; for my $e (sort { $a->[0] <=> $b->[0] } @ex) { $x .= substr($seq, $e->[0] - 1, $e->[1] - $e->[0] + 1) } $d eq '-' ? rc($x) : $x }
sub intact_info { my ($x, $g) = @_; my %r = (len => length $x, mod3 => length($x) % 3, start => substr($x, 0, 3)); my $p = translate($x); $r{term} = substr($p, -1) eq '*' ? 1 : 0; (my $body = $p) =~ s/\*$//; $r{internal} = ($body =~ tr/*//); $r{X} = ($p =~ tr/X//); $r{intact} = (!$r{mod3} && init_ok($r{start}, $g) && $r{term} && !$r{internal}) ? 1 : 0; $r{prot} = $p; \%r }
my $blast_genome; sub blast_genome { return $blast_genome if $blast_genome; require File::Temp; my ($h, $f) = File::Temp::tempfile(UNLINK => 1); (my $s = $seq) =~ tr/ACGT/N/c; print $h ">genome\n$s\n"; close $h; $blast_genome = $f }
my $tblastn = $o{'blast-bin'} ? "$o{'blast-bin'}/tblastn" : $ENV{BLAST_BIN} ? "$ENV{BLAST_BIN}/tblastn" : (grep { -x $_ } map { "$_/tblastn" } split /:/, $ENV{PATH})[0];
die "no tblastn (set --blast-bin or BLAST_BIN)\n" unless $tblastn && -x $tblastn;

# ---- canonical names, lineage expectations -----------------------------------------
my %CANON; { my $f = "$ANNOBTD/tables/plastid_cds_names.txt"; if (open my $ch, "<", $f) { while (<$ch>) { next if /^#/; $CANON{$_} = 1 for split } close $ch } }
my %cell;   # gene -> type(F/G/O) -> [n, median, cons, alt, altfrac]
if (open my $eh, "<", "$o{refdb}/gene_expect.tsv") { while (<$eh>) { chomp; my @f = split /\t/; next if $f[0] eq 'name';
    my $t = $f[1]; my $want = $t eq 'F' ? $o{family} : $t eq 'G' ? $o{genus} : $t eq 'O' ? $o{order} : undef; next unless $want && $f[2] eq $want; $cell{$f[0]}{$t} = [@f[3..7]] } close $eh }
# the settled expectation: family first, then genus, then order, each with n >= 10 and consensus >= 0.9 (the post filter's rule)
sub expect { my $g = shift; for my $t (qw(F G O)) { my $c = $cell{$g}{$t} or next; return ($c->[1], "$t cell n=$c->[0] cons=$c->[2]") if $c->[0] >= 10 && $c->[2] >= 0.9 } return () }
sub profile_lengths { my $g = shift; my (@l, %seen); for my $t (qw(G F O)) { my $c = $cell{$g}{$t} or next; push @l, $c->[1] unless $seen{$c->[1]}++ } for my $t (qw(G F O)) { my $c = $cell{$g}{$t} or next; push @l, $c->[3] if $c->[3] && $c->[4] >= 0.1 && !$seen{$c->[3]}++ } @l }

# ---- guides: exon structure and lengths per gene, sequences for proteins --------------
my (%gd, $nguides, %gseq);   # gd{gene}{guide} = { d, ex => [[s,e],...], len, nex }
die "--guides-dir and --guide-seqs are required\n" unless $o{'guides-dir'} && $o{'guide-seqs'};
opendir(my $gdir, $o{'guides-dir'}) or die "$o{'guides-dir'}: $!";
for my $g (sort grep { !/^\./ } readdir $gdir) { open my $gh, "<", "$o{'guides-dir'}/$g" or next; $nguides++;
    my %rows; while (<$gh>) { chomp; my @f = split /\t/; next unless @f >= 4 && $f[1] =~ /^\d+$/; next if $f[0] =~ /~|_intron/;
        (my $base = $f[0]) =~ s/_exon\d+$//; my $k = "$base\t$f[3]"; push @{ $rows{$k} }, [$f[1], $f[2], $f[0], $f[4] // ''] } close $gh;
    for my $k (keys %rows) { my ($base, $d) = split /\t/, $k; my @ex = sort { $a->[0] <=> $b->[0] } @{ $rows{$k} };
        my @groups = ([shift @ex]); for my $x (@ex) { if ($x->[0] - $groups[-1][-1][1] > 5000) { push @groups, [$x] } else { push @{ $groups[-1] }, $x } }
        my ($best) = sort { @$b <=> @$a } @groups;   # the IR gives two copies; one is enough
        my $L = 0; $L += $_->[1] - $_->[0] + 1 for @$best;
        my (@alll, @starts); my $gs = guide_seq($g);
        for my $grp (@groups) { my $l = 0; $l += $_->[1] - $_->[0] + 1 for @$grp; push @alll, $l;
            if ($gs && $base !~ /^(trn|rrn)/) { my ($s0, $e0) = ($grp->[0][0], $grp->[-1][1]); push @starts, $d eq '-' ? rc(substr($gs, $e0 - 3, 3)) : substr($gs, $s0 - 1, 3) } }
        # the other strand's copies: keep their lengths and starts, and the longer copy as the guide's own (a ycf1
        # straddling the IR has a full copy on one strand and a fragment on the other; hash order used to pick one)
        if (my $old = $gd{$base}{$g}) { push @{ $old->{alllens} }, @alll; push @{ $old->{allstarts} }, @starts;
            @$old{qw(d ex len nex type)} = ($d, [map { [$_->[0], $_->[1]] } @$best], $L, scalar @$best, $best->[0][3]) if $L > $old->{len} || ($L == $old->{len} && $d eq '+' && $old->{d} eq '-'); next }
        $gd{$base}{$g} = { d => $d, ex => [map { [$_->[0], $_->[1]] } @$best], len => $L, nex => scalar @$best, type => $best->[0][3], alllens => \@alll, allstarts => \@starts } } }
closedir $gdir;
sub guide_seq { my $g = shift; return $gseq{$g} if exists $gseq{$g}; my $s = ''; if (open my $h, "<", "$o{'guide-seqs'}/$g") { while (<$h>) { next if /^>/; s/\s//g; $s .= uc $_ } close $h } $gseq{$g} = $s }
sub mode { my %c; $c{$_}++ for @_; my ($m) = sort { $c{$b} <=> $c{$a} || $a <=> $b } keys %c; $m }
sub carriers { my $g = shift; $gd{$g} ? scalar keys %{ $gd{$g} } : 0 }
sub guide_len  { my $g = shift; $gd{$g} ? mode(map { $_->{len} } values %{ $gd{$g} }) : undef }
sub guide_nex  { my $g = shift; $gd{$g} ? mode(map { $_->{nex} } values %{ $gd{$g} }) : undef }
my %gprot;
# Smith-Waterman (match 2, mismatch -1, gap -2) of a guide's N-terminal peptide against a translated window of
# the target: the score, and where the guide's first residue falls in the window (as in post_filter_annotation.pl)
sub nterm_anchor { my ($g, $t) = @_; my ($best, $bi, $bj) = (0, 0, 0); my (@H, @P);
    for my $i (0 .. length $g) { $H[$i][0] = 0 } for my $j (0 .. length $t) { $H[0][$j] = 0 }
    for my $i (1 .. length $g) { for my $j (1 .. length $t) { my $m = substr($g, $i - 1, 1) eq substr($t, $j - 1, 1) && substr($t, $j - 1, 1) ne '*' ? 2 : -1;
        my @o = (0, $H[$i - 1][$j - 1] + $m, $H[$i - 1][$j] - 2, $H[$i][$j - 1] - 2); my $k = 0; for (1 .. 3) { $k = $_ if $o[$_] > $o[$k] }
        $H[$i][$j] = $o[$k]; $P[$i][$j] = $k; ($best, $bi, $bj) = ($o[$k], $i, $j) if $o[$k] > $best } }
    my ($i, $j) = ($bi, $bj); while ($i > 0 && $j > 0 && $H[$i][$j] > 0) { my $k = $P[$i][$j]; last unless $k; if ($k == 1) { $i--; $j-- } elsif ($k == 2) { $i-- } else { $j-- } }
    return ($best, $j - $i) }
# the N-terminal 30 aa of every guide copy of gene g that begins ATG or GTG
sub guide_nterms { my $g = shift; my @p; for my $name (sort keys %{ $gd{$g} || {} }) { my $r = $gd{$g}{$name}; my $gs = guide_seq($name) or next;
        my $x = ''; for my $e (sort { $a->[0] <=> $b->[0] } @{ $r->{ex} }) { $x .= substr($gs, $e->[0] - 1, $e->[1] - $e->[0] + 1) } $x = rc($x) if $r->{d} eq '-';
        my $c = substr($x, 0, 3); push @p, 'M' . substr(translate(substr($x, 3, 87)), 0, 29) if $c eq 'ATG' || $c eq 'GTG' } @p }
# A start put out of frame by a homopolymer run (post_filter_annotation.pl section 4b; Liliaceae cemA,
# an ATG followed by 9 to 25 A). For a call in frame with its stop: (new s, new e, length, run length,
# run base, anchors) when the guides' N-termini put the start where the call's frame has no ATG/GTG,
# an ATG sits out of frame within 6 nt, a run of >= 7 identical bases starts within 45 nt after it,
# and the whole span is within 9 nt of a guide copy's length.
sub hp_shift { my ($s, $e, $d, $g) = @_; return () if ($e - $s + 1) % 3; my @peps = guide_nterms($g); return () unless @peps && 2 * @peps >= carriers($g);
    my $avail = $d eq '+' ? $s - 1 : $len - $e; my $up = $avail < 450 ? $avail - $avail % 3 : 450; return () if $up < 30;
    my ($xs, $xe) = $d eq '+' ? ($s - $up, $e) : ($s, $e + $up); my $x = spliced($d, [$xs, $xe]);
    return () unless substr($x, 0, $up + 90) =~ /A{7}|C{7}|G{7}|T{7}/;
    my $win = translate(substr($x, 0, $up + 90)); my @anch; for my $gp (@peps) { my ($sc, $at) = nterm_anchor($gp, $win); push @anch, $at if $sc >= 20 } return () unless @anch;
    my $A = (sort { $a <=> $b } @anch)[$#anch / 2];
    for my $k ($A - 3 .. $A + 3) { next if $k < 0; return () if init_ok(substr($x, 3 * $k, 3), $g) }
    # The ATG is found from the run, not from the anchor's exact position: a lysine run aligns at any offset, so the
    # guides' first residue lands up to several codons off when the target's run is longer or shorter than theirs (Lilium
    # fargesii cemA: 19 A against the guides' 9-15 put it 7 nt from the ATG, outside a +-6 nt search). Any out-of-frame ATG
    # within 15 codons of the anchor that a run of >= 7 identical bases follows within 45 nt qualifies; the nearest wins.
    my ($o, $run, $rp);
    for my $c (sort { abs($a - 3 * $A) <=> abs($b - 3 * $A) || $a <=> $b } grep { $_ >= 0 && ($_ - 3 * $A) % 3 } (3 * $A - 45 .. 3 * $A + 45)) {
        next unless substr($x, $c, 3) eq 'ATG' && substr($x, $c + 3, 52) =~ /(A{7,}|C{7,}|G{7,}|T{7,})/ && $-[0] <= 45; ($o, $run, $rp) = ($c, $1, $-[0]); last }
    return () unless defined $o;
    my $NL = length($x) - $o; my $tol = 9 + (length($run) > 12 ? length($run) - 12 : 0);
    return () unless grep { my $c = $_; grep { abs($NL - $_) <= $tol } @{ $c->{alllens} } } values %{ $gd{$g} };
    my ($ns, $ne) = $d eq '+' ? ($xs + $o, $e) : ($s, $xe - $o);
    return ($ns, $ne, $NL, length $run, substr($run, 0, 1), scalar @anch) }
# The same gene already called whole: an ATG start, a length that is not a multiple of 3, a run of >= 7
# within 45 nt after the ATG, the body open in the stop's frame after the run, a guide copy's length.
# Returns (run length, run base, run position in the CDS).
sub hp_framed { my ($x, $g) = @_; my $L = length $x; return () unless $L % 3 && substr($x, 0, 3) eq 'ATG';
    return () unless substr($x, 3, 52) =~ /(A{7,}|C{7,}|G{7,}|T{7,})/; my ($run, $rp) = ($1, $-[0] + 3); return () if $rp > 48;
    my $p = $rp + length $run; $p++ while ($L - $p) % 3; my $body = translate(substr($x, $p)); return () unless substr($body, -1) eq '*' && index($body, '*') == length($body) - 1;
    return () unless $gd{$g} && grep { my $c = $_; grep { abs($L - $_) <= 9 } @{ $c->{alllens} } } values %{ $gd{$g} };
    return (length $run, substr($run, 0, 1), $rp + 1) }
sub guide_protein { my $g = shift; return @{ $gprot{$g} } if $gprot{$g}; for my $name (sort keys %{ $gd{$g} || {} }) { my $r = $gd{$g}{$name}; my $gs = guide_seq($name) or next;
        my $x = ''; for my $e (sort { $a->[0] <=> $b->[0] } @{ $r->{ex} }) { $x .= substr($gs, $e->[0] - 1, $e->[1] - $e->[0] + 1) } $x = rc($x) if $r->{d} eq '-';
        my $p = translate($x); $p =~ s/\*$//; next if $p =~ /\*/ || length($p) < 20; $gprot{$g} = [$name, $p, $r->{len}]; return @{ $gprot{$g} } } $gprot{$g} = []; () }

# ---- the annotation -----------------------------------------------------------------
my (@rows, @ir, %region);
open my $ah, "<", $ann or die "$ann: $!"; while (<$ah>) { chomp; next unless /\S/; my ($n, $s, $e, $d) = split /\t/; next unless defined $d && $s =~ /^\d+$/ && $e =~ /^\d+$/; ($s, $e) = ($e, $s) if $e < $s;
    if ($n =~ /^(LSC|SSC|IRA|IRB)$/) { $region{$n} = [$s, $e]; next } next if $n =~ /~/ || $n eq 'FULL'; push @rows, [$n, $s, $e, $d] } close $ah;
for my $r (qw(IRA IRB)) { push @ir, @{ $region{$r} } if $region{$r} }
@ir = sort { $a <=> $b } @ir;
sub near_ir { my ($s, $e, $pad) = @_; for my $b (@ir) { return $b if $b >= $s - $pad && $b <= $e + $pad } 0 }
sub ambiguous { my ($s, $e) = @_; $s = 1 if $s < 1; $e = $len if $e > $len; my $x = substr($seq, $s - 1, $e - $s + 1); my $n = ($x =~ tr/ACGT//c); $n }
# group CDS rows into genes (exons by name, IR copies by distance)
my (%feat, @order); our @rows_all = @rows;   # key gene\tstrand\tcopy -> { gene, d, ex => [[s,e,name]], pseudo }
{ my %by; for my $r (@rows) { my ($n, $s, $e, $d) = @$r; next if $n =~ /_intron\d*$/; (my $base = $n) =~ s/_exon\d+$//; push @{ $by{"$base\t$d"} }, [$s, $e, $n] }
  for my $k (sort keys %by) { my ($base, $d) = split /\t/, $k; my @ex = sort { $a->[0] <=> $b->[0] } @{ $by{$k} };
      my @groups = ([shift @ex]); for my $x (@ex) { if ($x->[0] - $groups[-1][-1][1] > 5000) { push @groups, [$x] } else { push @{ $groups[-1] }, $x } }
      my $c = 0; for my $grp (@groups) { $c++; my $key = "$base\t$d\t$c"; $feat{$key} = { gene => $base, d => $d, ex => $grp, pseudo => ($base =~ /_pseudogene$/ ? 1 : 0), rna => ($base =~ /^(trn|rrn)/ ? 1 : 0) }; push @order, $key } } }
my %present; for my $k (keys %feat) { (my $g = $feat{$k}{gene}) =~ s/_pseudogene$//; $present{$g} = 1 }
# same-strand CDS overlap of [s,e] with other genes
# an overlap the guides carry themselves: at least half of the guides with both genes overlap them by >= ov - 15 nt
# (ndhK/ndhC by about 121 nt in the Liliaceae and other lineages, psbD/psbC by 53, atpE/atpB by 4)
sub guides_overlap { my ($g, $who, $ov) = @_; my ($n, $ok) = (0, 0); for my $name (keys %{ $gd{$g} || {} }) { my $b = $gd{$who}{$name} or next; my $a = $gd{$g}{$name}; next unless $a->{d} eq $b->{d};
        my ($as, $ae) = ($a->{ex}[0][0], $a->{ex}[-1][1]); my ($bs, $be) = ($b->{ex}[0][0], $b->{ex}[-1][1]); my $lo = $as > $bs ? $as : $bs; my $hi = $ae < $be ? $ae : $be; $n++; $ok++ if $hi - $lo + 1 >= $ov - 15 }
    return $n && 2 * $ok >= $n ? 1 : 0 }
sub overlap_other { my ($key, $s, $e, $d) = @_; my ($best, $who) = (0, ''); for my $k (keys %feat) { next if $k eq $key; my $f = $feat{$k}; next if $f->{rna} || $f->{pseudo} || $f->{d} ne $d; for my $x (@{ $f->{ex} }) { my $lo = $s > $x->[0] ? $s : $x->[0]; my $hi = $e < $x->[1] ? $e : $x->[1]; if ($hi - $lo + 1 > $best) { $best = $hi - $lo + 1; $who = $f->{gene} } } } ($best, $who) }

# ---- tblastn of the guide protein, HSPs restricted to a window --------------------------
my $tmpdir; sub hsps { my ($g, $s, $e, $d) = @_; my ($gname, $prot, $glen) = guide_protein($g); return () unless $prot;
    require File::Temp; $tmpdir //= File::Temp::tempdir(CLEANUP => 1); open my $q, ">", "$tmpdir/q.fa" or die; print $q ">$g\n$prot\n"; close $q;
    my @h = map { chomp; my @f = split /\t/; { ss => $f[0] < $f[1] ? $f[0] : $f[1], se => $f[0] < $f[1] ? $f[1] : $f[0], frame => $f[2], pid => $f[3], alen => $f[4], qs => $f[5], qe => $f[6], ev => $f[7], d => ($f[2] < 0 ? '-' : '+') } }
        `$tblastn -query $tmpdir/q.fa -subject @{[blast_genome()]} -seg no -outfmt '6 sstart send sframe pident length qstart qend evalue' -evalue 1e-3 2>/dev/null`;
    @h = grep { $_->{pid} >= 50 && $_->{alen} >= 12 } @h;
    if (defined $s) { @h = grep { $_->{d} eq $d && $_->{se} >= $s - $o{window} && $_->{ss} <= $e + $o{window} } @h }
    return (\@h, length $prot, $gname, $glen) }
sub coverage { my ($h, $plen) = @_; my @c = (0) x ($plen + 1); for my $x (@$h) { $c[$_] = 1 for $x->{qs} .. $x->{qe} } my $n = 0; $n += $_ for @c; $plen ? $n / $plen : 0 }
# the open frame that a single-frame hit sits in: downstream to the stop, upstream to the initiator nearest the expected length
sub orf_from_hit { my ($hs, $he, $d, $want, $g) = @_; my ($s, $e);
    if ($d eq '+') { $e = $he; while ($e + 3 <= $len) { my $c = substr($seq, $e, 3); $e += 3; last if ($CODON{$c} // '') eq '*' } my @cand; my $p = $hs; my $c0 = substr($seq, $hs - 1, 3); push @cand, $hs if init_ok($c0, $g);
        while ($p >= 4 && $p > $hs - 600) { $p -= 3; my $c = substr($seq, $p - 1, 3); last if ($CODON{$c} // '') eq '*'; push @cand, $p if init_ok($c, $g) }
        $p = $hs; while ($p + 3 < $he - 60) { $p += 3; my $c = substr($seq, $p - 1, 3); push @cand, $p if init_ok($c, $g) }   # a start inside the hit, for an over-extended call
        return () unless @cand; ($s) = sort { abs(($e - $a + 1) - $want) <=> abs(($e - $b + 1) - $want) } @cand }
    else { $s = $hs; while ($s - 3 >= 1) { my $c = rc(substr($seq, $s - 4, 3)); $s -= 3; last if ($CODON{$c} // '') eq '*' } my @cand; my $p = $he; my $c0 = rc(substr($seq, $he - 3, 3)); push @cand, $he if init_ok($c0, $g);
        while ($p + 3 <= $len && $p < $he + 600) { $p += 3; my $c = rc(substr($seq, $p - 3, 3)); last if ($CODON{$c} // '') eq '*'; push @cand, $p if init_ok($c, $g) }
        $p = $he; while ($p - 3 > $hs + 60) { $p -= 3; my $c = rc(substr($seq, $p - 3, 3)); push @cand, $p if init_ok($c, $g) }
        return () unless @cand; ($e) = sort { abs(($a - $s + 1) - $want) <=> abs(($b - $s + 1) - $want) } @cand }
    ($s, $e) }
# an ATG/GTG on strand d within 15 nt of where the guide's first residue maps, in a frame other than the hit's
sub offframe_start { my ($h, $d) = @_; my ($first) = sort { $a->{qs} <=> $b->{qs} } @$h; return () unless $first;
    my $pos = $d eq '+' ? $first->{ss} - 3 * ($first->{qs} - 1) : $first->{se} + 3 * ($first->{qs} - 1);
    for my $k (-15 .. 15) { next unless $k % 3; my $p = $pos + $k; next if $p < 1 || $p + 2 > $len; my $c = $d eq '+' ? substr($seq, $p - 1, 3) : rc(substr($seq, $p - 3, 3)); return ($p, $c) if $INIT{$c} } () }
# A spliced gene with fewer exons than the guides: find the missing exon(s) by blastn of the
# guide's exon sequences in the window and propose the rows (exons renumbered in coding order,
# intron rows between). Returns (edits, exon rows) or () when an exon cannot be placed.
my $blastn = $tblastn; $blastn =~ s/tblastn$/blastn/;
sub place_missing_exons { my ($g, $gname, $d, $s, $e, $ex, $want) = @_; my $gr = $gd{$g}{$gname} or return (); my $gs = guide_seq($gname) or return ();
    my @gex = sort { $a->[0] <=> $b->[0] } @{ $gr->{ex} }; @gex = reverse @gex if $gr->{d} eq '-';   # coding order
    require File::Temp; $tmpdir //= File::Temp::tempdir(CLEANUP => 1); my @placed;
    for my $k (0 .. $#gex) { my ($gs0, $ge0) = @{ $gex[$k] }; my $x = substr($gs, $gs0 - 1, $ge0 - $gs0 + 1); $x = rc($x) if $gr->{d} eq '-'; return () if length($x) < 20;
        open my $q, ">", "$tmpdir/ex.fa" or die; print $q ">e$k\n$x\n"; close $q;
        my @h = map { chomp; [split /\t/] } `$blastn -query $tmpdir/ex.fa -subject @{[blast_genome()]} -outfmt '6 qstart qend sstart send length pident sstrand' -evalue 1e-10 2>/dev/null`;
        my ($best) = sort { $b->[4] <=> $a->[4] } grep { my ($a1, $a2) = sort { $a <=> $b } ($_->[2], $_->[3]); ($_->[6] eq ($d eq '+' ? 'plus' : 'minus')) && $a2 >= $s - $o{window} && $a1 <= $e + $o{window} } @h;
        return () unless $best && $best->[4] >= 0.95 * length($x) && $best->[5] >= 80;
        my ($a1, $a2) = sort { $a <=> $b } ($best->[2], $best->[3]);
        # extend a hit that stops short of the exon ends
        if ($d eq '+') { $a1 -= ($best->[0] - 1); $a2 += (length($x) - $best->[1]) } else { $a2 += ($best->[0] - 1); $a1 -= (length($x) - $best->[1]) }
        push @placed, [$a1, $a2, $k] }
    # which placed exons are already annotated
    my @rows = map { [@$_] } @$ex; my @new;
    for my $pl (@placed) { my $have = grep { my $lo = $pl->[0] > $_->[0] ? $pl->[0] : $_->[0]; my $hi = $pl->[1] < $_->[1] ? $pl->[1] : $_->[1]; $hi - $lo + 1 > 0.5 * ($pl->[1] - $pl->[0] + 1) } @rows; push @new, $pl unless $have }
    return () unless @new;
    my @all = sort { $a->[0] <=> $b->[0] } (@rows, map { [$_->[0], $_->[1], undef] } @new);
    my $ni = intact_info(spliced($d, map { [$_->[0], $_->[1]] } @all), $g); return () unless $ni->{intact} && abs($ni->{len} - $want) <= 0.1 * $want;
    my @coding = $d eq '+' ? @all : reverse @all; my @edits;
    for my $k (0 .. $#coding) { my $r = $coding[$k]; my $name = "${g}_exon" . ($k + 1); if (defined $r->[2]) { push @edits, "rename $r->[2] $r->[0] $r->[1] $name" if $r->[2] ne $name } else { push @edits, "add $name $r->[0] $r->[1] $d" } $r->[3] = $name }
    for my $k (0 .. $#coding - 1) { my ($x1, $x2) = ($coding[$k], $coding[$k + 1]); my ($is, $ie) = $x1->[0] < $x2->[0] ? ($x1->[1] + 1, $x2->[0] - 1) : ($x2->[1] + 1, $x1->[0] - 1); my $iname = "${g}_intron" . ($k + 1);
        push @edits, "add $iname $is $ie +" unless grep { $_->[0] eq $iname } @rows_all }
    return (join(' ; ', @edits), $ni->{len}) }
sub span_union { my ($h, $s, $e) = @_; my ($lo, $hi) = ($s, $e); for my $x (@$h) { $lo = $x->{ss} if !defined $lo || $x->{ss} < $lo; $hi = $x->{se} if !defined $hi || $x->{se} > $hi } ($lo, $hi) }

# ---- table2asn validator messages --------------------------------------------------------
my @val; if ($o{val} && open my $vh, "<", $o{val}) { while (<$vh>) { chomp; next unless /^(ERROR|WARNING|REJECT|INFO|FATAL):\s*valid\s*\[([^\]]+)\]\s*(.*?)\s*FEATURE:\s*(\w+):\s*(.*?)\s*<[^:>]*:([^>]*)>/; my ($lev, $code, $msg, $type, $prod, $loc) = ($1, $2, $3, $4, $5, $6);
        my @n = ($loc =~ /(\d+)/g); next unless @n; my ($lo, $hi) = (sort { $a <=> $b } @n)[0, -1]; push @val, { lev => $lev, code => $code, msg => $msg, type => $type, prod => $prod, lo => $lo, hi => $hi } } close $vh }
sub val_for { my ($s, $e) = @_; join('; ', map { "$_->{lev} $_->{code}: $_->{msg}" } grep { $_->{lo} <= $e && $_->{hi} >= $s && ($_->{hi} - $_->{lo}) < 2 * ($e - $s + 1) + 100 } @val) }

# ---- verdicts -----------------------------------------------------------------------
print $out join("\t", qw(feature start end strand status flags evidence edit table2asn)), "\n";
my %count;
sub emit { my ($name, $s, $e, $d, $st, $flags, $ev, $edit) = @_; $count{$st}++; (my $ed = $edit // '') =~ s/\t/ /g; print $out join("\t", $name, $s, $e, $d, $st, $flags, $ev, $ed, val_for($s, $e)), "\n" }   # edits are space-separated inside the cell; turn spaces into tabs for verdant_edit.pl
sub rowname { my $f = shift; my @ex = @{ $f->{ex} }; @ex > 1 ? "$f->{gene} (" . scalar(@ex) . " exons)" : $f->{gene} }

for my $key (@order) { my $f = $feat{$key}; my ($g, $d) = ($f->{gene}, $f->{d}); my @ex = sort { $a->[0] <=> $b->[0] } @{ $f->{ex} }; my ($s, $e) = ($ex[0][0], $ex[-1][1]);
    if ($f->{pseudo}) { emit(rowname($f), $s, $e, $d, 'PSEUDOGENE', '', 'annotated as a pseudogene; not re-tested'); next }
    if ($f->{rna}) { my @fl; my $L = 0; $L += $_->[1] - $_->[0] + 1 for @ex; my $gl = guide_len($g);
        my $anyg = $gd{$g} ? grep { my $c = $_; grep { abs($L - $_) <= 3 } @{ $c->{alllens} } } values %{ $gd{$g} } : 0;   # a length some guide copy carries (the 57 nt second trnfM-CAU of the grasses)
        push @fl, sprintf("length %d vs guides %d", $L, $gl) if $gl && !$anyg && ($g =~ /^trn/ ? (abs($L - $gl) > 0.2 * $gl && ($L < 68 || $L > 95)) : abs($L - $gl) > 0.05 * $gl);
        my $gmax = $gd{$g} ? (sort { $b <=> $a } map { $_->{nex} } values %{ $gd{$g} })[0] : 1;
        my @dupe; for my $i (0 .. $#ex) { for my $j ($i + 1 .. $#ex) { my $lo = $ex[$i][0] > $ex[$j][0] ? $ex[$i][0] : $ex[$j][0]; my $hi = $ex[$i][1] < $ex[$j][1] ? $ex[$i][1] : $ex[$j][1]; push @dupe, [$i, $j] if $hi - $lo + 1 > 0 } }
        push @fl, scalar(@ex) . " pieces" if (@ex > 1 && $gmax == 1) || (@dupe && @ex > $gmax);   # more rows than any guide splits the gene into, or rows that overlap each other
        push @fl, "ambiguous bases" if ambiguous($s, $e);
        if (grep { /pieces$/ } @fl) { my @drop; if (@dupe) { for my $pr (@dupe) { my ($i, $j) = @$pr; my $odd = ($ex[$i][2] =~ /_exon\d+$/ && $ex[$j][2] !~ /_exon\d+$/) ? $ex[$j] : ($ex[$j][2] =~ /_exon\d+$/ && $ex[$i][2] !~ /_exon\d+$/) ? $ex[$i] : (($ex[$i][1] - $ex[$i][0]) < ($ex[$j][1] - $ex[$j][0]) ? $ex[$i] : $ex[$j]); push @drop, $odd unless grep { $_ == $odd } @drop } }
            else { my ($keep) = sort { abs(($a->[1] - $a->[0] + 1) - ($gl // 72)) <=> abs(($b->[1] - $b->[0] + 1) - ($gl // 72)) } @ex; @drop = grep { $_ != $keep } @ex }
            emit(rowname($f), $s, $e, $d, 'ARTIFACT', join('; ', @fl), (@dupe ? "rows of $g overlap each other: the unsuffixed or shorter row duplicates an exon call" : "the guides carry one $g of " . ($gl // '?') . " nt; the other row is a duplicate call"), join(' ; ', map { "delete $_->[2] $_->[0] $_->[1]" } @drop)) }
        elsif (@fl) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', join('; ', @fl), 'RNA gene: checked by length and composition only') } else { emit(rowname($f), $s, $e, $d, 'PASS', '', '') } next }
    next if $g =~ /^rps12/;   # trans-spliced; its parts are checked by the submission tool
    my $cds = spliced($d, map { [$_->[0], $_->[1]] } @ex); my $ii = intact_info($cds, $g); my $L = $ii->{len};
    my @fl; my ($E, $why) = expect($g); my $gl = guide_len($g); my $gnex = guide_nex($g); my $nc = carriers($g);
    my $want = $E // $gl // $L; my $tol = $E ? 0.05 : 0.10;
    push @fl, "length $L not a multiple of 3" if $ii->{mod3};
    push @fl, "start codon $ii->{start}" unless init_ok($ii->{start}, $g);
    push @fl, "no terminal stop" unless $ii->{term};
    push @fl, "$ii->{internal} internal stop" . ($ii->{internal} > 1 ? 's' : '') if $ii->{internal};
    my @okL = grep { defined } ($E, $gl, profile_lengths($g), ($gd{$g} ? map { $_->{len} } values %{ $gd{$g} } : ()));
    my $lenok = !($E || $gl) || grep { abs($L - $_) <= ($_ == ($E // -1) ? $tol * $_ : 6) } @okL; $lenok = 1 if $gl && abs($L - $gl) <= $tol * $gl;
    push @fl, sprintf("length %d vs expected %d (%s)", $L, $want, $E ? $why : "guides' mode") unless $lenok;
    push @fl, sprintf("%d exon%s vs guides' %d", scalar @ex, @ex > 1 ? 's' : '', $gnex) if $gnex && @ex != $gnex;
    my ($ov, $who) = overlap_other($key, $s, $e, $d); my $ovok = $ov > 0 && guides_overlap($g, $who, $ov); push @fl, "overlaps $who by $ov nt" if !$ovok && ($ov > 60 || ($ov > 30 && !$lenok));
    push @fl, "not a canonical gene and carried by $nc of $nguides guides" if !$CANON{$g} && $nguides && $nc < $nguides / 2;
    # a second, much shorter copy of a gene another copy carries in full: the IR-boundary remnant
    my $other_full = 0; for my $k (keys %feat) { next if $k eq $key; my $o2 = $feat{$k}; next unless $o2->{gene} eq $g; my $L2 = 0; $L2 += $_->[1] - $_->[0] + 1 for @{ $o2->{ex} }; $other_full = 1 if $L2 > 1.5 * $L && abs($L2 - $want) <= 0.1 * $want }
    push @fl, "partial second copy" if $other_full && $L < 0.5 * $want;
    # a start put out of frame by a homopolymer run: called whole it is UNRESOLVED; called anywhere else
    # (an internal ATG, an upstream TTG/ATA/AAT in the body's frame) it moves to the whole gene
    if (@ex == 1 && $gd{$g} && !$f->{pseudo}) {
        if (my ($rl, $rb, $rp) = hp_framed($cds, $g)) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', join('; ', (grep { !/terminal stop|internal stop/ } @fl), "start frameshifted by a run of $rl $rb"),   # the stops are the start's frame read past the run
                "the CDS runs from the ATG the guides start $g at to the stop: a run of $rl $rb at CDS position $rp puts that ATG out of frame with the body (the run's length differs between plastomes of the lineage; when it is a multiple of 3 the gene is intact). Whether slippage on the run restores the frame, translation begins at a downstream ATG, or $g is a pseudogene cannot be told from the sequence: confirm the run's length in the reads; for GenBank the CDS needs /pseudo or a /ribosomal_slippage join"); next }
        if (my ($ns, $ne, $NL, $rl, $rb, $na) = hp_shift($s, $e, $d, $g)) { emit(rowname($f), $s, $e, $d, 'ARTIFACT', join('; ', @fl, "start frameshifted by a run of $rl $rb"),
                "the ATG at " . ($d eq '+' ? $ns : $ne) . ", where $na guides' N-termini align, is out of frame with the body: a run of $rl $rb follows it. The call took " . ($NL > $L ? "an internal $ii->{start}, " . ($NL - $L) . " nt short" : "the $ii->{start} " . ($L - $NL) . " nt upstream, in the body's frame") . "; the whole gene is $NL nt from that ATG (not a multiple of 3, reported UNRESOLVED once moved)", "move\t$ex[0][2]\t$s\t$e\t$ns\t$ne"); next } }
    # a start the guides carry further upstream (post_filter_annotation.pl section 4c): most guides' N-termini align
    # more than 10 codons upstream of the call, and an in-frame ATG sits within 10 codons of where they align
    if (@ex == 1 && $gd{$g} && !$ii->{mod3}) { my @gl = map { $_->{len} } values %{ $gd{$g} }; my $longer = grep { $_ > $L + 30 } @gl; my @peps = guide_nterms($g);
        if (@gl >= 2 && 2 * $longer > @gl && @peps >= 2) { my ($gmax) = sort { $b <=> $a } @gl; my $best;
            my $avail = $d eq '+' ? $s - 1 : $len - $e; my $up = $gmax - $L + 60; $up = $avail if $up > $avail; $up -= $up % 3;
            if ($up >= 33) { my ($xs, $xe) = $d eq '+' ? ($s - $up, $e) : ($s, $e + $up); my $x = spliced($d, [$xs, $xe]); my $k0 = $up / 3;
                my $win = translate(substr($x, 0, $up + 90)); my @upa; for my $gp (@peps) { my ($sc, $at) = nterm_anchor($gp, $win); push @upa, $at if $sc >= 20 && $at < $k0 - 10 }
                if (2 * @upa > @peps) { my $med = (sort { $a <=> $b } @upa)[$#upa / 2];
                    for my $c (sort { abs($a - $med) <=> abs($b - $med) || $a <=> $b } grep { $_ >= 0 && $_ < $k0 } ($med - 10 .. $med + 10)) { my $cod = substr($x, 3 * $c, 3);
                        next unless $cod eq 'ATG' || ($cod eq 'ACG' && $EDITABLE{$g}); next if translate(substr($x, 3 * $c, 3 * ($k0 - $c))) =~ /\*/;
                        my ($ns, $ne) = $d eq '+' ? ($xs + 3 * $c, $e) : ($s, $xe - 3 * $c); $best = [$ns, $ne, $ne - $ns + 1, $cod, scalar @upa, scalar @peps]; last } } }
            if ($best && intact_info(spliced($d, [$best->[0], $best->[1]]), $g)->{intact}) {
                emit(rowname($f), $s, $e, $d, 'ARTIFACT', join('; ', @fl, "start downstream of the guides'"), "$best->[4] of $best->[5] guides' N-termini align more than 10 codons upstream of this $L nt call (guides carry " . join(',', sort { $a <=> $b } @gl) . " nt); the in-frame $best->[3] there, at " . ($d eq '+' ? $best->[0] : $best->[1]) . ", gives an intact $best->[2] nt CDS", "move\t$ex[0][2]\t$s\t$e\t$best->[0]\t$best->[1]"); next } } }
    unless (@fl) { emit(rowname($f), $s, $e, $d, 'PASS', '', ''); next }
    my $flags = join('; ', @fl);
    if (!init_ok($ii->{start}, $g) && !$ii->{mod3} && $ii->{term} && !$ii->{internal} && $gd{$g}) { my @same = grep { my $c = $_; grep { $_ eq $ii->{start} } @{ $c->{allstarts} || [] } } values %{ $gd{$g} };
        # a declared start that fewer than half of the guides share gives way to an in-frame ATG/GTG at a length
        # the lineage or the guides carry (one Lilium guide declares TTG for cemA; the other seven start at ATG)
        if (@same && 2 * @same < scalar(keys %{ $gd{$g} }) && @ex == 1) { my %okl = map { $_ => 1 } ((defined $E ? $E : ()), map { @{ $_->{alllens} || [] } } values %{ $gd{$g} }); my @cand;
            my ($xs, $xe) = $d eq '+' ? ($s - 90, $e) : ($s, $e + 90); my @anch;   # window index 30 is the current start
            if ($xs >= 1 && $xe <= $len) { my $win = translate(substr(spliced($d, [$xs, $xe]), 0, 183)); for my $gp (guide_nterms($g)) { my ($sc, $at) = nterm_anchor($gp, $win); push @anch, $at if $sc >= 20 } }
            my $anchor = @anch ? (sort { $a <=> $b } @anch)[$#anch / 2] : undef;
            for my $k (-30 .. 30) { next unless $k; my ($ns, $ne) = $d eq '+' ? ($s + 3 * $k, $e) : ($s, $e - 3 * $k); next if $ns < 1 || $ne > $len || $ne - $ns + 1 < 30; my $nl = $ne - $ns + 1;
                my $near = defined $anchor ? abs($k + 30 - $anchor) : 99; next unless $near <= 3 || $okl{$nl};
                my $x = spliced($d, [$ns, $ne]); my $c = substr($x, 0, 3); next unless $c eq 'ATG' || $c eq 'GTG'; my $p = translate($x); next unless index($p, '*') == length($p) - 1;
                push @cand, [ $near <= 3 ? $near : 9, (defined $E && $nl == $E) ? 0 : 1, $c eq 'ATG' ? 0 : 1, abs $k, $ns, $ne, $nl, $c ] }
            if (@cand) { my ($b) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @cand;
                emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "only " . scalar(@same) . " of " . scalar(keys %{ $gd{$g} }) . " guides carrying $g begin it with $ii->{start}; the in-frame $b->[7] at " . ($d eq '+' ? $b->[4] : $b->[5]) . " gives $b->[6] nt, intact, " . ($b->[0] <= 3 ? "where " . scalar(@anch) . " guides' N-termini align" : $b->[1] ? "a guide copy's length" : "the lineage expectation"), "move\t$ex[0][2]\t$s\t$e\t$b->[4]\t$b->[5]"); next }
            # no normal start in frame: the minority's declared start may be right, but it is not taken on trust
            my $oof = 0; for my $q (-45 .. 45) { next unless $q % 3; my $p = $d eq '+' ? $s + $q : $e - $q; next if $p < 3 || $p + 2 > $len; $oof = 1 if spliced($d, $d eq '+' ? [$p, $p + 2] : [$p - 2, $p]) eq 'ATG' }
            emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "only " . scalar(@same) . " of " . scalar(keys %{ $gd{$g} }) . " guides carrying $g begin it with $ii->{start} (declared in their records); no in-frame ATG/GTG where the others' N-termini align" . ($oof ? "; an ATG lies out of frame within 15 codons, as a 1-nt indel in a run near the start would leave it: check the assembly there" : "")); next }
        if (@same) { my $lenfit = grep { my $c = $_; grep { abs($L - $_) <= 6 } @{ $c->{alllens} } } @same;
            emit(rowname($f), $s, $e, $d, ($lenfit ? 'PASS' : 'UNRESOLVED'), $flags, scalar(@same) . " of $nguides guides begin their own $g with $ii->{start} (a declared start, transl_except in those records)" . ($lenfit ? "; the length matches a guide copy" : "; but no guide copy has this length")); next } }
    if ($other_full && $L < 0.5 * $want) { emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "a full-length $g (" . ($want) . " nt expected) is annotated elsewhere; this is the remnant copy at the IR boundary", "delete\t$ex[0][2]\t$s\t$e"); next }
    if (my $n = ambiguous($s - 30, $e + 30)) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$n ambiguous base" . ($n > 1 ? 's' : '') . " in or beside the feature"); next }
    if (!$gd{$g}) { emit(rowname($f), $s, $e, $d, ($CANON{$g} ? 'UNRESOLVED' : 'ARTIFACT'), $flags, ($CANON{$g} ? "no guide carries $g; cannot re-derive" : "no guide carries $g and it is not a canonical plastid gene"), ($CANON{$g} ? undef : "delete\t$ex[0][2]\t$s\t$e")); next }
    my ($h, $plen, $gname, $glen) = hsps($g, $s, $e, $d);
    unless ($h && @$h) { emit(rowname($f), $s, $e, $d, 'REAL', $flags, "guide protein ($gname) finds nothing near $s-$e: no homologous sequence", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene"); next }
    my $cov = coverage($h, $plen); my %frames = map { $_->{frame} => 1 } @$h; my $nfr = scalar keys %frames;
    my ($us, $ue) = span_union($h, $s, $e); my $bir = near_ir($us, $ue, 0);
    my $evb = sprintf("guide %s: %d HSP%s, %.0f%% of %d aa, %d frame%s", $gname, scalar @$h, @$h > 1 ? 's' : '', 100 * $cov, $plen, $nfr, $nfr > 1 ? 's' : '');
    my $multi = ($gnex && $gnex > 1) || @ex > 1;
    if (!$multi) {
        if ($ii->{term} && !$ii->{internal} && !$ii->{mod3} && !init_ok($ii->{start}, $g) && $gd{$g} && grep { my $c = $_; grep { abs($L - $_) <= 6 } @{ $c->{alllens} } } values %{ $gd{$g} }) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; the CDS is open, complete and of a length the guides carry, but begins $ii->{start}: a declared start (transl_except) in lineage records is possible; not called a pseudogene on the start alone"); next }
        # HSPs in the frame of the call itself, and whether the others sit N-terminal of them
        my @main = grep { $d eq '+' ? (($_->{ss} - $s) % 3 == 0) : (($e - $_->{se}) % 3 == 0) } @$h; my $mcov = coverage(\@main, $plen);
        my ($mqs) = sort { $a <=> $b } map { $_->{qs} } @main; my ($mqe) = sort { $b <=> $a } map { $_->{qe} } @main;
        my $nterm_only = @main && !grep { $_->{qe} > ($mqs // 0) + 5 } grep { my $x = $_; !grep { $_ == $x } @main } @$h;
        if ($ii->{intact} && $mcov >= 0.7 && ($mqe // 0) >= $plen - 3 && $nterm_only && $nfr > 1) { emit(rowname($f), $s, $e, $d, 'PASS', $flags, "$evb; the call is intact and carries the protein from residue $mqs to the end in frame; the guide's first " . ($mqs - 1) . " residues are not conserved in frame, so the guide's own start is the over-extended one"); next }
        if ($ii->{intact} && $ov > 60 && !$ovok) { my $gotL = 0; my ($ns, $ne) = ($s, $e);
            PL: for my $PL (profile_lengths($g)) { next if $PL >= $L - 30; my $p = $d eq '+' ? $s : $e;
                for my $k (1 .. 200) { if ($d eq '+') { $p += 3; last if $p > $e - 90; my $c = substr($seq, $p - 1, 3); last if ($CODON{$c} // '') eq '*'; next unless $INIT{$c}; my $NL = $e - $p + 1; if (abs($NL - $PL) <= 6) { ($ns, $gotL) = ($p, $NL); last PL } }
                    else { $p -= 3; last if $p < $s + 90; my $c = rc(substr($seq, $p - 3, 3)); last if ($CODON{$c} // '') eq '*'; next unless $INIT{$c}; my $NL = $p - $s + 1; if (abs($NL - $PL) <= 6) { ($ne, $gotL) = ($p, $NL); last PL } } } }
            if ($gotL) { my ($ov2) = overlap_other($key, $ns, $ne, $d); if ($ov2 <= 30) { emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "$evb; the call overlaps $who by $ov nt; the in-frame ATG at " . ($d eq '+' ? $ns : $ne) . " gives $gotL nt, the lineage profile's length, with no overlap: over-extended start", "move\t$ex[0][2]\t$s\t$e\t$ns\t$ne"); next } } }
        if ($nfr == 1 && $cov >= 0.9) { my ($hs, $he) = ($h->[0]{ss}, $h->[0]{se}); for my $x (@$h) { $hs = $x->{ss} if $x->{ss} < $hs; $he = $x->{se} if $x->{se} > $he }
            my ($ns, $ne) = orf_from_hit($hs, $he, $d, $want, $g);
            if (defined $ns) { my $ni = intact_info(spliced($d, [$ns, $ne]), $g); my $NL = $ne - $ns + 1;
                if ($ni->{intact} && abs($NL - $want) <= 0.1 * $want) {
                    if ($ns == $s && $ne == $e) { emit(rowname($f), $s, $e, $d, 'PASS', $flags, "$evb; the call is the intact ORF the guide protein maps to; the deviation is the gene's own") }
                    else { emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "$evb; an intact $NL nt ORF at $ns-$ne $d (start " . $ni->{start} . ") matches the expectation; the annotation is " . ($L > $NL ? "over-extended" : "short"), "move\t$ex[0][2]\t$s\t$e\t$ns\t$ne") } next }
                if ($ni->{internal} && $ii->{intact} && $ns <= $s && $ne >= $e) { emit(rowname($f), $s, $e, $d, 'PASS', $flags, "$evb; the frame the guide protein maps to carries a stop upstream of the call's start, so the guide's N-terminal extension is not conserved and the call, which is intact, starts after it"); next }
                if ($ni->{internal} && !$bir) { emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb; the frame the protein maps to ($ns-$ne) stops early: $ni->{internal} internal stop" . ($ni->{internal} > 1 ? 's' : ''), "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene ; move\t${g}_pseudogene\t$ex[0][0]\t$ex[0][1]\t$us\t$ue"); next } }
            if ($bir) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; an IR boundary at $bir lies inside the homologous span $us-$ue"); next }
            if (my ($op, $oc) = offframe_start($h, $d)) { emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb; the only canonical start near the aligned N-terminus ($oc at $op) is out of frame with the body: a frameshift lies between it and the protein", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene ; move\t${g}_pseudogene\t$ex[0][0]\t$ex[0][1]\t$us\t$ue"); next }
            emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb over $us-$ue but no intact ORF of the expected length can be built there", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene ; move\t${g}_pseudogene\t$ex[0][0]\t$ex[0][1]\t$us\t$ue"); next }
        if ($bir) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; the homologous span $us-$ue crosses the IR boundary at $bir"); next }
        if ($nfr > 1 && $cov >= 0.5) { emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb: the protein aligns in $nfr frames over $us-$ue (frameshift)", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene ; move\t${g}_pseudogene\t$ex[0][0]\t$ex[0][1]\t$us\t$ue"); next }
        emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb: only part of the protein is present over $us-$ue", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene ; move\t${g}_pseudogene\t$ex[0][0]\t$ex[0][1]\t$us\t$ue"); next }
    # a start-only failure: an in-frame ATG/GTG (ACG for the editable genes) within 10 codons of the current start
    if (!$ii->{mod3} && $ii->{term} && !$ii->{internal} && !init_ok($ii->{start}, $g)) { my $first = $d eq '+' ? $ex[0] : $ex[-1]; my ($fs, $fe) = @$first[0, 1]; my @cand;
        for my $k (-10 .. 10) { next unless $k; my ($ns, $ne) = ($fs, $fe); if ($d eq '+') { $ns = $fs + 3 * $k; next if $ns < 1 || $ns >= $fe - 3; my $c = substr($seq, $ns - 1, 3); next unless init_ok($c, $g); push @cand, [$ns, $ne, $k, $c] }
            else { $ne = $fe - 3 * $k; next if $ne > $len || $ne <= $fs + 3; my $c = rc(substr($seq, $ne - 3, 3)); next unless init_ok($c, $g); push @cand, [$ns, $ne, $k, $c] } }
        my ($bestc) = grep { my @nx = map { [@$_] } @ex; my $i2 = $d eq '+' ? 0 : $#ex; ($nx[$i2][0], $nx[$i2][1]) = ($_->[0], $_->[1]); my $ni = intact_info(spliced($d, map { [$_->[0], $_->[1]] } @nx), $g); $ni->{intact} && abs($ni->{len} - $want) <= 0.1 * $want } sort { abs($a->[2]) <=> abs($b->[2]) } @cand;
        if ($bestc) { emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "$evb; the spliced CDS is intact once the start moves " . abs($bestc->[2]) . " codon" . (abs($bestc->[2]) > 1 ? 's' : '') . ($bestc->[2] > 0 ? ' downstream' : ' upstream') . " to $bestc->[3]", "move\t$first->[2]\t$fs\t$fe\t$bestc->[0]\t$bestc->[1]"); next }
        my $lenfit = $gd{$g} && grep { my $c = $_; grep { abs($L - $_) <= 6 } @{ $c->{alllens} } } values %{ $gd{$g} };
        if ($lenfit) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; the CDS is open, complete and of a length the guides carry, but begins $ii->{start}: lineage records declare such starts with a transl_except (the LSC rpl23 copy of the grasses); not called a pseudogene on the start alone"); next }
        if (my ($op, $oc) = offframe_start($h, $d)) { emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb; no in-frame canonical start within 10 codons, but $oc at $op is out of frame with the body: a frameshift lies between it and the protein", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene"); next }
        emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; the CDS is open and complete but begins $ii->{start} with no canonical start within 10 codons in any frame; GenBank records of the lineage may carry a transl_except here"); next }
    # spliced genes: try every junction shift within +-3 nt for an intact CDS of the expected length
    if (@ex > 1 && !$ii->{intact}) { my @best; my @j = (0 .. $#ex - 1); my @sh = (-3 .. 3);
        my $try; $try = sub { my ($i, @cur) = @_; if ($i > $#j) { my @nx = map { [@$_] } @ex; for my $k (0 .. $#j) { $nx[$k][1] += $cur[$k]; $nx[$k + 1][0] += $cur[$k] } my $ni = intact_info(spliced($d, map { [$_->[0], $_->[1]] } @nx), $g); my $NL = $ni->{len};
                push @best, [\@nx, $NL, $cur[0] // 0] if $ni->{intact} && abs($NL - $want) <= 0.1 * $want; return } for my $v (@sh) { $try->($i + 1, @cur, $v) } };
        $try->(0) if @j <= 2; if (@best) { my ($b) = sort { abs($a->[1] - $want) <=> abs($b->[1] - $want) } @best; my @edits; for my $k (0 .. $#ex) { my ($os, $oe, $on) = @{ $ex[$k] }; my ($ns, $ne) = @{ $b->[0][$k] }; push @edits, "move\t$on\t$os\t$oe\t$ns\t$ne" if $ns != $os || $ne != $oe }
            emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "$evb; shifting an exon boundary within 3 nt gives an intact $b->[1] nt CDS", join(" ; ", @edits)); next } }
    if ($gnex && @ex < $gnex && $cov < 0.9) { my @lost = grep { my $q = $_; !grep { $_->{qs} <= $q && $_->{qe} >= $q } @$h } (1 .. $plen); emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb: " . scalar(@lost) . " of $plen aa have no homologous sequence near $s-$e; an exon is gone", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene"); next }
    if ($gnex && @ex < $gnex) { my ($ed, $NL) = place_missing_exons($g, $gname, $d, $s, $e, [map { [$_->[0], $_->[1], $_->[2]] } @ex], $want);
        if ($ed) { emit(rowname($f), $s, $e, $d, 'ARTIFACT', $flags, "$evb; the guide's missing exon maps by blastn beside the call and the completed splice is an intact $NL nt CDS", $ed); next }
        emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; the protein is present but an exon was not placed (expected $gnex exons)"); next }
    if ($bir) { emit(rowname($f), $s, $e, $d, 'UNRESOLVED', $flags, "$evb; the IR boundary at $bir lies in the span $us-$ue"); next }
    if ($ii->{intact}) { emit(rowname($f), $s, $e, $d, 'PASS', $flags, "$evb; the spliced CDS is intact; the deviation is the gene's own"); next }
    emit(rowname($f), $s, $e, $d, 'REAL', $flags, "$evb; no boundary shift within 3 nt makes the spliced CDS intact", "rename\t$ex[0][2]\t$ex[0][0]\t$ex[0][1]\t${g}_pseudogene") }

# ---- canonical genes the guides carry but the annotation lacks ------------------------
for my $g (sort keys %gd) { next if $present{$g} || !$CANON{$g} || $g =~ /^rps12/; my $nc = carriers($g); next unless $nguides && $nc >= $nguides / 2;
    my ($h, $plen, $gname, $glen) = hsps($g); my $want = (expect($g))[0] // $glen;
    unless ($h && @$h) { emit($g, 0, 0, '', 'REAL', "absent; carried by $nc of $nguides guides", "guide protein ($gname) finds no homologous sequence anywhere in the genome: gene lost"); next }
    my ($best) = sort { $b->{alen} <=> $a->{alen} } @$h; my $d = $best->{d}; my @near = grep { $_->{d} eq $d && $_->{se} >= $best->{ss} - $o{window} && $_->{ss} <= $best->{se} + $o{window} } @$h;
    my $cov = coverage(\@near, $plen); my %frames = map { $_->{frame} => 1 } @near; my $nfr = keys %frames; my ($us, $ue) = span_union(\@near); my $bir = near_ir($us, $ue, 30);
    my $evb = sprintf("absent; carried by %d of %d guides; guide %s: %d HSP%s, %.0f%% of %d aa, %d frame%s over %d-%d %s", $nc, $nguides, $gname, scalar @near, @near > 1 ? 's' : '', 100 * $cov, $plen, $nfr, $nfr > 1 ? 's' : '', $us, $ue, $d);
    if ($nfr == 1 && $cov >= 0.9 && !(guide_nex($g) && guide_nex($g) > 1)) { my ($ns, $ne) = orf_from_hit($best->{ss}, $best->{se}, $d, $want, $g); if (defined $ns) { my $ni = intact_info(spliced($d, [$ns, $ne]), $g); if ($ni->{intact} && abs(($ne - $ns + 1) - $want) <= 0.1 * $want) { emit($g, $us, $ue, $d, 'ARTIFACT', 'absent', "$evb: an intact ORF of " . ($ne - $ns + 1) . " nt exists at $ns-$ne", "add\t$g\t$ns\t$ne\t$d"); next } } }
    if ($bir) { emit($g, $us, $ue, $d, 'UNRESOLVED', 'absent', "$evb; the span crosses the IR boundary at $bir"); next }
    if ($cov < 0.3) { emit($g, $us, $ue, $d, 'REAL', 'absent', "$evb: only a trace remains"); next }
    emit($g, $us, $ue, $d, 'REAL', 'absent', "$evb: " . ($nfr > 1 ? 'frameshifted' : 'interrupted') . " remnant", "add\t${g}_pseudogene\t$us\t$ue\t$d") }

# ---- tRNAs the guides carry but the annotation lacks ----------------------------------
for my $g (sort grep { /^trn/ } keys %gd) { next if $present{$g}; my $nc = carriers($g); next unless $nguides && $nc >= $nguides / 2; emit($g, 0, 0, '', 'UNRESOLVED', "absent; carried by $nc of $nguides guides", "tRNA not annotated; run ARAGORN on the region to decide") }

printf STDERR "%s: %s\n", $ann, join(', ', map { "$_ $count{$_}" } sort keys %count);
