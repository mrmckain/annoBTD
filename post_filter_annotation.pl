#!/usr/bin/perl
# Post-annotation filter, run on the final annotation before the checkers.
#
# 1. Duplicate RNA calls. With many guides, one tRNA locus is transferred under
#    several names ("trnK" from an older record, "trnK-UUU" from a newer one,
#    trnI-CAU vs trnM-CAU where GenBank itself is confused) and each name gets its
#    own call a few nucleotides apart. tRNA/rRNA features overlapping >= 50% of the
#    shorter one on the same strand collapse to one: the name carrying an
#    anticodon wins, then the name more guides use, then the longer feature.
#
# 2. Presence prior. A CDS gene annotated in fewer than --min-presence of the
#    target lineage's species (family with >= 10 species, else order) is dropped
#    WHEN the guides are split on it too (no more than half of them annotate it):
#    ycf15, ycf68, orfNN and the like that a minority of labs annotate and the
#    curated truth omits. Share = species with the gene (gene_expect, one vote per
#    species) / species in the lineage (lineage_species_counts.tsv). The guide
#    condition is what separates those from a real gene a lineage under-annotates:
#    Asparagaceae infA sits at 0.47 of species, yet 13 of 16 guides carry it, and
#    the curated genomes all have it.
#
# 3. tRNA window refinement (--aragorn BIN). Guides carry their submitters' tRNA
#    window conventions (modern GenBank records sit 4 nt from where the gene
#    folds). ARAGORN finds the gene from the sequence; each tRNA call is snapped
#    to the overlapping ARAGORN gene (same strand, >= 50% of the call covered,
#    same intron status; ARAGORN's identity is ignored - it misnames the spliced
#    plastid tRNAs) and the per-identity house convention (--trna-convention,
#    truth-minus-ARAGORN offsets) is applied. Spliced calls keep their exon
#    lengths anchored to the corrected gene ends.
#
# Removing a feature merges its flanking intergenic rows (A~B, B~C -> A~C).
#
# USAGE: post_filter_annotation.pl <plastome.fsa> <annotation.txt> [--guides-dir DIR]
#        [--profile gene_expect.tsv --lineage-counts lineage_species_counts.tsv
#         --family F --order O --min-presence 0.5] [--no-rna-dedupe] [--log FILE]
#        [--aragorn BIN --trna-convention trna_window_convention.tsv]
#        [--species-units species_units.tsv --species-genes species_lengths.tsv --presence-mode family|local]
#    presence-mode local: the guide-count condition is replaced by the share of the
#    guides' SPECIES (their one-vote-per-group consensus gene sets) that carry the gene,
#    so one guide's odd record no longer decides; falls back to the guide count when no
#    guide resolves to a species.
use strict; use warnings; use FindBin;
my ($fsa, $ann, @rest) = @ARGV; die "usage: $0 <plastome.fsa> <annotation.txt> [options]\n" unless $ann;
my %o = ('min-presence' => 0.5);
while (@rest) { my $a = shift @rest; if ($a =~ /^--(no-rna-dedupe)$/) { $o{$1} = 1 } elsif ($a =~ /^--(.+)$/) { $o{$1} = shift @rest } }
my ($len, $seq) = (0, ''); if (open my $f, "<", $fsa) { while (<$f>) { next if /^>/; s/\s//g; $seq .= uc $_ } close $f; $len = length $seq }
# How well a window's ends pair as a tRNA acceptor stem (7 bp, G:U allowed).
# Annotation conventions shift tRNA windows by a few nucleotides (seven modern Zea
# guides put trnH-GUG 4 nt downstream of the classic record; the classic window
# pairs 7/7, the shifted one 4/7), and this is the arbiter a vote cannot be.
my %pair = map { $_ => 1 } qw(AT TA GC CG GT TG);
sub stem_score { my ($s0, $e0, $dir) = @_; return 0 unless $seq && $s0 >= 1 && $e0 <= $len; my $t = substr($seq, $s0 - 1, $e0 - $s0 + 1);
    if ($dir eq '-') { $t = reverse $t; $t =~ tr/ACGT/TGCA/ } my $L = length $t; my $n = 0; for my $i (0 .. 6) { $n++ if $pair{ substr($t, $i, 1) . substr($t, $L - 1 - $i, 1) } } $n }

my (@rows, @feat);
open my $h, "<", $ann or die "$ann: $!";
while (<$h>) { chomp; next unless /\S/; my @f = split /\t/; push @rows, \@f;
    next if $f[0] =~ /~/ || $f[0] =~ /^(LSC|SSC|IRA|IRB|FULL)$/; push @feat, $#rows }
my @ir = map { [$_->[1], $_->[2]] } grep { $_->[0] =~ /^IR[AB]$/ } @rows;
sub in_ir { my ($s, $e) = @_; for my $r (@ir) { return 1 if $s >= $r->[0] && $e <= $r->[1] } 0 }
close $h;
my $log = $o{log} ? do { open my $l, ">", $o{log} or die; $l } : \*STDERR;
my %drop;

# --- 1. RNA de-duplication ---------------------------------------------------
unless ($o{'no-rna-dedupe'}) {
    my %support;
    if ($o{'guides-dir'} && opendir(my $d, $o{'guides-dir'})) {
        for my $g (grep { !/^\./ } readdir $d) { my %seen;
            open my $gh, "<", "$o{'guides-dir'}/$g" or next;
            while (<$gh>) { chomp; my @f = split /\t/; next unless ($f[4] // '') =~ /^(tRNA|rRNA)$/; (my $n = $f[0]) =~ s/_exon\d+$//; $seen{$n} = 1 }
            close $gh; $support{$_}++ for keys %seen }
        closedir $d;
    }
    my @rna = grep { $rows[$_][0] =~ /^(trn|rrn)/ } @feat;
    # An exon of a spliced tRNA holds only half the acceptor stem, so score the
    # spliced product: this exon joined to its nearest sibling exon of the same
    # full name (IR copies of a name are told apart by distance).
    my %stem_cache;
    sub row_stem { my $k = shift; return $stem_cache{$k} if exists $stem_cache{$k}; my $r = $rows[$k];
        my $v;
        if ($r->[0] =~ /^(.+)_exon(\d+)$/) { my ($base, $ex) = ($1, $2);
            my @sib = sort { abs($rows[$a][1] - $r->[1]) <=> abs($rows[$b][1] - $r->[1]) } grep { $_ != $k && $rows[$_][0] =~ /^\Q$base\E_exon\d+$/ && $rows[$_][3] eq $r->[3] && abs($rows[$_][1] - $r->[1]) < 5000 } @rna;
            if (@sib) { my @parts = sort { $a->[1] <=> $b->[1] } ($r, $rows[$sib[0]]);
                my $t = join('', map { substr($seq, $_->[1] - 1, $_->[2] - $_->[1] + 1) } @parts);
                if ($r->[3] eq '-') { $t = reverse $t; $t =~ tr/ACGT/TGCA/ } my $L = length $t; my $n = 0;
                for my $i (0 .. 6) { $n++ if $pair{ substr($t, $i, 1) . substr($t, $L - 1 - $i, 1) } } $v = $n } else { $v = 0 } }
        else { $v = stem_score($r->[1], $r->[2], $r->[3]) }
        $stem_cache{$k} = $v }
    for my $i (0 .. $#rna) { for my $j ($i + 1 .. $#rna) {
        my ($a, $b) = ($rows[$rna[$i]], $rows[$rna[$j]]); next if $drop{$rna[$i]} || $drop{$rna[$j]};
        next unless $a->[3] eq $b->[3];
        my $lo = $a->[1] > $b->[1] ? $a->[1] : $b->[1]; my $hi = $a->[2] < $b->[2] ? $a->[2] : $b->[2];
        my $ov = $hi - $lo + 1; next if $ov <= 0;
        my ($la, $lb) = ($a->[2] - $a->[1] + 1, $b->[2] - $b->[1] + 1); next if $ov < 0.5 * ($la < $lb ? $la : $lb);
        (my $na = $a->[0]) =~ s/_exon\d+$//; (my $nb = $b->[0]) =~ s/_exon\d+$//;
        # exons of one spliced tRNA never overlap each other; same-name exon pairs are IR copies, skip
        next if $na eq $nb && $a->[0] ne $b->[0];
        # an unspliced call overlapping an exon of a spliced tRNA is a guide's intron-less
        # annotation of that gene: the exon wins outright (a half tRNA cannot win a stem test)
        my ($xa, $xb) = ($a->[0] =~ /_exon\d+$/ ? 1 : 0, $b->[0] =~ /_exon\d+$/ ? 1 : 0);
        if ($xa != $xb) { my ($w, $l) = $xa ? ($rna[$i], $rna[$j]) : ($rna[$j], $rna[$i]); $drop{$l} = 1;
            printf $log "RNA_DUP\tdrop %s %d-%d %s\tkeep %s %d-%d (unspliced call inside a spliced tRNA)\n", $rows[$l][0], $rows[$l][1], $rows[$l][2], $rows[$l][3], $rows[$w][0], $rows[$w][1], $rows[$w][2]; next }
        # coordinates: acceptor-stem pairing first, then guide support, then length;
        # name: decided separately below. An exon feature can only belong to one of
        # the six intron-containing plastid tRNAs; a guide's "trnE-UUC_exon1" over a
        # trnI-GAU exon is a mislabel and loses outright.
        my %spliced = map { $_ => 1 } qw(trnI-GAU trnA-UGC trnK-UUU trnG-UCC trnL-UAA trnV-UAC trnI trnA trnK trnG trnL trnV);
        my $bad_a = ($a->[0] =~ /_exon\d+$/ && !$spliced{$na}) ? 1 : 0; my $bad_b = ($b->[0] =~ /_exon\d+$/ && !$spliced{$nb}) ? 1 : 0;
        my ($sta, $stb) = (row_stem($rna[$i]), row_stem($rna[$j]));
        my $keep_a = ($bad_b <=> $bad_a || $sta <=> $stb || ($support{$na} // 0) <=> ($support{$nb} // 0) || $la <=> $lb || $b->[0] cmp $a->[0]) >= 0;
        my @sa = ($sta, $support{$na} // 0); my @sb = ($stb, $support{$nb} // 0);
        # trnI-CAU / trnM-CAU: GenBank names the two CAU tRNAs inconsistently and a
        # vote just repeats the confusion. Position decides: the CAU tRNA inside the
        # inverted repeat is trnI-CAU (lysidine-modified isoleucine), the one in the
        # LSC is the elongator trnM-CAU.
        if ({ map { $_ => 1 } $na, $nb }->{'trnI-CAU'} && { map { $_ => 1 } $na, $nb }->{'trnM-CAU'} && @ir) {
            my $want = in_ir($lo, $hi) ? 'trnI-CAU' : 'trnM-CAU'; $keep_a = ($na eq $want) }
        my $loser = $keep_a ? $rna[$j] : $rna[$i]; my $winner = $keep_a ? $rna[$i] : $rna[$j];
        $drop{$loser} = 1;
        # the NAME is decided separately from the coordinates: the name more guides
        # use, then the one carrying an anticodon (a lone guide's trnT-CGU, a tRNA no
        # other record has, must not rename everyone else's trnG-UCC)
        unless ({ map { $_ => 1 } $na, $nb }->{'trnI-CAU'} && { map { $_ => 1 } $na, $nb }->{'trnM-CAU'}) {
            my ($la_, $lb_) = map { /^trn([A-Za-z]{1,2})/ ? $1 : $_ } ($na, $nb);
            my ($best) = ($bad_a != $bad_b) ? ($bad_a ? $nb : $na) : ($la_ eq $lb_)
                ? (sort { ($b =~ /-[ACGU]{3}$/ ? 1 : 0) <=> ($a =~ /-[ACGU]{3}$/ ? 1 : 0) || ($support{$b} // 0) <=> ($support{$a} // 0) || $a cmp $b } ($na, $nb))   # same tRNA: keep the anticodon
                : (sort { ($support{$b} // 0) <=> ($support{$a} // 0) || ($b =~ /-[ACGU]{3}$/ ? 1 : 0) <=> ($a =~ /-[ACGU]{3}$/ ? 1 : 0) || $a cmp $b } ($na, $nb));
            my $ex = $rows[$winner][0] =~ /(_exon\d+)$/ ? $1 : ''; $rows[$winner][0] = $best . $ex }
        printf $log "RNA_DUP\tdrop %s %d-%d %s\tkeep %s %d-%d (stem %d/%d, guides %d/%d)\n", $rows[$loser][0], $rows[$loser][1], $rows[$loser][2], $rows[$loser][3],
            $rows[$winner][0], $rows[$winner][1], $rows[$winner][2], ($keep_a ? ($sa[0], $sb[0], $sa[1], $sb[1]) : ($sb[0], $sa[0], $sb[1], $sa[1]));
    } }
}

# --- 1b. tRNA identities no plastid carries -------------------------------------
# A tRNA name whose letter+anticodon appears in none of the curated genomes
# (annoBTD/plastid_trna_names.txt) and that only one guide supports is a guide's
# mistake (a "trnT-CGU" from one mislabeled record), not a discovery.
{
    my $list = (grep { -s $_ } map { "$_/plastid_trna_names.txt" } grep { defined } $ENV{ANNOBTD_DIR}, $FindBin::Bin)[0];
    if ($list && !$o{'no-rna-dedupe'}) { my %ok; open my $fh, "<", $list or die; while (<$fh>) { chomp; $ok{$_} = 1 if /\S/ } close $fh;
        my %support; if ($o{'guides-dir'} && opendir(my $d, $o{'guides-dir'})) {
            for my $g (grep { !/^\./ } readdir $d) { my %seen; open my $gh, "<", "$o{'guides-dir'}/$g" or next;
                while (<$gh>) { chomp; my @f = split /\t/; next unless ($f[4] // '') eq 'tRNA'; (my $n = $f[0]) =~ s/_exon\d+$//; $seen{$n} = 1 } close $gh; $support{$_}++ for keys %seen } closedir $d }
        for my $i (@feat) { next if $drop{$i}; my $n = $rows[$i][0]; next unless $n =~ /^trn/; (my $b = $n) =~ s/_exon\d+$//; next unless $b =~ /-[ACGU]{3}$/; next if $ok{$b}; next if ($support{$b} // 0) > 1;
            $drop{$i} = 1; printf $log "RNA_IDENTITY\tdrop %s %d-%d %s\tno curated plastid has %s; %d guide(s)\n", $n, $rows[$i][1], $rows[$i][2], $rows[$i][3], $b, $support{$b} // 0 }
    }
}

# --- 2. presence prior --------------------------------------------------------
if ($o{profile} && $o{'lineage-counts'} && -s $o{profile} && -s $o{'lineage-counts'}) {
    my %nsp; open my $c, "<", $o{'lineage-counts'} or die; while (<$c>) { chomp; my ($t, $l, $n) = split /\t/; $nsp{"$t:$l"} = $n } close $c;
    my %with; open my $p, "<", $o{profile} or die; <$p>;
    while (<$p>) { chomp; my ($g, $lt, $lin, $n) = split /\t/; next if $g =~ /_exon\d+$/; $with{$g}{"$lt:$lin"} = $n } close $p;
    # local prior: gene sets of the guides' species consensus forms
    my (%spgene, %nsp_local); if (($o{'presence-mode'} // 'family') eq 'local' && $o{'species-units'} && $o{'species-genes'} && $o{'guides-dir'}) {
        my (%unit_of, %rep_of); if (open my $uh, "<", $o{'species-units'}) { while (<$uh>) { chomp; my @f = split /\t/; next if $f[0] eq 'accession'; my $sp = $f[1]; for my $a ($f[0], ($f[8] ne '-' ? split(/,/, $f[8]) : ())) { $unit_of{$a} = $sp; (my $b = $a) =~ s/\.\d+$//; $unit_of{$b} = $sp }   # versioned and base accessions
            $rep_of{$sp} = $f[0] if $f[7] eq '1' } close $uh }
        my %want_rep; if (opendir(my $d, $o{'guides-dir'})) { for my $g (grep { !/^\./ } readdir $d) { my ($acc) = $g =~ /([A-Z]{1,2}_?\d+(?:\.\d+)?)$/; next unless $acc; my $sp = $unit_of{$acc}; next unless $sp && $rep_of{$sp}; $want_rep{ $rep_of{$sp} } = $sp } closedir $d }
        if (%want_rep && open my $lh, "<", $o{'species-genes'}) { my %have; while (<$lh>) { my ($a, $name) = split /\t/; next unless $want_rep{$a}; next if $name =~ /_exon\d+$/; $have{$a}{$name} = 1 } close $lh;
            $nsp_local{n} = scalar keys %have; for my $a (keys %have) { $spgene{$_}++ for keys %{ $have{$a} } } }
    }
    # guide support per CDS gene: how many guide annotations carry it
    my (%gsup, $ng); if ($o{'guides-dir'} && opendir(my $d, $o{'guides-dir'})) {
        for my $g (grep { !/^\./ } readdir $d) { my %seen; open my $gh, "<", "$o{'guides-dir'}/$g" or next; $ng++;
            while (<$gh>) { chomp; my @f = split /\t/; next unless ($f[4] // '') eq 'CDS'; (my $n = $f[0]) =~ s/_exon\d+$//; $seen{$n} = 1 } close $gh; $gsup{$_}++ for keys %seen } closedir $d }
    my @lin = grep { $nsp{$_} && $nsp{$_} >= 10 } ("F:" . ($o{family} // 'NA'), "O:" . ($o{order} // 'NA'));
    if (@lin) { my $lin = $lin[0];
        for my $i (@feat) { next if $drop{$i}; my $n = $rows[$i][0]; next if $n =~ /^(trn|rrn)/ || $n =~ /_intron\d*$/; (my $g = $n) =~ s/_exon\d+$//;
            my $share = ($with{$g}{$lin} // 0) / $nsp{$lin};
            my $gfrac = $ng ? ($gsup{$g} // 0) / $ng : 0;
            my $local = ($nsp_local{n} // 0) >= 2 ? ($spgene{$g} // 0) / $nsp_local{n} : undef;   # share of the guides' species that carry the gene
            my $support_ok = defined $local ? $local > 0.5 : $gfrac > 0.5;
            if ($share < $o{'min-presence'} && !$support_ok) { $drop{$i} = 1; printf $log "PRESENCE\tdrop %s %d-%d %s\t%s: %d of %d species (%.2f); %s\n", $n, $rows[$i][1], $rows[$i][2], $rows[$i][3], $lin, $with{$g}{$lin} // 0, $nsp{$lin}, $share, (defined $local ? sprintf("nearest species %d of %d", $spgene{$g} // 0, $nsp_local{n}) : sprintf("guides %d of %d", $gsup{$g} // 0, $ng // 0)) }
            elsif ($share < $o{'min-presence'}) { printf $log "PRESENCE\tkeep %s %d-%d %s\t%s: %d of %d species (%.2f) but %s\n", $n, $rows[$i][1], $rows[$i][2], $rows[$i][3], $lin, $with{$g}{$lin} // 0, $nsp{$lin}, $share, (defined $local ? sprintf("%d of %d nearest species carry it", $spgene{$g} // 0, $nsp_local{n}) : sprintf("guides %d of %d carry it", $gsup{$g} // 0, $ng // 0)) } }
    } else { print $log "PRESENCE\tskipped: no lineage with >= 10 species\n" }
}

# --- 3. tRNA window refinement with ARAGORN ---------------------------------------
my %moved;
if ($o{aragorn} && -x $o{aragorn} && $seq) {
    my %conv; if ($o{'trna-convention'} && open my $cf, "<", $o{'trna-convention'}) { while (<$cf>) { next if /^#/; my ($id, $off) = split; $conv{$id} = [$1, $2] if $off && $off =~ m{^(-?\d+)/(-?\d+)$} } close $cf }   # residuals file
    # curated exon lengths per spliced identity (trna_exon_lengths.tsv beside the convention file)
    my %exlen; if ($o{'trna-convention'}) { (my $ef = $o{'trna-convention'}) =~ s{[^/]*$}{trna_exon_lengths.tsv}; if (open my $eh, "<", $ef) { while (<$eh>) { next if /^#/; my ($id, $l) = split; $exlen{$id} = [$1, $2] if $l && $l =~ m{^(\d+)/(\d+)$} } close $eh } }
    my @ar; if (open my $ah, "-|", $o{aragorn}, '-t', '-i', '-w', '-gcbact', $fsa) {
        my %aa3 = (Ala=>'A',Arg=>'R',Asn=>'N',Asp=>'D',Cys=>'C',Gln=>'Q',Glu=>'E',Gly=>'G',His=>'H',Ile=>'I',Leu=>'L',Lys=>'K',Met=>'M',Phe=>'F',Pro=>'P',Ser=>'S',Thr=>'T',Trp=>'W',Tyr=>'Y',Val=>'V',fMet=>'fM');
        while (<$ah>) { next unless /^\d+\s+tRNA-(\w+)\s+(c?)\[(\d+),(\d+)\]\s+\d+\s+\((\w+)\)(.*)/; my ($aa, $c, $gs, $ge, $ac, $rest) = ($1, $2, $3, $4, $5, $6);
            (my $acu = uc $ac) =~ tr/T/U/; my $id = 'trn' . ($aa3{$aa} // $aa) . "-$acu";
            push @ar, { s => $gs, e => $ge, d => ($c ? '-' : '+'), intron => ($rest =~ /i\(\d+,\d+\)/ ? 1 : 0), id => $id } } close $ah }
    # group tRNA rows into genes: exon rows of one name and strand within 3.5 kb
    # exons pair by tRNA LETTER (de-duplication can leave one exon named trnG-UCC and
    # its partner trnG); the pair then takes the anticodon-bearing name if either has it
    my @tg; for my $i (@feat) { next if $drop{$i}; my $r = $rows[$i]; next unless $r->[0] =~ /^trn/; (my $n = $r->[0]) =~ s/_exon\d+$//; my $ex = $r->[0] =~ /_exon\d+$/ ? 1 : 0; (my $letter = $n) =~ s/-[ACGU]{3}$//; my $grp;
        if ($ex) { for my $g (@tg) { next unless $g->{ex} && $g->{letter} eq $letter && $g->{d} eq $r->[3] && @{ $g->{rows} } < 2; if (grep { abs($rows[$_][1] - $r->[1]) < 3500 } @{ $g->{rows} }) { $grp = $g; last } } }
        unless ($grp) { $grp = { n => $n, letter => $letter, d => $r->[3], rows => [], ex => $ex }; push @tg, $grp } push @{ $grp->{rows} }, $i; $grp->{n} = $n if $n =~ /-[ACGU]{3}$/ }
    for my $g (grep { $_->{ex} && @{ $_->{rows} } == 2 && $_->{n} =~ /-[ACGU]{3}$/ } @tg) { for my $i (@{ $g->{rows} }) { my $ex = $rows[$i][0] =~ /(_exon\d+)$/ ? $1 : ''; $rows[$i][0] = $g->{n} . $ex; $moved{$i} = 1 } }
    my %spliced_id_fallback = (trnI => 'trnI-GAU', trnA => 'trnA-UGC', trnK => 'trnK-UUU', trnG => 'trnG-UCC', trnL => 'trnL-UAA', trnV => 'trnV-UAC');
    for my $g (@tg) { my @r = map { $rows[$_] } @{ $g->{rows} }; next if $g->{ex} && @r != 2;
        my ($s0, $e0) = (1e12, 0); for my $f (@r) { $s0 = $f->[1] if $f->[1] < $s0; $e0 = $f->[2] if $f->[2] > $e0 }
        my ($best, $bo) = (undef, 0); for my $x (@ar) { next unless $x->{d} eq $g->{d} && $x->{intron} == ($g->{ex} ? 1 : 0); my $lo = $s0 > $x->{s} ? $s0 : $x->{s}; my $hi = $e0 < $x->{e} ? $e0 : $x->{e}; my $ov = $hi - $lo + 1; ($best, $bo) = ($x, $ov) if $ov > $bo }
        # no ARAGORN gene (it misses some spliced tRNAs): a spliced call still gets its
        # exon lengths normalised, anchored to its own outer ends
        if (!($best && $bo >= 0.5 * ($e0 - $s0 + 1))) { next unless $g->{ex}; $best = { s => $s0, e => $e0, id => $g->{n} }; $bo = $e0 - $s0 + 1; my $nm = $g->{n} =~ /-[ACGU]{3}$/ ? $g->{n} : ($spliced_id_fallback{ $g->{n} } // $g->{n}); $best->{s} = $s0; $best->{e} = $e0; $best->{noaragorn} = 1 }
        # identity for the convention: the call's own name when it carries an anticodon;
        # else ARAGORN's (reliable for unspliced genes), else the one spliced identity
        # a letter can have (trnI-GAU, trnA-UGC, trnK-UUU, trnG-UCC, trnL-UAA, trnV-UAC)
        my %spliced_id = (trnI => 'trnI-GAU', trnA => 'trnA-UGC', trnK => 'trnK-UUU', trnG => 'trnG-UCC', trnL => 'trnL-UAA', trnV => 'trnV-UAC');
        my $id = $g->{n} =~ /-[ACGU]{3}$/ ? $g->{n} : $g->{ex} ? ($spliced_id{ $g->{n} } // $g->{n}) : $best->{id};
        my ($ns, $ne) = $best->{noaragorn} ? ($s0, $e0) : trna_place(\$seq, $best->{s}, $best->{e}, $g->{d}, $id, \%conv);
        my $keep_ends = ($ns == $s0 && $ne == $e0); next if $keep_ends && !$g->{ex};
        if ($g->{ex}) { my ($r1, $r2) = $g->{d} eq '+' ? (sort { $a->[1] <=> $b->[1] } @r) : (sort { $b->[1] <=> $a->[1] } @r);
            my ($l1, $l2) = $exlen{$id} ? @{ $exlen{$id} } : ($r1->[2] - $r1->[1] + 1, $r2->[2] - $r2->[1] + 1);
            next if $keep_ends && $l1 == $r1->[2] - $r1->[1] + 1 && $l2 == $r2->[2] - $r2->[1] + 1;
            my @new = $g->{d} eq '+' ? ([$ns, $ns + $l1 - 1], [$ne - $l2 + 1, $ne]) : ([$ne - $l1 + 1, $ne], [$ns, $ns + $l2 - 1]);
            printf $log "TRNA_WINDOW\t%s %d-%d,%d-%d %s\t-> %d-%d,%d-%d\n", $g->{n}, $r1->[1], $r1->[2], $r2->[1], $r2->[2], $g->{d}, $new[0][0], $new[0][1], $new[1][0], $new[1][1];
            ($r1->[1], $r1->[2]) = @{ $new[0] }; ($r2->[1], $r2->[2]) = @{ $new[1] }; $moved{$_} = 1 for @{ $g->{rows} } }
        else { printf $log "TRNA_WINDOW\t%s %d-%d %s\t-> %d-%d\n", $r[0][0], $s0, $e0, $g->{d}, $ns, $ne; ($r[0][1], $r[0][2]) = ($ns, $ne); $moved{ $g->{rows}[0] } = 1 } }
}

# --- rewrite, merging intergenic rows around removed features ------------------
unless (%drop || %moved) { close $log if $o{log}; exit 0 }
my @keep = grep { !$drop{$_} } @feat;
my @out;
# feature rows in file order, skipping dropped ones and all old intergenic rows
for my $i (0 .. $#rows) { next if $drop{$i}; next if $rows[$i][0] =~ /~/; push @out, $rows[$i] }
# regenerate intergenic rows between consecutive kept features (by start)
my @sorted = sort { $rows[$a][1] <=> $rows[$b][1] || $rows[$a][2] <=> $rows[$b][2] } @keep;
my @gaps; my $prev_name = 'Start'; my $prev_end = 0;
for my $i (@sorted) { my ($n, $s, $e) = @{ $rows[$i] }[0 .. 2];
    push @gaps, ["$prev_name~$n", $prev_end + 1, $s - 1, '+'] if $s > $prev_end + 1;
    if ($e > $prev_end) { $prev_end = $e; $prev_name = $n } }
push @gaps, ["$prev_name~End", $prev_end + 1, $len, '+'] if $len && $len > $prev_end;
# place gap rows before the region rows, in coordinate order with the features
my @regions = grep { $_->[0] =~ /^(LSC|SSC|IRA|IRB|FULL)$/ } @out; my @features = grep { $_->[0] !~ /^(LSC|SSC|IRA|IRB|FULL)$/ } @out;
my @merged = sort { $a->[1] <=> $b->[1] || ($a->[0] =~ /~/ ? 0 : 1) <=> ($b->[0] =~ /~/ ? 0 : 1) } (@features, @gaps);
open my $w, ">", $ann or die; print $w join("\t", @$_), "\n" for @merged, @regions; close $w;
printf $log "removed %d feature rows, moved %d tRNA rows\n", scalar keys %drop, scalar keys %moved; close $log if $o{log};

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
