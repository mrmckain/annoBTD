#!/usr/bin/perl
use strict;
use warnings;

# Unit tests for the boundary-finding helpers in match_orfs_to_blast_v2.2.pl.
#
# The subs are lifted out of the source at run time and eval'd here, so production
# code needs no refactor and the tests cannot drift from a stale copy.
#
# USAGE: perl unit_boundary.pl
# Exit status is the number of failures.
#
# Tests marked KNOWN-BAD document defects found in the diagnostic. They are expected
# to fail today; when the boundary rewrite lands, flip them to expect the right answer.

my $SRC = "../match_orfs_to_blast_v2.2.pl";
open my $fh, "<", $SRC or die "Cannot open $SRC: $!";
my $src = do { local $/; <$fh> };
close $fh;

our %codons = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G','GCN'=>'A','CGN'=>'R','GGN'=>'G','CCN'=>'P','TCN'=>'S','ACN'=>'T','GTN'=>'V');

# Brace-counting extraction: some subs (match_end) close an inner block at column 0,
# so a non-greedy /^\}/ match would truncate them.
sub extract_sub {
    my ($text, $name) = @_;
    return undef unless $text =~ /^sub \Q$name\E\s*\{/m;
    my $start = $-[0];
    my $depth = 0;
    my $in_pod = 0;
    my @lines = split /\n/, substr($text, $start);
    my @out;
    for my $line (@lines) {
        $in_pod = 1 if $line =~ /^=\w+/;
        push @out, $line;
        unless ($in_pod) {
            my $code = $line;
            $code =~ s/#.*$//;
            $depth++ while $code =~ /\{/g;
            $depth-- while $code =~ /\}/g;
        }
        $in_pod = 0 if $line =~ /^=cut/;
        last if $depth == 0 && @out > 1;
    }
    return join("\n", @out);
}

for my $name (qw(match_start alt_start match_end translate best_match best_match_frames extend_rna_hit splice_score splice_donor_ok pick_short_exon locate_ref_start)) {
    my $body = extract_sub($src, $name);
    die "Could not extract sub $name from $SRC\n" unless $body;
    $body =~ s/\bmy \%codons\b/our %codons/g;
    eval "package main; no strict; no warnings; $body; 1" or die "eval of $name failed: $@";
}

my ($pass, $fail, @failures) = (0, 0);
sub check {
    my ($desc, $got, $want, $known_bad) = @_;
    my $ok = defined $got && defined $want && $got eq $want;
    if ($ok) { $pass++; printf "  ok    %s\n", $desc }
    else {
        $fail++;
        push @failures, $desc;
        printf "  %s %s\n      got %s, want %s\n",
            ($known_bad ? "KNOWN-BAD" : "FAIL "), $desc,
            (defined $got ? "'$got'" : 'undef'), (defined $want ? "'$want'" : 'undef');
    }
    return $ok;
}

print "translate()\n";
check("frame 0 of ATGAAATAA", translate("ATGAAATAA"), "MK_");
check("unknown codon becomes X", translate("ATGNNNTAA"), "MX_");

print "\nmatch_start()  -- picks the start codon for a CDS\n";
# True start is in frame at offset 0; a decoy ATG sits out of frame at offset 7.
my $orf_decoy = "ATG" . "AAACATGCC" . ("GGGTTTAAACCC" x 5);
check("prefers the in-frame ATG over an out-of-frame decoy",
      match_start($orf_decoy, 5), 0, 1);

# Two in-frame ATGs, at nt 0 and nt 9. Which one is right depends on where the
# reference N-terminus aligns, so both directions are specified.
my $orf_two = "ATG" . "AAAGGG" . "ATG" . ("CCCTTTAAAGGG" x 4);
check("reference aligns at nt 0  -> picks the ATG at 0",
      match_start($orf_two, 0), 0);
check("reference aligns at nt 15 -> picks the nearer ATG at 9",
      match_start($orf_two, 5), 9);

# An out-of-frame ATG nearer the reference position must still lose to an in-frame one.
my $orf_near = "ATG" . "AAAGG" . "ATG" . ("CCCTTTAAAGGG" x 4);   # decoy ATG at nt 8
check("in-frame ATG wins even when an out-of-frame ATG is closer to the target",
      match_start($orf_near, 4), 0);

print "\nalt_start()  -- non-ATG start codons (ACG/GTG are real in plastid CDS)\n";
my $gene_acg = "ACG" . "AAAGGGCCC" x 3;
my $orf_acg  = "ACG" . "AAAGGGCCC" x 3;
check("finds an in-frame ACG start", alt_start($orf_acg, 5, $gene_acg), 0);

print "\nmatch_end()  -- 3' boundary, measured on the nucleotide scale\n";
# Protein ends ...WY then stop. match_end returns an offset back from the ORF end.
my $orf_prot  = "MKAILVPQWY_";
my $gene_prot = "MKAILVPQWY_";
{
    my $got = match_end($orf_prot, $gene_prot, 5);
    check("identical ORF and reference give a zero-length 3' trim", $got, 0, 1);
}

{
    # match_end measures on 3*residues, so its result never includes a trailing
    # partial codon. The callers add length($orf_seq)%3 to put it on the same
    # scale as their stop-codon path, which measures on the nucleotide sequence.
    my $prot = "MKAILVPQWY";        # 10 residues -> 30 nt of full codons
    my $got = match_end($prot, $prot, 5);
    check("the result is a whole number of codons back from the protein end",
          $got % 3, 0);
    check("identical protein and reference put the end at the last residue",
          $got, 0);
}

print "\nbest_match_frames()  -- frame selection\n";
{
    # Reference protein, and an ORF whose correct reading is frame 1.
    my $gene_prot2 = "MKAILVPQWYCDEFGH";
    my $f0 = "XXXXXXXXXXXXXXXX";
    my $f1 = $gene_prot2;
    my $f2 = "ZZZZZZZZZZZZZZZZ";
    my ($fortot, $for_frame, $revtot, $rev_frame) =
        best_match_frames($gene_prot2, $f0, $f2, $f1, $f2, $f2, $f2);
    check("selects forward frame 1 when only that frame matches", $for_frame, 1);
    check("forward score beats reverse", ($fortot > $revtot ? 1 : 0), 1);
}

print "\nextend_rna_hit()  -- walk a tRNA/rRNA HSP out to the full reference gene\n";
{
    # Reference gene is 100 nt. Plus-strand hit covering reference 3..98,
    # sitting at plastome 1003..1098. True gene is plastome 1001..1100.
    my ($s,$e) = extend_rna_hit(1003, 1098, 3, 98, 100, 150000);
    check("plus strand: both ends walked out", "$s-$e", "1001-1100");

    # Minus-strand hit: qs pairs with the HIGH subject coord.
    # Reference 98..3 across plastome 1003..1098 -> true gene 1001..1100.
    ($s,$e) = extend_rna_hit(1003, 1098, 98, 3, 100, 150000);
    check("minus strand: ends walked out the other way", "$s-$e", "1001-1100");

    # A full-length hit must not move.
    ($s,$e) = extend_rna_hit(1001, 1100, 1, 100, 100, 150000);
    check("full-length hit is left alone", "$s-$e", "1001-1100");
    ($s,$e) = extend_rna_hit(1001, 1100, 100, 1, 100, 150000);
    check("full-length minus hit is left alone", "$s-$e", "1001-1100");

    # Asymmetric: only the 5' end of the reference is missing.
    ($s,$e) = extend_rna_hit(1003, 1100, 3, 100, 100, 150000);
    check("only the unaligned end moves", "$s-$e", "1001-1100");

    # Clamping at the ends of a linear molecule.
    ($s,$e) = extend_rna_hit(2, 97, 3, 98, 100, 150000);
    check("clamps at position 1", "$s-$e", "1-99");
    ($s,$e) = extend_rna_hit(149903, 149998, 3, 98, 100, 150000);
    check("clamps at the genome end", "$s-$e", "149901-150000");
}

print "\nsplice_score()  -- plastid group II introns are GT...AY, not GT...AG\n";
{
    our $plastome;
    # 96 nt body, not self-similar under reverse complement; intron spans 101..200.
    my $body = "ACGGTTCAAG" x 9 . "ACGGTT";
    my $flank = "N" x 100;

    $plastome = $flank . "GT" . $body . "AC" . $flank;
    check("GT...AC scores full marks", splice_score(101, 200, "+"), 4);

    $plastome = $flank . "GT" . $body . "AT" . $flank;
    check("GT...AT scores full marks", splice_score(101, 200, "+"), 4);

    # Tobacco rpl2 is genuinely GT...AA: still a real intron, scored below AY.
    $plastome = $flank . "GT" . $body . "AA" . $flank;
    check("GT...AA still scores (tobacco rpl2 is this)", splice_score(101, 200, "+"), 3);

    $plastome = $flank . "GT" . $body . "AG" . $flank;
    check("GT...AG scores donor plus a weak acceptor", splice_score(101, 200, "+"), 3);

    $plastome = $flank . "TT" . $body . "AC" . $flank;
    check("no donor scores acceptor only", splice_score(101, 200, "+"), 2);

    $plastome = $flank . "TT" . $body . "GG" . $flank;
    check("no signal at either end scores 0", splice_score(101, 200, "+"), 0);

    # Strand: GT...AT discriminates (GT...AC is palindromic at the ends).
    my $intron = "GT" . $body . "AT";
    my $rc = reverse $intron; $rc =~ tr/ACGT/TGCA/;
    $plastome = $flank . $rc . $flank;
    check("minus strand is judged on the transcribed strand",
          splice_score(101, 200, "-"), 4);
    check("...and that same span scores worse read as plus",
          (splice_score(101, 200, "+") < 4 ? 1 : 0), 1);

    $plastome = "GTAAAC";
    check("a span under 10 nt is not an intron", splice_score(1, 6, "+"), -1);

    # The donor is the only signal allowed to move a boundary.
    $plastome = $flank . "GT" . $body . "GG" . $flank;
    check("donor test passes on GT even with a dead acceptor",
          splice_donor_ok(101, 200, "+"), 1);
    $plastome = $flank . "TT" . $body . "AC" . $flank;
    check("donor test fails on TT even with a perfect acceptor",
          splice_donor_ok(101, 200, "+"), 0);
}

print "\npick_short_exon()  -- a 6 nt exon matches dozens of places\n";
{
    our $plastome;
    # Exon at 1-based 1001-1006 followed by a GT intron running to 1800, where
    # exon 2 starts. A decoy copy of the exon sits at 1051, closer to exon 2 but
    # followed by CT - the real Spinacia petB situation.
    my $filler = "AAACCCGGGTTT";
    $plastome  = ($filler x 83) . "";                        # 996 nt of filler
    $plastome .= "ATGAGT";                                   # exon at 997..1002 (1-based)
    $plastome .= "GT" . ("ACGTTGCA" x 5);                    # GT donor
    $plastome .= "ATGAGT";                                   # decoy at 1045..1050
    $plastome .= "CT" . ("ACGTTGCA" x 90);                   # CT, no donor
    my $real  = 996;    # 0-based start of the true exon
    my $decoy = 1044;   # 0-based start of the decoy
    my $anchor = 1780;  # 0-based start of exon 2, further from the real exon
    my $got = pick_short_exon([$real, $decoy], 6, $anchor, "+");
    check("prefers the GT-donor candidate over the nearer one", $got, $real);

    # With no donor anywhere, fall back to plain proximity.
    $plastome =~ s/^(.{1002})GT/$1CT/s;
    my $got2 = pick_short_exon([$real, $decoy], 6, $anchor, "+");
    check("falls back to nearest when no candidate has a donor", $got2, $decoy);

    check("a single candidate is returned unchanged",
          pick_short_exon([$real], 6, $anchor, "+"), $real);
}

print "\nlocate_ref_start()  -- anchor the reference N-terminus specifically\n";
{
    # The reference's residues 1.. appear late in the ORF, but its first three
    # residues also occur early by chance. The long probe must win.
    my $gene_prot = "MKAILVPQWYCDEFGH";
    my $orf_prot  = "ZZKAI" . ("Z" x 20) . $gene_prot . "ZZZZ";
    my $want = index($orf_prot, $gene_prot);
    check("a 3-mer decoy does not capture the anchor",
          locate_ref_start($orf_prot, $gene_prot), $want);
    check("no match at all returns -1",
          locate_ref_start("QQQQQQQQQQQQ", $gene_prot), -1);
}

print "\n";
printf "%d passed, %d failed\n", $pass, $fail;
if (@failures) {
    print "failing:\n";
    print "  - $_\n" for @failures;
}
exit($fail > 255 ? 255 : $fail);
