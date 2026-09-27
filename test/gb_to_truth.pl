#!/usr/bin/perl
use strict;
use FindBin;
use warnings;

# Build a ground-truth annotation table from a GenBank flatfile.
#
# USAGE: gb_to_truth.pl <file.gb> > truth.tsv
# OUT:   gene <TAB> start <TAB> end <TAB> strand <TAB> type <TAB> flags
#        1-based inclusive forward-strand coordinates, one row per exon.
#        Multi-exon features become <gene>_exon1..N numbered in TRANSCRIPTION order.
#        flags is a comma list, currently "trans_spliced" or "-".
#
# Handles the cases genbank2verdant.pl does not:
#   * strand from complement()
#   * complement(join(A,B,C))   -> transcription order is C,B,A
#   * join(complement(A),B,C)   -> per-segment strands, file order is transcription order
#   * rRNA features carrying only /product, no /gene
#   * wrapped locations and wrapped qualifier values

my $gb = shift or die "usage: $0 <file.gb>\n";
open my $fh, "<", $gb or die "Cannot open $gb: $!";

# RefSeq plastid rRNAs are usually annotated with /product only.
my @RRNA_MAP = (
    [ qr/\b16S\b/i        => 'rrn16'  ],
    [ qr/\b23S\b/i        => 'rrn23'  ],
    [ qr/\b4\.?5S\b/i     => 'rrn4.5' ],
    [ qr/\b5S\b/i         => 'rrn5'   ],
);

my %AA3 = (
    Ala=>'A', Arg=>'R', Asn=>'N', Asp=>'D', Cys=>'C', Gln=>'Q', Glu=>'E',
    Gly=>'G', His=>'H', Ile=>'I', Leu=>'L', Lys=>'K', Met=>'M', Phe=>'F',
    Pro=>'P', Ser=>'S', Thr=>'T', Trp=>'W', Tyr=>'Y', Val=>'V',
    fMet=>'fM', Sec=>'U', OTHER=>'X', Xle=>'X', Undet=>'X',
);

# annoBTD expects trnX-YYY. Records vary: /gene="trnA-UGC", /gene="trnA", or no
# /gene at all with /product="tRNA-Ala (UGC)". Normalise all of them, and never
# emit a name containing whitespace - downstream parsers split on it.
sub normalize_gene {
    my ($name, $type, $product) = @_;
    $name = '' unless defined $name;
    $product = '' unless defined $product;

    # Plastid gene renames. Records deposited at different times use different
    # names for the same gene, and a guide named pafI cannot match a target named
    # ycf3 - the reference is simply never found.
    # The table lives in annoBTD/gene_synonyms.tsv (shared with extract_guide_lengths.pl
    # so profile cells merge the same way); the built-in list is the fallback.
    $name = canonical_gene($name);

    # rrn16 / rrn16S / rrn4.5S are all the same gene across records.
    if ($name =~ /^rrn([\d.]+)S?$/i) {
        return "rrn$1";
    }

    if ($type eq 'tRNA') {
        my $anticodon;
        $anticodon = uc $1 if $product =~ /\(\s*([ACGUTacgut]{3})\s*\)/;
        $anticodon = uc $1 if !$anticodon && $name =~ /-\s*([ACGUTacgut]{3})\s*$/;
        $anticodon =~ tr/T/U/ if $anticodon;

        my $letter;
        if ($name =~ /^trn([A-Za-z]{1,2})(?:-|$)/ && $name !~ /^tRNA/i) { $letter = $1 }
        elsif (($name =~ /tRNA-(\w+)/i ? $1 : ($product =~ /tRNA-(\w+)/i ? $1 : '')) =~ /^(\w+)$/) {
            my $aa = $1;
            $letter = $AA3{$aa} || $AA3{ucfirst lc $aa};
        }
        if ($letter) {
            return $anticodon ? "trn$letter-$anticodon" : "trn$letter";
        }
    }
    $name =~ s/\s+//g;
    return $name;
}

my (@rows, %seen);
my $skipped_unnamed = 0;
my ($key, $loc, @quals);

sub qual {
    my $want = shift;
    for my $q (@quals) { return $1 if $q =~ /^\/\Q$want\E="?([^"]*)"?/ }
    return undef;
}

sub flush {
    return unless defined $key;
    my $type = { CDS => 'CDS', tRNA => 'tRNA', rRNA => 'rRNA' }->{$key};
    return unless $type;
    return if grep { m{^/pseudo} } @quals;

    my $gene = qual('gene');
    if (!defined $gene && $type eq 'rRNA') {
        my $product = qual('product') || '';
        for my $m (@RRNA_MAP) { if ($product =~ $m->[0]) { $gene = $m->[1]; last } }
    }
    $gene = normalize_gene($gene, $type, qual('product'));
    # Features carrying only a locus_tag are hypothetical ORFs (/note="ORF105").
    # No reference set can transfer a name to them, so they are not a fair part of
    # the benchmark - and leaving them in a guide file injects junk reference names.
    unless (defined $gene && length $gene) {
        $skipped_unnamed++;
        return;
    }

    my $trans = (grep { m{^/trans_splicing} } @quals) ? 'trans_spliced' : '-';

    # Outer complement() wrapping the whole location reverses transcription order.
    my $l = $loc;
    $l =~ s/\s+//g;
    my $outer_rc = 0;
    if ($l =~ /^complement\((.*)\)$/) {
        my $inner = $1;
        # Only an outer wrapper if the parens balance across the whole inner string.
        my $d = 0; my $ok = 1;
        for my $c (split //, $inner) {
            $d++ if $c eq '('; $d-- if $c eq ')';
            if ($d < 0) { $ok = 0; last }
        }
        if ($ok && $d == 0) { $outer_rc = 1; $l = $inner }
    }

    # Segments in file order, each with its own strand.
    my @seg;
    while ($l =~ /(complement\()?[<>]?(\d+)\s*\.\.\s*[<>]?(\d+)\)?/g) {
        push @seg, { start => $2, end => $3, strand => $1 ? '-' : '+' };
    }
    return unless @seg;

    if ($outer_rc) {
        $_->{strand} = '-' for @seg;
        @seg = reverse @seg;   # GenBank lists ascending; transcription runs the other way
    }

    my $n = 0;
    for my $s (@seg) {
        $n++;
        my $name = @seg > 1 ? "${gene}_exon$n" : $gene;
        my $rowkey = join("\t", $name, $s->{start}, $s->{end}, $s->{strand});
        next if $seen{$rowkey}++;   # trans-spliced genes share an exon between two CDS records
        push @rows, [ $name, $s->{start}, $s->{end}, $s->{strand}, $type, $trans ];
    }
}

my $in_features = 0;
my $in_origin   = 0;
my $genome      = '';
while (my $line = <$fh>) {
    if ($in_origin) {
        last if $line =~ m{^//};
        $line =~ s/[^A-Za-z]//g;
        $genome .= $line;
        next;
    }
    if ($line =~ /^FEATURES/)             { $in_features = 1; next }
    if ($line =~ /^ORIGIN/)               { flush(); $in_features = 0; $in_origin = 1; next }
    if ($line =~ m{^(CONTIG|//)})         { flush(); last }
    next unless $in_features;

    if ($line =~ /^ {5}(\S+)\s+(\S.*?)\s*$/) {      # feature key in cols 6-20
        flush();
        ($key, $loc, @quals) = ($1, $2);
    }
    elsif (defined $key && $line =~ /^ {21}(\S.*?)\s*$/) {
        my $cont = $1;
        if    ($cont =~ m{^/}) { push @quals, $cont }
        elsif (@quals)         { $quals[-1] .= $cont }
        else                   { $loc .= $cont }
    }
}
flush();
close $fh;

# Flag same-named loci that carry identical sequence but differently sized
# annotations - two IR copies of one gene where the record annotates one a few
# nucleotides longer than the other. No annotator can satisfy both: a reference
# set holds one sequence per gene, so whichever extent it takes, the other copy
# scores as an error. Marking them keeps that ceiling visible instead of looking
# like an algorithm defect.
if (length $genome) {
    $genome = uc $genome;
    my %by_name;
    for my $r (@rows) {
        my $seq = substr($genome, $r->[1] - 1, $r->[2] - $r->[1] + 1);
        if ($r->[3] eq '-') { $seq = reverse $seq; $seq =~ tr/ACGT/TGCA/ }
        push @{ $by_name{ $r->[0] } }, [ $r, $seq ];
    }
    for my $name (keys %by_name) {
        my @v = @{ $by_name{$name} };
        next if @v < 2;
        for my $i (0 .. $#v) {
            for my $j ($i + 1 .. $#v) {
                my ($a, $b) = ($v[$i][1], $v[$j][1]);
                next if $a eq $b;
                my ($short, $long) = length($a) <= length($b) ? ($a, $b) : ($b, $a);
                next unless index($long, $short) == 0;
                next if length($long) - length($short) > 10;
                for my $k ($i, $j) {
                    my $f = $v[$k][0][5];
                    $v[$k][0][5] = ($f eq '-') ? 'inconsistent_copies'
                                               : "$f,inconsistent_copies";
                }
            }
        }
    }
}

die "No CDS/tRNA/rRNA features parsed from $gb - is it a GenBank flatfile?\n" unless @rows;

for my $r (sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @rows) {
    print join("\t", @$r), "\n";
}

my %genes = map { ($_->[0] =~ /^(.+)_exon\d+$/ ? $1 : $_->[0]) => 1 } @rows;
printf STDERR "%s: %d exon rows, %d distinct genes (%d unnamed ORF features skipped)\n",
    $gb, scalar(@rows), scalar(keys %genes), $skipped_unnamed;

# Canonical plastid gene symbol: synonym table (annoBTD/gene_synonyms.tsv) applied
# case-insensitively, then case-folding to the canonical spelling of a known gene
# (clpp, ClpP, CLPP -> clpP). Unknown names pass through unchanged.
{
    my (%syn, %canon, $loaded);
    sub canonical_gene {
        my $name = shift; return $name unless defined $name && length $name;
        unless ($loaded) {
            $loaded = 1;
            my %builtin = (clpp1 => 'clpP', pafi => 'ycf3', pafii => 'ycf4', pbf1 => 'psbN', lhba => 'psbZ', ycf9 => 'psbZ',
                           psb30 => 'ycf12', tic214 => 'ycf1', ycf10 => 'cemA', ycf5 => 'ccsA', ycf6 => 'petN');
            %syn = %builtin;
            for my $d (grep { defined } $ENV{ANNOBTD_DIR}, "$FindBin::Bin/..", $FindBin::Bin) {
                if (open my $fh, "<", "$d/gene_synonyms.tsv") {
                    while (<$fh>) { chomp; next if /^#/ || !/\S/; my ($a, $b) = split /\t/; $syn{lc $a} = $b if $b } close $fh; last }
            }
            $canon{lc $_} = $_ for qw(accD atpA atpB atpE atpF atpH atpI ccsA cemA chlB chlL chlN clpP infA matK ndhA ndhB ndhC ndhD
                ndhE ndhF ndhG ndhH ndhI ndhJ ndhK petA petB petD petG petL petN psaA psaB psaC psaI psaJ psaM psbA psbB psbC psbD
                psbE psbF psbH psbI psbJ psbK psbL psbM psbN psbT psbZ rbcL rpl2 rpl14 rpl16 rpl20 rpl22 rpl23 rpl32 rpl33 rpl36
                rpoA rpoB rpoC1 rpoC2 rps2 rps3 rps4 rps7 rps8 rps11 rps12 rps14 rps15 rps16 rps18 rps19 ycf1 ycf2 ycf3 ycf4 ycf12
                ycf15 ycf68 ycf66);
            $canon{lc $_} = $_ for values %syn;
        }
        my $k = lc $name; $k =~ s/\s+//g;
        return $syn{$k} if exists $syn{$k};
        return $canon{$k} if exists $canon{$k};
        # "ycf10/cemA", "ycf3/pafI": both names of one gene joined with a slash
        if ($k =~ m{/}) { my @p = map { $syn{$_} // $canon{$_} } grep { exists $syn{$_} || exists $canon{$_} } split m{/}, $k;
            return $p[0] if @p && !grep { $_ ne $p[0] } @p }
        return $name;
    }
}
