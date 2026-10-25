#!/usr/bin/perl
# annoBTD annotation table -> NCBI 5-column feature table (.tbl) and a FASTA with source
# modifiers, the pair table2asn / tbl2asn turns into a GenBank submission (.sqn).
#
# USAGE: annobtd2tbl.pl <annotation.txt> <genome.fsa> <seqid> --out-dir DIR
#            [--organism "Agave amica"] [--modifier isolate=BNM_2022_32] [--modifier ...]
#            [--products tables/plastid_products.tsv] [--check]
#
#   seqid        the sequence identifier (no spaces); files DIR/<seqid>.fsa and .tbl are written
#   --organism   [organism=...] on the FASTA header; [location=chloroplast] [topology=circular]
#                [completeness=complete] [moltype=genomic DNA] are always added
#   --modifier   any further source modifier, repeatable (isolate, specimen_voucher, country, ...)
#   --check      translate every CDS and report start codon, internal stops, length (what the
#                .val files would say), on stderr
#
# Features written: gene + CDS (exons as a join) with /product, /codon_start 1; a CDS whose
# start is not ATG or GTG gets /exception "RNA editing" and /transl_except (pos:..,aa:Met)
# as the McKain-lab records do; rps12 is two CDS joining the LSC exon with each IR copy's
# exons 2 and 3, /exception "trans-splicing"; gene + tRNA or rRNA with /product derived
# from the name; a *_pseudogene row becomes a gene with /pseudo and a note. Introns and
# spacer rows are not features.
use strict; use warnings; use FindBin;
my ($ann, $fsa, $seqid, @rest) = @ARGV; die "usage: $0 <annotation.txt> <genome.fsa> <seqid> --out-dir DIR [--organism X] [--modifier k=v]... [--products FILE] [--check]\n" unless $seqid;
my (%o, @mods); while (@rest) { my $a = shift @rest; if ($a eq '--modifier') { push @mods, shift @rest } elsif ($a eq '--check') { $o{check} = 1 } elsif ($a =~ /^--(.+)$/) { $o{$1} = shift @rest } }
$o{'out-dir'} //= '.'; mkdir $o{'out-dir'} unless -d $o{'out-dir'};
my $prodf = $o{products} // "$FindBin::Bin/../tables/plastid_products.tsv"; my %prod; if (open my $ph, "<", $prodf) { while (<$ph>) { chomp; next if /^#/ || !/\S/; my ($g, $p) = split /\t/; $prod{$g} = $p } close $ph }
my %aa3 = (A=>'Ala',C=>'Cys',D=>'Asp',E=>'Glu',F=>'Phe',G=>'Gly',H=>'His',I=>'Ile',K=>'Lys',L=>'Leu',M=>'Met',N=>'Asn',P=>'Pro',Q=>'Gln',R=>'Arg',S=>'Ser',T=>'Thr',V=>'Val',W=>'Trp',Y=>'Tyr',fM=>'Met');
my %CODON; { my @b = qw(T C A G); my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; my $i = 0; for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1) } } } }
my $seq = ''; open my $fh, "<", $fsa or die "$fsa: $!"; while (<$fh>) { next if /^>/; s/\s//g; $seq .= uc $_ } close $fh; my $len = length $seq;
sub cds_seq { my (@iv) = @_; my $d = $iv[0][2]; my @s = sort { $a->[0] <=> $b->[0] } @iv; my $x = join '', map { substr($seq, $_->[0] - 1, $_->[1] - $_->[0] + 1) } @s; if ($d eq '-') { $x = reverse $x; $x =~ tr/ACGT/TGCA/ } $x }
sub translate { my $s = shift; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr($s, $i, 3)} // 'X' } $p }
# read features: name -> list of [start, end, strand]; exon rows grouped by base name
my (%rows, @order); open my $ah, "<", $ann or die "$ann: $!";
while (<$ah>) { chomp; next unless /\S/; my @f = split /\t/; next if $f[0] =~ /~/ || $f[0] =~ /^(LSC|SSC|IRA|IRB|FULL)$/ || $f[0] =~ /_intron\d*$/; (my $base = $f[0]) =~ s/_exon\d+$//; push @order, $base unless $rows{$base}; push @{ $rows{$base} }, [$f[1], $f[2], $f[3], ($f[0] =~ /_exon(\d+)$/ ? $1 : 0)] } close $ah;
# group a name's rows into genes: same strand, exons within 6 kb of each other; IR copies are separate genes
my @genes;   # [name, [intervals in genome order], strand, is_exonic]
for my $name (@order) { my @r = sort { $a->[0] <=> $b->[0] } @{ $rows{$name} };
    if ($name eq 'rps12') { push @genes, map { ['rps12', $_] } @r; next }   # handled below
    my @grp; for my $x (@r) { my $g = $grp[-1]; if ($g && $g->[0][2] eq $x->[2] && $x->[3] && $g->[-1][3] && abs($x->[0] - $g->[-1][1]) < 6000) { push @$g, $x } else { push @grp, [$x] } }
    push @genes, [$name, $_, $_->[0][2], ($_->[0][3] ? 1 : 0)] for @grp }
my $loc = sub { my ($iv, $d) = @_; my @s = sort { $a->[0] <=> $b->[0] } @$iv; @s = reverse @s if $d eq '-'; join "\n", map { $d eq '+' ? "$_->[0]\t$_->[1]" : "$_->[1]\t$_->[0]" } @s };
my @out; my @check;
my $cds = sub { my ($name, $iv, $d, $exc) = @_; my $x = cds_seq(map { [@$_[0,1], $d] } @$iv); my $start = substr($x, 0, 3); my $p = translate($x);
    my ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\tgene", (map { $_ } @more), "\t\t\tgene\t$name";
    ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\tCDS", @more, "\t\t\tgene\t$name", "\t\t\tproduct\t" . ($prod{$name} // "hypothetical protein"), "\t\t\tcodon_start\t1";
    push @out, "\t\t\texception\ttrans-splicing" if $exc;
    if ($start !~ /^(ATG|GTG)$/) { my @s = sort { $a->[0] <=> $b->[0] } @$iv; my ($ps, $pe) = $d eq '+' ? ($s[0][0], $s[0][0] + 2) : ($s[-1][1] - 2, $s[-1][1]); push @out, "\t\t\texception\tRNA editing", "\t\t\ttransl_except\t(pos:" . ($d eq '+' ? "$ps..$pe" : "complement($ps..$pe)") . ",aa:Met)" }
    my $stops = ($p =~ tr/\*//) - (substr($p, -1) eq '*' ? 1 : 0); push @check, sprintf("%s\t%s\t%d nt\tstart %s%s%s", $name, $loc->($iv, $d) =~ s/\n/,/gr, length $x, $start, ($stops ? "\tINTERNAL STOPS $stops" : ""), (substr($p, -1) ne '*' ? "\tNO TERMINAL STOP" : "")) };
for my $g (@genes) { my ($name, $iv, $d, $ex) = @$g; next if $name eq 'rps12';
    if ($name =~ /^(.+)_pseudogene$/) { my $base = $1; my ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\tgene", @more, "\t\t\tgene\t$base", "\t\t\tpseudo", "\t\t\tnote\t$base pseudogene"; next }
    if ($name =~ /^trn/) { my ($letter) = $name =~ /^trn(fM|[A-Z])/; my ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\tgene", @more, "\t\t\tgene\t$name"; ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\ttRNA", @more, "\t\t\tgene\t$name", "\t\t\tproduct\ttRNA-" . ($aa3{$letter} // 'Xxx'); next }
    if ($name =~ /^rrn(\S+)/) { my $u = $1; my ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\tgene", @more, "\t\t\tgene\t$name"; ($first, @more) = split /\n/, $loc->($iv, $d); push @out, "$first\trRNA", @more, "\t\t\tgene\t$name", "\t\t\tproduct\t${u}S ribosomal RNA"; next }
    $cds->($name, $iv, $d, 0) }
# rps12: exon 1 (one copy, usually in the LSC) joined with exons 2+3 of each IR copy
if ($rows{rps12}) { my @r = @{ $rows{rps12} }; my @e1 = grep { $_->[3] == 1 } @r; my @e23 = grep { $_->[3] >= 2 } @r;
    my @copies; for my $x (sort { $a->[0] <=> $b->[0] } @e23) { my $c = $copies[-1]; if ($c && $c->[0][2] eq $x->[2] && abs($x->[0] - $c->[-1][1]) < 3000) { push @$c, $x } else { push @copies, [$x] } }
    for my $c (@copies) { my $d = $c->[0][2]; my @iv = (@e1, @$c); @iv = grep { defined } @iv; my $trans = (@e1 && $e1[0][2] ne $d) ? 1 : 1;   # always trans-spliced
        # exon 1 may lie on the other strand: list it first in transcript order, with its own strand
        my ($first, @more) = split /\n/, $loc->([@e1], $e1[0][2]); my @loc23 = split /\n/, $loc->($c, $d);
        push @out, "$first\tgene", @more, @loc23, "\t\t\tgene\trps12", "\t\t\texception\ttrans-splicing";
        push @out, "$first\tCDS", @more, @loc23, "\t\t\tgene\trps12", "\t\t\tproduct\t" . ($prod{rps12} // 'ribosomal protein S12'), "\t\t\tcodon_start\t1", "\t\t\texception\ttrans-splicing";
        my $x = cds_seq([@{$e1[0]}[0,1], $e1[0][2]]) . cds_seq(map { [@$_[0,1], $d] } @$c); my $p = translate($x); my $stops = ($p =~ tr/\*//) - (substr($p,-1) eq '*' ? 1 : 0);
        push @check, sprintf("rps12\t%s+%s\t%d nt\tstart %s%s%s", $loc->([@e1], $e1[0][2]) =~ s/\n/,/gr, $loc->($c, $d) =~ s/\n/,/gr, length $x, substr($x,0,3), ($stops ? "\tINTERNAL STOPS $stops" : ""), (substr($p,-1) ne '*' ? "\tNO TERMINAL STOP" : "")) } }
# write
my $dir = $o{'out-dir'}; open my $tb, ">", "$dir/$seqid.tbl" or die; print $tb ">Feature $seqid\n", join("\n", @out), "\n"; close $tb;
my $hdr = ">$seqid" . ($o{organism} ? " [organism=$o{organism}]" : "") . " [location=chloroplast] [topology=circular] [completeness=complete] [moltype=genomic DNA]" . join('', map { " [$_]" } @mods);
open my $fo, ">", "$dir/$seqid.fsa" or die; print $fo "$hdr\n"; for (my $i = 0; $i < $len; $i += 70) { print $fo substr($seq, $i, 70), "\n" } close $fo;
if ($o{check}) { my $bad = grep { /INTERNAL STOPS|NO TERMINAL STOP/ } @check; printf STDERR "%s: %d CDS checked, %d with problems\n", $seqid, scalar @check, $bad; print STDERR "  $_\n" for grep { /INTERNAL STOPS|NO TERMINAL STOP|start (?!ATG|GTG)/ } @check }
