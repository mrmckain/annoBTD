#!/usr/bin/perl
# relabel_guides.pl - correct guide genes that GenBank records name wrongly, by vote
# of the other guides.
#
# USAGE: relabel_guides.pl <filedir> <seqdir> [--blast-bin DIR]
#
# Every single-exon CDS of every guide is translated and compared (blastp, no
# low-complexity filter) with the single-exon CDSs of the other guides. For each
# copy, a gene name's support is the number of other guides carrying a copy of
# that name whose protein matches it (>= 70% identity over >= 80% of both
# proteins). A copy whose own name has less support than another name is renamed
# to the best-supported one, in place in <filedir>. Ties are left alone.
#
# The case it was written for: 40 Liliaceae records (Tulipa, Cardiocrinum,
# Calochortus, Clintonia, Gagea, Medeola, Notholirion) annotate psbJ as a second
# "psaJ", so a guide pool with three of them named the target's psbJ psaJ and the
# duplicate filter then dropped it. When most of the guides share a mislabel, the
# vote follows them; nothing here can see that.
#
# Each rename is reported on stderr: RELABEL <guide> <old> -> <new> <start>-<end> <strand> (support).
use strict; use warnings;
use File::Temp qw(tempdir);

my ($filedir, $seqdir, @rest) = @ARGV;
die "usage: $0 <filedir> <seqdir> [--blast-bin DIR]\n" unless $seqdir && -d $filedir && -d $seqdir;
my %o; while (@rest) { my $k = shift @rest; $k =~ s/^--//; $o{$k} = shift @rest }
my $blastp = $o{'blast-bin'} ? "$o{'blast-bin'}/blastp" : $ENV{BLAST_BIN} ? "$ENV{BLAST_BIN}/blastp" : (grep { -x $_ } map { "$_/blastp" } split /:/, $ENV{PATH})[0];
unless ($blastp && -x $blastp) { warn "relabel_guides.pl: no blastp; guide names not checked\n"; exit 0 }

my %CODON; { my @b = qw(T C A G); my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; my $i = 0;
    for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1) } } } }
sub translate { my $s = shift; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr($s, $i, 3)} // 'X' } $p =~ s/\*$//; $p }

opendir my $dh, $filedir or die "$filedir: $!"; my @guides = sort grep { !/^\./ && -f "$filedir/$_" && -f "$seqdir/$_" } readdir $dh; closedir $dh;
exit 0 if @guides < 3;
my (%rows, %prot);   # guide -> [lines]; "guide\tline" -> protein
for my $g (@guides) {
    open my $f, "<", "$seqdir/$g" or die; my $seq = ''; while (<$f>) { next if /^>/; s/\s//g; $seq .= uc } close $f;
    open my $a, "<", "$filedir/$g" or die; my @l = <$a>; close $a; $rows{$g} = \@l;
    for my $i (0 .. $#l) { my @c = split /\t/, $l[$i]; next unless @c >= 5 && $c[4] eq 'CDS' && $c[0] !~ /_exon\d+$/ && $c[1] =~ /^\d+$/ && $c[2] =~ /^\d+$/;
        my ($s, $e) = $c[1] < $c[2] ? @c[1, 2] : @c[2, 1]; next if $e > length $seq || $e - $s + 1 < 60;
        my $n = substr($seq, $s - 1, $e - $s + 1); if ($c[3] eq '-') { $n = reverse $n; $n =~ tr/ACGT/TGCA/ }
        my $p = translate($n); next if length $p < 20 || ($p =~ tr/*//) > 1;
        $prot{"$g\t$i"} = $p }
}
my $dir = tempdir(CLEANUP => 1); my @ids = sort keys %prot; my %idx; @idx{0 .. $#ids} = @ids;
open my $q, ">", "$dir/p.fa" or die; print $q ">$_\n$prot{$idx{$_}}\n" for 0 .. $#ids; close $q;
my %support;   # copy id -> name -> {guide => 1}
for (`$blastp -query $dir/p.fa -subject $dir/p.fa -seg no -evalue 1e-5 -outfmt '6 qseqid sseqid pident length qlen slen' 2>/dev/null`) {
    chomp; my ($qi, $si, $pid, $len, $ql, $sl) = split /\t/; next if $qi eq $si || $pid < 70 || $len < 0.8 * $ql || $len < 0.8 * $sl;
    my ($qg) = split /\t/, $idx{$qi}; my ($sg, $sl_i) = split /\t/, $idx{$si}; next if $qg eq $sg;
    my ($name) = split /\t/, $rows{$sg}[$sl_i]; $support{$qi}{$name}{$sg} = 1 }
my %changed;
for my $id (sort { $a <=> $b } keys %support) {
    my ($g, $i) = split /\t/, $idx{$id}; my @c = split /\t/, $rows{$g}[$i]; my $own = $c[0];
    my %n = map { $_ => scalar keys %{ $support{$id}{$_} } } keys %{ $support{$id} }; my $mine = $n{$own} // 0;
    my ($best, $second) = sort { $n{$b} <=> $n{$a} || $a cmp $b } grep { $_ ne $own } keys %n; next unless $best;
    next unless $n{$best} > $mine && $n{$best} > ($second ? $n{$second} : 0);   # a clear winner only
    $c[0] = $best; $rows{$g}[$i] = join "\t", @c; $changed{$g} = 1;
    printf STDERR "RELABEL\t%s\t%s -> %s\t%s-%s %s\t(%d other guides call this protein %s, %d call it %s)\n", $g, $own, $best, $c[1], $c[2], $c[3], $n{$best}, $best, $mine, $own }
for my $g (sort keys %changed) { open my $o, ">", "$filedir/$g" or die "$filedir/$g: $!"; print $o @{ $rows{$g} }; close $o }
