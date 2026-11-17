#!/usr/bin/perl
# relabel_guide_trnas.pl - name and orient the guides' tRNA genes from their sequence.
#
# USAGE: relabel_guide_trnas.pl <filedir> <seqdir> <trna_library.fasta> [--blast-bin DIR] [--exclude name,acc,...]
#
# GenBank records swap tRNA names and, now and then, strands. Every Liliaceae guide of a
# Lilium genome names the unspliced trnG-GCC "trnG-UCC" and the spliced trnG-UCC
# "trnG-GCC"; Fritillaria walujewii NC_037215 annotates the spliced trnV-UAC on the wrong
# strand, with its exons numbered for that strand. The post filter renames the TARGET's
# calls by sequence, but the guides' names still decide which calls are duplicates, what
# the validator expects (it proposed deleting an exon of a correct trnG-UCC because "the
# guides carry one trnG-UCC of 71 nt"), and a wrong-strand guide exon displaces the real
# one (trnV-UAC exon 1 was missing from 14 of 20 Lilium annotations).
#
# Each guide tRNA gene (a single row, or the exon rows of one name and strand within
# 3.5 kb, joined) is BLASTed against the curated library (tables/trna_library.fasta, one
# entry per identity per curated genome; --exclude leaves a benchmark genome's own entries
# out). With the post filter's rule - alignment over >= 80% of the gene at >= 75% identity,
# the best other identity at <= 90% of the top bit score - the gene takes the top identity's
# name, and when the top hit is on the minus strand the gene is turned over: strand flipped,
# and for a spliced gene the exons renumbered in the new transcript order.
#
# Each change is reported on stderr: TRNA_RELABEL / TRNA_STRAND <guide> <old> -> <new> <coords>.
use strict; use warnings; use File::Temp qw(tempdir);
my ($filedir, $seqdir, $lib, @rest) = @ARGV;
die "usage: $0 <filedir> <seqdir> <trna_library.fasta> [--blast-bin DIR] [--exclude a,b]\n" unless $lib && -d $filedir && -d $seqdir;
my %o; while (@rest) { my $k = shift @rest; $k =~ s/^--//; $o{$k} = shift @rest }
my $blastn = $o{'blast-bin'} ? "$o{'blast-bin'}/blastn" : $ENV{BLAST_BIN} ? "$ENV{BLAST_BIN}/blastn" : (grep { -x $_ } map { "$_/blastn" } split /:/, $ENV{PATH})[0];
unless ($blastn && -x $blastn && -s $lib) { warn "relabel_guide_trnas.pl: no blastn or no library; guide tRNA names not checked\n"; exit 0 }
my $dir = tempdir(CLEANUP => 1);
my %excl = map { $_ => 1 } grep { length } split /,/, ($o{exclude} // $ENV{ANNOBTD_TRNA_LIBRARY_EXCLUDE} // '');
{ open my $lf, "<", $lib or die "$lib: $!"; open my $lo, ">", "$dir/lib.fa" or die; my $skip = 0;
  while (<$lf>) { if (/^>(\S+)/) { my @h = split /\|/, $1; $skip = (grep { $excl{$_} } @h) ? 1 : 0 } print $lo $_ unless $skip } close $lf; close $lo }
sub rc { my $x = reverse shift; $x =~ tr/ACGT/TGCA/; $x }

opendir my $dh, $filedir or die "$filedir: $!"; my @guides = sort grep { !/^\./ && -f "$filedir/$_" && -f "$seqdir/$_" } readdir $dh; closedir $dh;
for my $g (@guides) {
    open my $f, "<", "$seqdir/$g" or next; my $seq = ''; while (<$f>) { next if /^>/; s/\s//g; $seq .= uc } close $f;
    open my $ah, "<", "$filedir/$g" or next; my @l = <$ah>; close $ah; my @row = map { chomp(my $x = $_); [split /\t/, $x, -1] } @l;
    # genes: an unspliced row, or exon rows of one base name and strand within 3.5 kb
    my @gene; for my $i (0 .. $#row) { my $r = $row[$i]; next unless @$r >= 5 && $r->[4] eq 'tRNA' && $r->[0] =~ /^trn/ && $r->[1] =~ /^\d+$/ && $r->[2] =~ /^\d+$/ && $r->[2] <= length $seq;
        if ($r->[0] =~ /^(.+)_exon\d+$/) { my $base = $1; my ($grp) = grep { $_->{ex} && $_->{base} eq $base && $_->{d} eq $r->[3] && @{ $_->{rows} } < 2 && abs($row[ $_->{rows}[0] ][1] - $r->[1]) < 3500 } @gene;
            if ($grp) { push @{ $grp->{rows} }, $i } else { push @gene, { base => $base, d => $r->[3], ex => 1, rows => [$i] } } }
        else { push @gene, { base => $r->[0], d => $r->[3], ex => 0, rows => [$i] } } }
    open my $q, ">", "$dir/q.fa" or die; my $nq = 0;
    for my $k (0 .. $#gene) { my $ge = $gene[$k]; next if $ge->{ex} && @{ $ge->{rows} } != 2;
        my $s = join '', map { substr($seq, $_->[1] - 1, $_->[2] - $_->[1] + 1) } sort { $a->[1] <=> $b->[1] } map { $row[$_] } @{ $ge->{rows} }; $s = rc($s) if $ge->{d} eq '-';
        next if length($s) < 50 || length($s) > 110; print $q ">$k\n$s\n"; $nq++ } close $q;
    next unless $nq; my %hit;   # gene -> identity -> [bits, strand]
    for (`$blastn -task blastn -word_size 7 -dust no -evalue 1e-3 -query $dir/q.fa -subject $dir/lib.fa -outfmt '6 qseqid sseqid pident length qlen bitscore sstrand' 2>/dev/null`) {
        chomp; my ($k, $sid, $pid, $alen, $qlen, $bits, $str) = split /\t/; next if $alen < 0.8 * $qlen || $pid < 75; my ($id) = split /\|/, $sid;
        $hit{$k}{$id} = [$bits, $str] if !$hit{$k}{$id} || $bits > $hit{$k}{$id}[0] }
    my $changed = 0;
    for my $k (sort { $a <=> $b } keys %hit) { my $ge = $gene[$k]; my @ids = sort { $hit{$k}{$b}[0] <=> $hit{$k}{$a}[0] || $a cmp $b } keys %{ $hit{$k} };
        my ($best, $second) = @ids; next if $second && $hit{$k}{$second}[0] > 0.9 * $hit{$k}{$best}[0];
        my @r = map { $row[$_] } @{ $ge->{rows} }; my $where = join(',', map { "$_->[1]-$_->[2]" } @r) . " $ge->{d}";
        if ($hit{$k}{$best}[1] eq 'minus') { my $nd = $ge->{d} eq '+' ? '-' : '+'; $_->[3] = $nd for @r;
            if ($ge->{ex}) { my @o = sort { $a->[1] <=> $b->[1] } @r; @o = reverse @o if $nd eq '-'; for my $n (0 .. $#o) { $o[$n][0] =~ s/_exon\d+$/'_exon' . ($n + 1)/e } }
            print STDERR "TRNA_STRAND\t$g\t$ge->{base} $where\t-> $nd (the gene matches $best on the other strand)\n"; $ge->{d} = $nd; $changed = 1 }
        next if $ge->{base} eq $best;
        for my $x (@r) { my $ex = $x->[0] =~ /(_exon\d+)$/ ? $1 : ''; $x->[0] = $best . $ex }
        printf STDERR "TRNA_RELABEL\t%s\t%s -> %s\t%s\t(%.0f bits; runner-up %s %.0f)\n", $g, $ge->{base}, $best, $where, $hit{$k}{$best}[0], $second // '-', $second ? $hit{$k}{$second}[0] : 0; $changed = 1 }
    next unless $changed; open my $w, ">", "$filedir/$g" or die "$filedir/$g: $!"; print $w join("\t", @$_), "\n" for @row; close $w }
