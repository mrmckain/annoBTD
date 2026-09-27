#!/usr/bin/perl
use strict; use warnings;
# Build a taxonomy table for every complete chloroplast genome in GenBank, in the
# same 7-column layout as taxonomy_cache.tsv (tribe/subfamily left NA):
#   accession organism taxid tribe subfamily family order
# USAGE: build_taxonomy_all.pl <out.tsv> [existing.tsv]   (existing rows are kept)
my ($out, $existing) = @ARGV; die "usage\n" unless $out;
my $E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $TERM = '(chloroplast[Title] OR plastid[Title]) AND complete genome[Title] AND 80000:300000[SLEN]';
sub get { my $u = shift; my $r; for (1..4) { $r = `curl -g -sS -m 120 "$u"`; last if length $r; sleep 2 } select(undef,undef,undef,0.4); $r }
my %have;
if ($existing && -s $existing) { open my $h, "<", $existing or die; while (<$h>) { chomp; s/\r//g; my @f = split /\t/; next if $. == 1; (my $a = $f[0]) =~ s/\.\d+$//; $have{$a} = $_ } close $h }
if (-s $out) { open my $h, "<", $out or die; while (<$h>) { chomp; s/\r//g; my @f = split /\t/; next if $. == 1; (my $a = $f[0]) =~ s/\.\d+$//; $have{$a} //= $_ } close $h }
# 1. search with history
my $term = $TERM; $term =~ s/([^A-Za-z0-9_.\-])/sprintf("%%%02X",ord($1))/ge;
my $x = get("$E/esearch.fcgi?db=nuccore&term=$term&usehistory=y&retmax=0");
my ($cnt) = $x =~ /<Count>(\d+)/; my ($qk) = $x =~ /<QueryKey>(\d+)/; my ($we) = $x =~ /<WebEnv>(\S+?)</;
die "esearch failed\n" unless $cnt && $we;
print STDERR "esearch: $cnt records\n";
# 2. accession list
my @acc;
for (my $rs = 0; $rs < $cnt; $rs += 10000) {
    my $t = get("$E/efetch.fcgi?db=nuccore&query_key=$qk&WebEnv=$we&rettype=acc&retmode=text&retstart=$rs&retmax=10000");
    $t =~ tr/\r//d; push @acc, grep { /\S/ } split /\n/, $t;
}
print STDERR scalar(@acc), " accessions listed\n";
my @need = grep { (my $b = $_) =~ s/\.\d+$//; !$have{$b} } @acc;
print STDERR scalar(@need), " need taxonomy\n";
# 3. esummary -> taxid, title
my (%tax, %title);
for (my $i = 0; $i < @need; $i += 300) {
    my @b = @need[$i .. ($i + 299 < $#need ? $i + 299 : $#need)];
    my $t = get("$E/esummary.fcgi?db=nuccore&id=" . join(',', @b) . "&retmode=xml");
    while ($t =~ /<DocSum>(.*?)<\/DocSum>/sg) {
        my $d = $1; my ($av) = $d =~ /Name="AccessionVersion"[^>]*>([^<]+)/; my ($ti) = $d =~ /Name="TaxId"[^>]*>(\d+)/; my ($tt) = $d =~ /Name="Title"[^>]*>([^<]+)/;
        next unless $av && $ti; $tax{$av} = $ti; $title{$av} = $tt // '';
    }
    print STDERR "  esummary $i/", scalar(@need), "\n" if $i % 6000 == 0;
}
# 4. taxonomy -> organism, family, order
my %lin; my %uniq = map { $_ => 1 } values %tax; my @tids = sort keys %uniq;
for (my $i = 0; $i < @tids; $i += 300) {
    my @b = @tids[$i .. ($i + 299 < $#tids ? $i + 299 : $#tids)];
    my $t = get("$E/efetch.fcgi?db=taxonomy&id=" . join(',', @b) . "&retmode=xml");
    my ($depth, $cur, $inlin, $rank, $name) = (0);
    $t =~ tr/\r//d; for my $line (split /\n/, $t) {
        if ($line =~ /<Taxon>/) { $depth++; if ($depth == 1) { $cur = { fam => 'NA', ord => 'NA' } } }
        if ($depth == 1 && !$inlin) {
            if ($line =~ /<TaxId>(\d+)</ && !$cur->{id}) { $cur->{id} = $1 }
            if ($line =~ /<ScientificName>([^<]+)</ && !$cur->{org}) { $cur->{org} = $1 }
        }
        $inlin = 1 if $line =~ /<LineageEx>/; $inlin = 0 if $line =~ /<\/LineageEx>/;
        if ($inlin) { if ($line =~ /<ScientificName>([^<]+)</) { $name = $1 } if ($line =~ /<Rank>(\w+)</) { $cur->{fam} = $name if $1 eq 'family'; $cur->{ord} = $name if $1 eq 'order' } }
        if ($line =~ /<\/Taxon>/) { $depth--; if ($depth == 0 && $cur && $cur->{id}) { $lin{ $cur->{id} } = $cur; $cur = undef } }
    }
}
print STDERR scalar(keys %lin), " taxids resolved\n";
open my $o, ">>", $out or die; print $o join("\t", qw(accession organism taxid tribe subfamily family order)), "\n" unless -s $out;
my $n = 0;
for my $av (@need) { my $ti = $tax{$av} or next; my $l = $lin{$ti} || {}; my $org = $l->{org} // $title{$av}; $org =~ s/\s+chloroplast.*//i;
    print $o join("\t", $av, $org, $ti, 'NA', 'NA', $l->{fam} // 'NA', $l->{ord} // 'NA'), "\n"; $n++ }
close $o; print STDERR "wrote $n new rows to $out\n";
