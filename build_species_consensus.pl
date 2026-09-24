#!/usr/bin/perl
# One vote per independent research group, then one form per species.
#
# GenBank holds many records of the same species, and the big series come from
# single labs (Hordeum vulgare: 1,432 records, 1,395 from one institute). Counting
# records lets one group's annotation methodology define the "usual" form of a
# gene. So: RefSeq copies are merged with the INSDC record they are identical to;
# records are grouped by submitter (last author + institution); each group gets
# one vote per feature (its own modal value); the species' form of each feature
# is the value most groups voted for. Every unit is then scored against its
# species' multi-group consensus, and the unit agreeing best is the species'
# representative.
#
# USAGE: build_species_consensus.pl <record_meta.tsv> <lengths.tsv> <out_prefix>
#   record_meta.tsv: extract_record_meta.pl output
#   lengths.tsv:     extract_guide_lengths.pl output (accession, name, exons, nt)
# OUT: <prefix>_units.tsv     accession species group units_in_species groups_in_species compared mismatches is_rep merged_accessions
#      <prefix>_lengths.tsv   rep_accession name exons nt   (the species' modal form; one row per species per feature)
#      <prefix>_consensus.tsv species feature modal groups_for groups_total units_for units_total
use strict; use warnings;
my ($meta, $lens, $prefix, @rest) = @ARGV; die "usage: $0 <record_meta.tsv> <lengths.tsv> <out_prefix> [--exclude taxonomy_consistency.tsv]\n" unless $prefix;
# --exclude: check_taxonomy_consistency.pl output; MISLABELED records (sequence from
# another order than the declared family) neither vote nor represent a species.
my %excl; while (@rest) { my $a = shift @rest; if ($a eq "--exclude") { my $f = shift @rest; open my $x, "<", $f or die; while (<$x>) { chomp; my @c = split /\t/; $excl{$c[0]} = 1 if ($c[10] // "") eq "MISLABELED" } close $x } }
printf STDERR "excluding %d mislabeled records\n", scalar keys %excl if %excl;

my (%rec, %base2acc);
open my $m, "<", $meta or die; <$m>;
while (<$m>) { chomp; my ($acc, $org, $sp, $src, $date, $first, $last, $aff) = split /\t/; next if $excl{$acc};
    (my $b = $acc) =~ s/\.\d+$//; $base2acc{$b} = $acc;
    $rec{$acc} = { sp => $sp, src => $src, date => $date, last => $last, aff => $aff } }
close $m;

# units: a RefSeq copy joins its cached source's unit
my (%unit_of, %unit);
for my $acc (sort keys %rec) {
    my $r = $rec{$acc}; my $uid = $acc;
    if ($r->{src} ne '-') { (my $sb = $r->{src}) =~ s/\.\d+$//; $uid = $base2acc{$sb} if $base2acc{$sb} }
    $unit_of{$acc} = $uid; push @{ $unit{$uid}{members} }, $acc; $unit{$uid}{sp} //= $r->{sp};
}
sub institution { my $a = shift; return '' if !defined $a || $a eq '-'; $a =~ s/^Contact:\S+\s+\S+\s+//; my ($i) = split /,/, $a; $i //= ''; $i =~ s/^\s+|\s+$//g; return lc $i }
my %months = (JAN=>1,FEB=>2,MAR=>3,APR=>4,MAY=>5,JUN=>6,JUL=>7,AUG=>8,SEP=>9,OCT=>10,NOV=>11,DEC=>12);
sub datekey { my $d = shift; return 0 unless $d && $d =~ /^(\d+)-([A-Z]{3})-(\d{4})$/; return $3 * 10000 + ($months{$2} || 0) * 100 + $1 }
for my $uid (keys %unit) {
    my $u = $unit{$uid}; my @mem = @{ $u->{members} };
    my ($last, $inst, $date) = ('-', '', 0);
    for my $a (@mem) { my $r = $rec{$a}; $last = $r->{last} if $last eq '-' && $r->{last} ne '-'; my $i = institution($r->{aff}); $inst = $i if $inst eq '' && $i ne ''; my $dk = datekey($r->{date}); $date = $dk if $dk > $date }
    $u->{group} = $last eq '-' ? "acc:$uid" : ($inst ne '' ? "$last|$inst" : $last);
    $u->{date} = $date;
    # the accession whose annotation represents the unit: the RefSeq copy when there is one
    my ($nc) = grep { /^NC_/ } @mem; $u->{acc} = $nc // $mem[0];
    $u->{merged} = join(",", grep { $_ ne $u->{acc} } @mem) || '-';
}
my %acc2unit; $acc2unit{ $unit{$_}{acc} } = $_ for keys %unit;

# features of each unit's representative accession
my %F;
open my $l, "<", $lens or die;
while (<$l>) { chomp; my ($acc, $name, $exons, $nt) = split /\t/; next unless $acc2unit{$acc};
    my $f = $F{ $acc2unit{$acc} } ||= {};
    $f->{"len:$name"} = $nt;                     # IR copies: same length, last one wins
    $f->{"ex:$name"} = $exons if $name !~ /_exon\d+$/;
}
close $l;

my %by_sp; push @{ $by_sp{ $unit{$_}{sp} } }, $_ for keys %unit;
open my $U, ">", "${prefix}_units.tsv" or die; open my $L, ">", "${prefix}_lengths.tsv" or die; open my $C, ">", "${prefix}_consensus.tsv" or die;
print $U join("\t", qw(accession species group units_in_species groups_in_species compared mismatches is_rep merged_accessions)), "\n";
print $L join("\t", qw(accession name exons nt)), "\n";
print $C join("\t", qw(species feature modal groups_for groups_total units_for units_total)), "\n";
sub modal { my ($h, $tie) = @_; my ($v) = sort { $h->{$b} <=> $h->{$a} || ($tie ? ($tie->{$b} // 0) <=> ($tie->{$a} // 0) : 0) || $a <=> $b } keys %$h; return $v }
for my $sp (sort keys %by_sp) {
    my @units = @{ $by_sp{$sp} }; my %groups; push @{ $groups{ $unit{$_}{group} } }, $_ for @units;
    my $ng = scalar keys %groups;
    # per group, per feature: the group's modal value
    my (%gvote, %ucount);
    for my $g (keys %groups) {
        my %vals; for my $uid (@{ $groups{$g} }) { my $f = $F{$uid} or next; for my $k (keys %$f) { $vals{$k}{ $f->{$k} }++; $ucount{$k}{ $f->{$k} }++ } }
        for my $k (keys %vals) { $gvote{$k}{ modal($vals{$k}) }++ }
    }
    my %modal;
    for my $k (keys %gvote) {
        my $gt = 0; $gt += $_ for values %{ $gvote{$k} }; my $ut = 0; $ut += $_ for values %{ $ucount{$k} };
        my $v = modal($gvote{$k}, $ucount{$k}); $modal{$k} = [$v, $gt];
        printf $C "%s\t%s\t%s\t%d\t%d\t%d\t%d\n", $sp, $k, $v, $gvote{$k}{$v}, $gt, $ucount{$k}{$v} // 0, $ut;
    }
    # score units against the multi-group consensus (features with >= 2 groups voting)
    my @scored;
    for my $uid (@units) { my $f = $F{$uid} || {}; my ($n, $bad) = (0, 0);
        for my $k (keys %$f) { my $mm = $modal{$k} or next; next if $mm->[1] < 2; $n++; $bad++ if $f->{$k} != $mm->[0] }
        push @scored, [$uid, $n, $bad] }
    my ($rep) = sort { $a->[2] <=> $b->[2] || ($unit{$b->[0]}{acc} =~ /^NC_/) <=> ($unit{$a->[0]}{acc} =~ /^NC_/) || $unit{$b->[0]}{date} <=> $unit{$a->[0]}{date} || $a->[0] cmp $b->[0] } @scored;
    for my $s (sort { $a->[0] cmp $b->[0] } @scored) { my $u = $unit{$s->[0]};
        printf $U "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n", $u->{acc}, $sp, $u->{group}, scalar @units, $ng, $s->[1], $s->[2], ($s->[0] eq $rep->[0] ? 1 : 0), $u->{merged} }
    # the species' form, keyed by the representative's accession (for taxonomy lookup)
    my $racc = $unit{ $rep->[0] }{acc};
    for my $k (sort keys %modal) { next unless $k =~ /^len:(.*)$/; my $name = $1;
        my $exons = $name =~ /_exon\d+$/ ? 1 : ($modal{"ex:$name"} ? $modal{"ex:$name"}[0] : 1);
        print $L join("\t", $racc, $name, $exons, $modal{$k}[0]), "\n" }
}
close $U; close $L; close $C;
printf STDERR "%d records -> %d units -> %d species\n", scalar keys %rec, scalar keys %unit, scalar keys %by_sp;
