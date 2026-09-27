#!/usr/bin/perl
# Does each record's SEQUENCE agree with its declared family? A mislabeled or
# contaminated GenBank record (MZ127828.1 "Bidens pilosa" is an Arabidopsis
# plastome) would be picked as a guide for the wrong lineage and would vote in
# the wrong profile cells.
#
# All-against-all Jaccard over 67k sketches is 2 billion comparisons, so:
#   1. family core = hashes present in >= 50% of the family's sketches but in
#      <= MAX_GLOBAL (10%) of all sketches - plastomes are so conserved that
#      universal k-mers would otherwise be every family's core (families with
#      >= MIN_MEMBERS records);
#   2. each record is scored against every core by containment,
#      |sketch ∩ core| / |core|; its best-scoring family is compared to the
#      declared one;
#   3. a record whose best family is not its declared family, with margin, is
#      verified by true Jaccard against up to 5 members of each family;
#   4. a record whose family has no core (fewer than MIN_MEMBERS records, or no
#      distinctive hashes) is tested the same way at ORDER level, against order
#      cores (hashes in >= 50% of the order, <= 15% of all sketches).
#
# USAGE: check_taxonomy_consistency.pl <sketch_db.tsv> <out.tsv> [--min-members 5] [--margin 2.0] [--min-best 0.4]
# OUT:   accession organism family order best_family best_order c_declared c_best
#        j_declared j_best verdict level
#        verdict: OK | UNTESTED | REVIEW | MISLABELED | ASSIGNABLE; level: family | order | none.
#        ASSIGNABLE: a record with NO declared lineage whose best family is unambiguous
#        (containment >= min-best, runner-up <= half of it) and verified by Jaccard
#        >= 0.15 to a member of that family; best_family/best_order are the proposed
#        lineage (assign_lineages.pl applies them to the taxonomy table).
#        At order level best_family/best_order hold the best ORDER and c_*/j_* are
#        order-level scores; REVIEW cannot occur (there is no finer rank to dispute).
use strict; use warnings;
my ($db, $out, @rest) = @ARGV; die "usage: $0 <sketch_db.tsv> <out.tsv> [--min-members N] [--margin X] [--min-best X]\n" unless $out;
my %o = ('min-members' => 5, margin => 2.0, 'min-best' => 0.4, 'max-global' => 0.10);
while (@rest) { my $a = shift @rest; $o{$1} = shift @rest if $a =~ /^--(.+)$/ }

my (@acc, @org, @fam, @ord, @sk);   # sketches packed as N* to fit 67k x 2000 in memory
open my $h, "<", $db or die; my $i = 0;
while (<$h>) { chomp; my ($a, $o_, $f, $r, $len, $hl) = split /\t/; next unless defined $hl;
    push @acc, $a; push @org, $o_; push @fam, $f; push @ord, $r; push @sk, pack("N*", split /,/, $hl); $i++ }
close $h; printf STDERR "%d sketches\n", scalar @acc;

# global document frequency of every hash, then family cores of distinctive hashes
my %df; for my $j (0 .. $#acc) { $df{$_}++ for unpack "N*", $sk[$j] }
my $gmax = $o{'max-global'} * @acc; printf STDERR "%d distinct hashes\n", scalar keys %df;
my (%members, %forder); push @{ $members{ $fam[$_] } }, $_ for 0 .. $#acc; $forder{ $fam[$_] } //= $ord[$_] for 0 .. $#acc;
my %core;   # family -> packed core hashes
for my $f (sort keys %members) { my @m = @{ $members{$f} }; next if $f eq 'NA' || @m < $o{'min-members'};
    my %cnt; for my $j (@m) { $cnt{$_}++ for unpack "N*", $sk[$j] }
    my $need = 0.5 * @m; my @c = sort { $a <=> $b } grep { $cnt{$_} >= $need && $df{$_} <= $gmax } keys %cnt;
    next if @c < 20;   # a family with no distinctive hashes cannot be tested
    @c = (sort { $df{$a} <=> $df{$b} || $a <=> $b } @c)[0 .. 199] if @c > 200;   # the 200 most distinctive
    $core{$f} = pack "N*", @c }
printf STDERR "%d family cores, mean size %.0f\n", scalar keys %core, eval { my $t = 0; $t += length($core{$_}) / 4 for keys %core; $t / (keys %core || 1) };
my @cf = sort keys %core; my %csize = map { $_ => length($core{$_}) / 4 } @cf;
# order cores, for records whose family cannot be tested
my %omembers; push @{ $omembers{ $ord[$_] } }, $_ for 0 .. $#acc;
my %ocore; my $omax = 0.15 * @acc;
for my $o_ (sort keys %omembers) { my @m = @{ $omembers{$o_} }; next if $o_ eq 'NA' || @m < $o{'min-members'};
    my %cnt; for my $j (@m) { $cnt{$_}++ for unpack "N*", $sk[$j] }
    my $need = 0.5 * @m; my @c = sort { $a <=> $b } grep { $cnt{$_} >= $need && $df{$_} <= $omax } keys %cnt;
    next if @c < 20; @c = (sort { $df{$a} <=> $df{$b} || $a <=> $b } @c)[0 .. 199] if @c > 200;
    $ocore{$o_} = pack "N*", @c }
my @co = sort keys %ocore; my %osize = map { $_ => length($ocore{$_}) / 4 } @co;
printf STDERR "%d order cores\n", scalar @co;

sub jaccard { my ($x, $y) = @_; my %s; $s{$_} = 1 for unpack "N*", $x; my $n = 0; $n++ for grep { $s{$_} } unpack "N*", $y; my $u = 4000 - $n; return $u ? $n / $u : 0 }
my (%sample, %osample);   # family/order -> up to 5 member indexes spread through the list
for my $f (@cf) { my @m = @{ $members{$f} }; my $step = int(@m / 5) || 1; my @s; for (my $k = 0; $k < @m && @s < 5; $k += $step) { push @s, $m[$k] } $sample{$f} = \@s }
for my $o_ (@co) { my @m = @{ $omembers{$o_} }; my $step = int(@m / 5) || 1; my @s; for (my $k = 0; $k < @m && @s < 5; $k += $step) { push @s, $m[$k] } $osample{$o_} = \@s }

open my $w, ">", $out or die;
print $w join("\t", qw(accession organism family order best_family best_order c_declared c_best j_declared j_best verdict level)), "\n";
my ($ntest, $nrev, $nmis, $notest, $nomis) = (0, 0, 0, 0, 0);
for my $r (0 .. $#acc) {
    my %q; $q{$_} = 1 for unpack "N*", $sk[$r];
    my ($best, $cb, $cd, $second) = ('-', 0, 0, 0);
    for my $f (@cf) { my $n = 0; $n++ for grep { $q{$_} } unpack "N*", $core{$f}; my $c = $n / $csize{$f};
        $cd = $c if $f eq $fam[$r]; if ($c > $cb) { ($second, $best, $cb) = ($cb, $f, $c) } elsif ($c > $second) { $second = $c } }
    my ($verdict, $jd, $jb, $level) = ('UNTESTED', 'NA', 'NA', 'family');
    if (!$core{ $fam[$r] }) {
        # order level
        $level = 'none';
        if ($ocore{ $ord[$r] }) { $level = 'order'; ($best, $cb, $cd) = ('-', 0, 0);
            for my $o_ (@co) { my $n = 0; $n++ for grep { $q{$_} } unpack "N*", $ocore{$o_}; my $c = $n / $osize{$o_};
                $cd = $c if $o_ eq $ord[$r]; if ($c > $cb) { ($best, $cb) = ($o_, $c) } }
            if ($best eq $ord[$r] || $cb < $o{'min-best'} || $cb < $o{margin} * $cd) { $verdict = 'OK' }
            else { my ($mb, $md) = (0, 0);
                for my $j (@{ $osample{$best} }) { next if $j == $r; my $v = jaccard($sk[$r], $sk[$j]); $mb = $v if $v > $mb }
                for my $j (@{ $osample{ $ord[$r] } }) { next if $j == $r; my $v = jaccard($sk[$r], $sk[$j]); $md = $v if $v > $md }
                ($jd, $jb) = (sprintf("%.3f", $md), sprintf("%.3f", $mb));
                $verdict = ($mb >= $o{margin} * $md && $mb >= 0.05) ? 'MISLABELED' : 'OK'; $nomis++ if $verdict eq 'MISLABELED' }
            $notest++ }
        $verdict = 'UNTESTED' if $level eq 'none';
        # no declared lineage at all: propose one when the sequence says so unambiguously
        if ($fam[$r] eq 'NA' && $ord[$r] eq 'NA' && $best ne '-' && $cb >= $o{'min-best'} && $second <= 0.5 * $cb) {
            my $mb = 0; for my $j (@{ $sample{$best} }) { next if $j == $r; my $v = jaccard($sk[$r], $sk[$j]); $mb = $v if $v > $mb }
            $jb = sprintf("%.3f", $mb); if ($mb >= 0.15) { $verdict = 'ASSIGNABLE'; $level = 'none' } } }
    elsif ($best eq $fam[$r] || $cb < $o{'min-best'} || $cb < $o{margin} * $cd) { $verdict = 'OK' }
    else {
        # verify by true Jaccard: nearest of 5 members of each family (excluding self)
        my ($mb, $md) = (0, 0);
        for my $j (@{ $sample{$best} }) { next if $j == $r; my $v = jaccard($sk[$r], $sk[$j]); $mb = $v if $v > $mb }
        for my $j (@{ $sample{ $fam[$r] } }) { next if $j == $r; my $v = jaccard($sk[$r], $sk[$j]); $md = $v if $v > $md }
        ($jd, $jb) = (sprintf("%.3f", $md), sprintf("%.3f", $mb));
        $verdict = ($mb >= $o{margin} * $md && $mb >= 0.05) ? ($forder{$best} eq $ord[$r] ? 'REVIEW' : 'MISLABELED') : 'OK';
        $nrev++ if $verdict eq 'REVIEW'; $nmis++ if $verdict eq 'MISLABELED';
    }
    $ntest++ if $verdict ne 'UNTESTED';
    printf $w "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%s\t%s\t%s\t%s\n", $acc[$r], $org[$r], $fam[$r], $ord[$r], $best, ($level eq 'order' ? $best : $forder{$best} // '-'), $cd, $cb, $jd, $jb, $verdict, $level;
    printf STDERR "  %d/%d\n", $r + 1, scalar @acc if ($r + 1) % 5000 == 0;
}
close $w; printf STDERR "tested %d (family), REVIEW %d, MISLABELED %d; order-level tested %d, MISLABELED %d\n", $ntest - $notest, $nrev, $nmis, $notest, $nomis;
