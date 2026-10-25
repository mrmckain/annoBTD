#!/usr/bin/perl
# Score validate_annotation.pl verdicts against hand edits.
# USAGE: score_validation.pl <validation_dir> <manual_edits.tsv>
use strict; use warnings;
my ($dir, $edits) = @ARGV; die "usage: $0 <validation_dir> <manual_edits.tsv>\n" unless $edits;
# truth: genome -> list of { gene, s, e, want (REAL|ARTIFACT), action, ns, ne }
my (%truth, %chain);
open my $eh, "<", $edits or die; while (<$eh>) { chomp; next if /^#/ || !/\S/; my ($g, $act, $name, $s, $e, @x) = split /\t/; (my $base = $name) =~ s/_exon\d+$//; $base =~ s/_pseudogene$//;
    if ($act eq 'move' && $name =~ /_pseudogene$/) { $chain{"$g\t$base"} = [@x[0, 1]]; next }   # the span given to a pseudogene just renamed
    my $want = ($act eq 'rename' && $x[0] =~ /_pseudogene$/) || ($act eq 'add' && $name =~ /_pseudogene$/) ? 'REAL' : $act eq 'delete' || $act eq 'move' ? 'ARTIFACT' : 'OTHER';
    next if $act eq 'delete' && ($name eq 'infA' || $name =~ /_pseudogene$/);   # paired with the add of a pseudogene, or the withdrawal of an earlier call
    push @{ $truth{$g} }, { gene => $base, s => $s, e => $e, want => $want, act => $act, ns => $x[0], ne => $x[1], line => $_ } } close $eh;
my ($hit, $miss, %fp, %fpl, $n) = (0, 0);
my @missed;
opendir(my $dh, $dir) or die; for my $f (sort grep { /\.tsv$/ } readdir $dh) { (my $g = $f) =~ s/\.tsv$//; open my $vh, "<", "$dir/$f" or next; my @v; while (<$vh>) { chomp; my @c = split /\t/; next if $c[0] eq 'feature'; push @v, { name => $c[0], s => $c[1], e => $c[2], st => $c[4], flags => $c[5], ev => $c[6], edit => $c[7] // '' } } close $vh; $n++;
    my %used;
    for my $t (@{ $truth{$g} || [] }) { my $ok;
        for my $i (0 .. $#v) { my $r = $v[$i]; (my $rb = $r->{name}) =~ s/ \(.*//; next unless $rb eq $t->{gene}; next if $r->{st} eq 'PASS';
            my $over = $t->{act} eq 'add' ? 1 : ($r->{s} <= $t->{e} && $r->{e} >= $t->{s}); next unless $over;
            $used{$i} = 1; my $coord = ''; if ($t->{act} eq 'move' && $r->{edit} =~ /^move \S+ \d+ \d+ (\d+) (\d+)/) { $coord = ($1 == $t->{ns} && $2 == $t->{ne}) ? ' coords=same' : " coords=$1-$2 vs $t->{ns}-$t->{ne}" }
            elsif ($t->{act} eq 'add' && $r->{edit} =~ /^add \S+ (\d+) (\d+)/) { $coord = (abs($1 - $t->{s}) <= 30 && abs($2 - $t->{e}) <= 30) ? ' span~same' : " span=$1-$2 vs $t->{s}-$t->{e}" }
            $ok = ($r->{st} eq $t->{want}) ? "HIT  $g $t->{gene} $t->{act}: $r->{st}$coord" : "WRONG $g $t->{gene} $t->{act}: got $r->{st}, wanted $t->{want} ($r->{ev})"; last }
        if ($ok && $ok =~ /^HIT/) { $hit++; print "$ok\n" } else { $miss++; push @missed, $ok // "MISS  $g $t->{gene} $t->{act}: no verdict (PASS)" } }
    for my $i (0 .. $#v) { next if $used{$i} || $v[$i]{st} eq 'PASS'; (my $rb = $v[$i]{name}) =~ s/ \(.*//; $fp{"$rb\t$v[$i]{st}"}++; push @{ $fpl{"$rb\t$v[$i]{st}"} }, "$g: $v[$i]{flags} | $v[$i]{ev}" } }
print "$_\n" for @missed;
printf "\n%d genomes scored; hand edits found with the right verdict: %d of %d\n", $n, $hit, $hit + $miss;
printf "verdicts not matching a hand edit (feature, status, count):\n"; for my $k (sort { $fp{$b} <=> $fp{$a} || $a cmp $b } keys %fp) { printf "  %-28s %3d   e.g. %s\n", $k, $fp{$k}, substr($fpl{$k}[0], 0, 160) }
