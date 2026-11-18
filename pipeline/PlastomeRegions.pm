package PlastomeRegions;
use strict; use warnings;

# PlastomeRegions - find the inverted repeat of a plastome and name its four regions.
#
#   my $r = PlastomeRegions::find($sequence);              # undef when there is no inverted repeat
#   $r->{LSC}, {IRB}, {SSC}, {IRA}   each a list of [start, end] pieces, 1-based inclusive, ascending.
#                                    One piece, or two when the region crosses the record's origin
#                                    (the first piece runs to the last base, the second starts at 1).
#   $r->{identity}                   percent identity between the two copies over the aligned columns
#   $r->{note}                       what was unusual about the record, or ''
#
# What it replaces. match_orfs_to_blast took the reverse complement of the record's last 20 nt, looked
# for it once with index(), and walked a 20-mer inward until its reverse complement no longer occurred
# anywhere. That is right for a record laid out LSC-IRb-SSC-IRa whose last base is the last base of IRa
# and whose two copies are identical, and for nothing else:
#   - a record that ends a few bases after IRa, or before it (NC_062580.1: 2 nt; NC_026787.1: 64 nt of
#     LSC follow IRa), has no reverse complement of its last 20 nt, index() returned -1 and every row
#     came out as "IRB 1-1";
#   - one substitution between the copies ended the walk (NC_029937.1: IRB 130 nt of 26,540);
#   - a record that begins in another region (NC_054182.1 begins with IRb) has no IR at its end.
# Daucus carota in the regression set was in the first group too.
#
# Method. The record is treated as a circle. Every K-mer that occurs once is looked up as its reverse
# complement; a hit at i whose partner starts at j pairs base p with base (i + j + K - 1 - p), so all
# hits of one inverted repeat share the sum i + j + K - 1 (an indel between the copies moves it by the
# indel's length). The sum with the most hits is the inverted repeat; hits within BAND of it are
# chained into runs (a substitution removes K hits, so runs are joined across gaps up to MERGE_GAP),
# the longest run is one copy and its ends are extended base by base against the other copy, through
# an isolated mismatch but not into unrelated sequence. The second copy is the mirror image of the
# first, so the two always agree. The longer single-copy region is the LSC; IRb is the copy that
# follows it in the record's reading direction, as GenBank records are laid out.
#
# An inverted repeat shorter than MIN_IR (400 nt; ANNOBTD_MIN_IR) is not reported: Pinus keeps a 495 nt
# IR, but the dispersed inverted repeats of an IR-less plastome are not regions.

our $K         = 25;
our $BAND      = 200;
our $MERGE_GAP = 400;
our $MIN_IR    = defined $ENV{ANNOBTD_MIN_IR} ? $ENV{ANNOBTD_MIN_IR} : 400;

sub rc { my $s = reverse $_[0]; $s =~ tr/ACGT/TGCA/; $s }

sub find {
    my ($seq) = @_;
    $seq = uc $seq; $seq =~ s/\s+//g;
    my $L = length $seq;
    return undef if $L < 4 * $K;
    my $ext = $seq . substr($seq, 0, $K - 1);              # circular: a K-mer may start at any base

    my %pos;                                                # K-mer -> its start, or -1 when it occurs more than once
    for (my $i = 0; $i < $L; $i++) { my $w = substr($ext, $i, $K); next if $w =~ tr/ACGT//c; $pos{$w} = exists $pos{$w} ? -1 : $i }
    my (%votes, @hits);
    for (my $i = 0; $i < $L; $i++) {
        my $w = substr($ext, $i, $K); next if $w =~ tr/ACGT//c; next if $pos{$w} < 0;
        my $j = $pos{ rc($w) }; next unless defined $j && $j >= 0 && $j != $i;
        my $s = ($i + $j + $K - 1) % $L; $votes{$s}++; push @hits, [$i, $s] }
    return undef unless @hits;
    my ($s0) = sort { $votes{$b} <=> $votes{$a} || $a <=> $b } keys %votes;
    my $near = sub { my $d = abs($_[0] - $s0); $d = $L - $d if $d > $L / 2; $d <= $BAND };
    @hits = grep { $near->($_->[1]) } @hits;                # already in order of position

    # runs of hits: [first position, last position, sum at the first hit, sum at the last hit]
    my @runs = ([$hits[0][0], $hits[0][0], $hits[0][1], $hits[0][1]]);
    for my $h (@hits[1 .. $#hits]) {
        if ($h->[0] - $runs[-1][1] > $MERGE_GAP) { push @runs, [$h->[0], $h->[0], $h->[1], $h->[1]] }
        else { @{ $runs[-1] }[1, 3] = @$h } }
    if (@runs >= 2 && $runs[0][0] + $L - $runs[-1][1] <= $MERGE_GAP) {     # one run across the origin
        my $first = shift @runs; $runs[-1] = [$runs[-1][0], $first->[1] + $L, $runs[-1][2], $first->[3]] }
    my ($run) = sort { ($b->[1] - $b->[0]) <=> ($a->[1] - $a->[0]) || $a->[0] <=> $b->[0] } @runs;

    # copy 1 is [a, b) on the circle; base p of it pairs with base (sum - p) of copy 2
    my ($a, $b, $sa, $sb) = ($run->[0], $run->[1] + $K, $run->[2], $run->[3]);
    my $base = sub { substr($seq, $_[0] % $L, 1) };
    my $pairs = sub { my ($x, $y) = ($base->($_[0]), $base->($_[1])); $x =~ tr/ACGT/TGCA/ ? $x eq $y : 0 };
    my $extend = sub {                                      # $dir -1: left from a; +1: right from b
        my ($p, $sum, $dir, $max) = @_; my ($score, $best, $n) = (0, 0, 0);
        for (my $x = 1; $x <= $max; $x++, $p += $dir) {
            $score += $pairs->($p, $sum - $p) ? 1 : -5;
            if ($score > $best) { $best = $score; $n = $x }
            last if $score < $best - 15 }
        $n };
    my $room = int(($L - 2 * ($b - $a)) / 2); $room = 0 if $room < 0;      # the copies cannot meet
    my $xa = $extend->($a - 1, $sa, -1, $room); $a -= $xa;
    my $xb = $extend->($b, $sb, 1, $room - $xa); $b += $xb;
    my ($c, $d) = ($sb - $b + 1, $sa - $a + 1);             # copy 2, the mirror image
    return undef if $b - $a < $MIN_IR || $d - $c < $MIN_IR;

    # put copy 2 after copy 1 on the circle: a < b <= c < d <= a + L
    my ($len1, $len2) = ($b - $a, $d - $c);
    $a %= $L; $b = $a + $len1; $c = $b + ($c - $b) % $L; $d = $c + $len2;
    return undef if $d > $a + $L;                           # the copies overlap: a palindrome, not an IR
    my @x = ($b, $c); my @y = ($d, $a + $L);                # the two single-copy stretches
    my %r = $x[1] - $x[0] >= $y[1] - $y[0]
        ? (LSC => \@x, IRB => [$c, $d], SSC => \@y, IRA => [$a, $b])
        : (LSC => \@y, IRB => [$a, $b], SSC => \@x, IRA => [$c, $d]);

    my ($n, $same) = (0, 0);                                # identity of the copies, ungapped along each diagonal
    for (my $p = $a; $p < $b; $p++) { my $sum = $p - $a < ($b - $a) / 2 ? $sa : $sb; $n++; $same++ if $pairs->($p, $sum - $p) }
    my @note;
    my %out = (identity => $n ? 100 * $same / $n : 0, length => $L);
    for my $name (qw(LSC IRB SSC IRA)) {
        my ($s, $e) = @{ $r{$name} }; my $len = $e - $s; $s %= $L;
        if ($len <= 0) { $out{$name} = []; push @note, "no $name: the copies of the inverted repeat are adjacent"; next }
        if ($s + $len <= $L) { $out{$name} = [[$s + 1, $s + $len]] }
        else { $out{$name} = [[$s + 1, $L], [1, $len - ($L - $s)]]; push @note, "$name crosses the record's origin and is written as two rows" } }
    push @note, sprintf("the IR copies differ in length (%d and %d nt)", $b - $a, $d - $c) if $b - $a != $d - $c;
    push @note, "the record does not begin with the LSC" unless @{ $out{LSC} } && $out{LSC}[0][0] == 1;
    $out{note} = join('; ', @note);
    return \%out;
}

1;
