#!/usr/bin/perl
use strict;
use warnings;

# Rank candidate guide genomes by k-mer similarity to the genome being annotated.
#
# USAGE: rank_references_by_kmer.pl <target.fsa> <guide1.fsa> [guide2.fsa ...] > ranking.txt
#    or: rank_references_by_kmer.pl <target.fsa> --dir <directory of guide .fsa files>
#
# OUT: guide_name <TAB> jaccard <TAB> containment <TAB> rank      (closest first)
#
# Why this instead of picking guides taxonomically: exon structure is lineage
# specific. Grasses carry different exon counts from tobacco for several genes, so
# a guide chosen for convenience rather than for relatedness transfers the wrong
# splicing pattern - and annoBTD copies whatever the winning guide says. Ranking by
# k-mer similarity is a cheap proxy for relatedness that needs no taxonomy lookup
# and degrades gracefully for anything unplaced.
#
# Both strands are canonicalised (each k-mer is stored as the lexicographically
# smaller of itself and its reverse complement) so assembly orientation does not
# matter. Similarity is reported two ways:
#   jaccard     |A n B| / |A u B|  - symmetric, penalises genome size differences
#   containment |A n B| / |A|      - fraction of the TARGET covered by the guide,
#                                    which is the more useful number here because a
#                                    guide may be much larger or smaller than the
#                                    target and still annotate it well.
# Ranking is by containment, with jaccard as the tie-break.

my $K = $ENV{ANNOBTD_KMER_K} || 21;

my $target = shift @ARGV or die
	"usage: $0 <target.fsa> <guide.fsa ...>\n" .
	"       $0 <target.fsa> --dir <guide_directory>\n";

my @guides;
if (@ARGV && $ARGV[0] eq '--dir') {
	shift @ARGV;
	my $dir = shift @ARGV or die "--dir needs a directory\n";
	opendir(my $dh, $dir) or die "Cannot read $dir: $!";
	@guides = map { "$dir/$_" }
	          grep { /\.(fsa|fa|fasta|fna)$/i || -f "$dir/$_" }
	          grep { !/^\./ } readdir($dh);
	closedir $dh;
	@guides = grep { -f $_ } @guides;
}
else {
	@guides = @ARGV;
}
die "No guide genomes given\n" unless @guides;

sub read_seq {
	my $file = shift;
	open my $fh, "<", $file or die "Cannot open $file: $!";
	my $s = "";
	while (<$fh>) { chomp; next if /^>/; s/\s+//g; $s .= $_ }
	close $fh;
	return uc $s;
}

# Canonical k-mer set. Plastomes are ~150 kb, so the exact set is small enough to
# hold outright - no sketching needed, and the answer is exact.
sub kmer_set {
	my ($seq, $k) = @_;
	my %set;
	my $n = length($seq) - $k;
	for my $i (0 .. $n) {
		my $kmer = substr($seq, $i, $k);
		next if $kmer =~ /[^ACGT]/;          # skip ambiguity
		my $rc = reverse $kmer;
		$rc =~ tr/ACGT/TGCA/;
		$set{ $kmer le $rc ? $kmer : $rc } = 1;
	}
	return \%set;
}

my $tseq = read_seq($target);
die "No sequence in $target\n" unless length $tseq;
my $tset = kmer_set($tseq, $K);
my $tn = scalar keys %$tset;
die "Target $target yielded no usable ${K}-mers\n" unless $tn;

my @rows;
for my $g (@guides) {
	my $gseq = read_seq($g);
	next unless length $gseq;
	my $gset = kmer_set($gseq, $K);
	my $gn = scalar keys %$gset;
	next unless $gn;

	my $shared = 0;
	# Iterate the smaller set.
	my ($small, $large) = $tn <= $gn ? ($tset, $gset) : ($gset, $tset);
	for my $kmer (keys %$small) { $shared++ if exists $large->{$kmer} }

	my $union = $tn + $gn - $shared;
	my $name = $g; $name =~ s{.*/}{}; $name =~ s{\.(fsa|fa|fasta|fna)$}{}i;
	push @rows, {
		name        => $name,
		jaccard     => $union ? $shared / $union : 0,
		containment => $shared / $tn,
	};
}
die "No guide genome produced usable k-mers\n" unless @rows;

@rows = sort { $b->{containment} <=> $a->{containment}
            || $b->{jaccard}     <=> $a->{jaccard}
            || $a->{name} cmp $b->{name} } @rows;

my $rank = 0;
for my $r (@rows) {
	$rank++;
	printf "%s\t%.6f\t%.6f\t%d\n", $r->{name}, $r->{jaccard}, $r->{containment}, $rank;
}

printf STDERR "%s: ranked %d guides by %d-mer similarity (closest: %s, containment %.3f)\n",
	$target, scalar(@rows), $K, $rows[0]{name}, $rows[0]{containment};
