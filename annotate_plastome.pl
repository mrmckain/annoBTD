#!/usr/bin/perl
use strict;
use warnings;

#Plastome annotation - ORF discovery
#USAGE: 1: Full sequence 2: species [3: max overlap fraction, default 0.9]
#       [4: min ORF nt, default 19] [5: min length ratio vs coverer, default 0.8]
#
#Emits <species>_orffinder_coordinates.txt and <species>_orffinder_seqs.fsa
#Coordinates are 0-based, inclusive, on the forward strand.

my $MAX_OVERLAP = defined $ARGV[2] ? $ARGV[2] : 0.9;  #drop an ORF if this fraction of it is already claimed
my $MIN_ORF_NT  = defined $ARGV[3] ? $ARGV[3] : 19;   #keep ORFs strictly longer than this many nt
#...but never drop an ORF that is nearly as long as whatever already covers it.
#Overlapping ORFs on one strand are in different frames, so a covered ORF can be a
#different real gene rather than a redundant fragment. psbJ and psaI were being
#lost exactly this way: the ORF matching the gene's true coordinates was dropped
#in favour of a slightly longer neighbour in another frame.
my $MIN_LEN_RATIO = defined $ARGV[4] ? $ARGV[4] : 0.6;

die "usage: $0 <plastome.fsa> <species> [max_overlap] [min_orf_nt]\n" unless defined $ARGV[1];

my $plastome = "";
open my $cp_file, "<", $ARGV[0] or die "Cannot open $ARGV[0]: $!";
my $nrec = 0;
while(<$cp_file>){
	chomp;
	if(/^>/){
		$nrec++;
		next;
	}
	s/\s+//g;
	$plastome .= $_;
}
close $cp_file;
die "No sequence found in $ARGV[0]\n" unless length $plastome;
warn "WARNING: $ARGV[0] holds $nrec records; they were concatenated into one sequence.\n" if $nrec > 1;
$plastome = uc $plastome;

#Identify ORFs in sequence in both directions
my %forward_orfs = &orf_finder($plastome, $MIN_ORF_NT);
my $rc_plastome = reverse($plastome);
$rc_plastome =~ tr/ATCGatcg/TAGCtagc/;
my %reverse_orfs = &orf_finder($rc_plastome, $MIN_ORF_NT);

my $pl_len = length($plastome);

#Combine ORFs into single forward-based coordinates.
#Strand is tracked separately: an ORF on one strand must not mask an ORF on the other.
my %orf_strand;
$orf_strand{$_} = "+" for keys %forward_orfs;

for my $rid (sort keys %reverse_orfs){
	$reverse_orfs{$rid} =~ /(\d+)-(\d+)/;
	my $rstart = $1;
	my $rend = $2;
	my $nstart = $pl_len - $rend - 1;
	my $nend = $pl_len - $rstart - 1;
	my $orfid = $rid . "-r";
	$forward_orfs{$orfid} = $nstart . "-" . $nend;
	$orf_strand{$orfid} = "-";
}

#Filter ORFs to remove nested/redundant ORFs, largest first, independently per strand.
my %orf_size;
for my $orfid (sort keys %forward_orfs){
	my @tarray = split("-", $forward_orfs{$orfid});
	$orf_size{$orfid} = abs($tarray[1] - $tarray[0]);
}

#Array-backed coverage tracks are far cheaper than a per-base hash. Each claimed
#position stores the length of the ORF that claimed it, so a candidate can be
#compared against its coverer rather than just asking "is this covered".
my %covered_track = ("+" => [], "-" => []);

my %good_orfs;
my @sorted_size_orfs = sort { $orf_size{$b} <=> $orf_size{$a} || $a cmp $b } keys %orf_size;
for my $sorted_orfs (@sorted_size_orfs){
	my ($s, $e) = split("-", $forward_orfs{$sorted_orfs});
	my $track = $covered_track{ $orf_strand{$sorted_orfs} };
	my $span = $e - $s + 1;
	next if $span <= 0;

	my $covered = 0;
	my $max_claim = 0;
	for (my $k = $s; $k <= $e; $k++){
		next unless $track->[$k];
		$covered++;
		$max_claim = $track->[$k] if $track->[$k] > $max_claim;
	}
	#Previously compared against 2, which a fraction can never reach, so nothing was ever filtered.
	if(($covered / $span) > $MAX_OVERLAP){
		#Redundant only if it is also substantially shorter than its coverer.
		next if !$max_claim || ($span / $max_claim) < $MIN_LEN_RATIO;
	}

	$track->[$_] = $span for $s .. $e;
	$good_orfs{$sorted_orfs} = $forward_orfs{$sorted_orfs};
}

open my $outfile, ">", $ARGV[1] . "_orffinder_coordinates.txt" or die "Cannot write coordinates: $!";
open my $seqfile, ">", $ARGV[1] . "_orffinder_seqs.fsa" or die "Cannot write seqs: $!";
for my $ssiv (sort keys %good_orfs){
	$good_orfs{$ssiv} =~ /(\d+)-(\d+)/;
	my $start = $1;
	my $end = $2;
	my $tseq = substr($plastome, $start, ($end - $start + 1));
	print $outfile "$ssiv\t$forward_orfs{$ssiv}\n";
	print $seqfile ">$ssiv\n$tseq\n";
}
close $outfile;
close $seqfile;

printf STDERR "%s: %d ORFs found, %d kept (max_overlap=%s, min_orf_nt=%d, min_len_ratio=%s)\n",
	$ARGV[1], scalar(keys %forward_orfs), scalar(keys %good_orfs), $MAX_OVERLAP, $MIN_ORF_NT, $MIN_LEN_RATIO;


sub orf_finder
{
	my ($cur_plastome, $min_nt) = @_;

	my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');

	my %orfs;
	my $orf_count = 0;
	my $len = length($cur_plastome);

	for (my $i = 0; $i <= 2; $i++){
		my ($cur_start, $cur_end);
		my $orf_len = 0;

		for (my $j = $i; $j <= $len - 3; $j += 3){
			my $curcodon = substr($cur_plastome, $j, 3);
			my $aa = $codons{$curcodon};

			if(!defined $aa){
				#Ambiguity or non-ACGT: close the current ORF.
				if($orf_len > $min_nt){
					$orf_count++;
					$orfs{"orf$orf_count"} = $cur_start . "-" . $cur_end;
				}
				($cur_start, $cur_end, $orf_len) = (undef, undef, 0);
				next;
			}

			if(defined $cur_start){
				$orf_len += 3;
				$cur_end = $j + 2;
				if($aa eq "_"){
					if($orf_len > $min_nt){
						$orf_count++;
						$orfs{"orf$orf_count"} = $cur_start . "-" . $cur_end;
					}
					($cur_start, $cur_end, $orf_len) = (undef, undef, 0);
				}
			}
			elsif($aa ne "_"){
				$cur_start = $j;
				$cur_end   = $j + 2;
				$orf_len   = 3;
			}
		}

		#Flush the ORF still open at the end of this frame.
		if(defined $cur_start && $orf_len > $min_nt){
			$orf_count++;
			$orfs{"orf$orf_count"} = $cur_start . "-" . $cur_end;
		}
	}
	return %orfs;
}
