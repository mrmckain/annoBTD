#!/usr/bin/perl -w
use strict;

my %blastgenes;
my %blastrnas;
my %orfs;
my %idorfs_f;
my %short_exons;
my %short_last_exons;
#Profile expectation (third scoring signal) - loaded HERE, at startup. It was first
#placed beside expected_nt() near the bottom of the file, and since the main body
#runs before those lines execute, every call saw an empty table.
my (%EXPECT);
if ($ENV{ANNOBTD_PROFILE} && -s $ENV{ANNOBTD_PROFILE}) {
	open my $ph, "<", $ENV{ANNOBTD_PROFILE} or die;
	while (<$ph>) { chomp; my ($n, $lt, $lin, $cnt, $med, $cons, $alt) = split /\t/; next if $n eq "name"; $EXPECT{$n}{"$lt:$lin"} = [$med, $cnt, $cons // 1, $alt || 0] }
	close $ph;
}
my $EXP_FAM = $ENV{ANNOBTD_FAMILY} // ''; my $EXP_ORD = $ENV{ANNOBTD_ORDER} // ''; my $EXP_GEN = $ENV{ANNOBTD_GENUS} // '';
#Genus stratum ("G:" cells, one vote per species): consulted first when the target's genus
#has at least this many species in the cell. A family split between two conventions is
#often settled inside one genus.
my $EXP_MIN_N_GENUS = defined $ENV{ANNOBTD_EXPECT_MIN_N_GENUS} ? $ENV{ANNOBTD_EXPECT_MIN_N_GENUS} : 10;
#Validated on 17 genomes against curated truth: where >=95% of a family's records
#agree on a length the expectation matched truth in 98.6% of CDS; at 0.90-0.95,
#94%; below 0.70 it was a coin flip (49.6%). So an expectation is used only when
#the cell is SETTLED, and never borrowed from ALL lineages - the top-band misses
#were almost all Marchantia, a liverwort scored against an angiosperm-wide
#consensus because no Marchantiaceae cell existed yet.
#MIN_N counts VOTES in the cell. Since gene_expect6 a vote is a species (one vote per
#independent submitter group, then one per species: build_species_consensus.pl), so
#ten species is a stronger basis than the twenty records the gate first meant.
my $EXP_MIN_N = defined $ENV{ANNOBTD_EXPECT_MIN_N} ? $ENV{ANNOBTD_EXPECT_MIN_N} : 10;
my $EXP_MIN_CONS = defined $ENV{ANNOBTD_EXPECT_MIN_CONSENSUS} ? $ENV{ANNOBTD_EXPECT_MIN_CONSENSUS} : 0.90;
#Weak mode: a populated cell whose majority is below MIN_CONS but at least this
#may still drive an adjustment when the CURRENT call matches neither of its
#modes (Pinus petA: 86% at 960, 14% at 996, call at 690). Choosing between the
#modes is left alone; leaving a length no record supports is not. Measured on the
#17: +1 (Pinus petA) / -3 (Poaceae ndhK, whose true 684 form has ZERO share in
#GenBank - the convention there is 747/741/738). Off by default (threshold > 1).
my $EXP_WEAK_CONS = defined $ENV{ANNOBTD_EXPECT_WEAK_CONSENSUS} ? $ENV{ANNOBTD_EXPECT_WEAK_CONSENSUS} : 1.01;

my $sid;
my $torf_counter=1;
my %orf_pos;
my ($idorfs_f, $orf_pos);

open my $file, "<", $ARGV[3]; #orf positional file
while(<$file>){
		chomp;
		my @tarray = split /\s+/;
		$tarray[1] =~ /(\d+)-(\d+)/;
		$orf_pos{$tarray[0]}{"Start"}=$1;
		$orf_pos{$tarray[0]}{"End"}=$2;

}
close $file;

%blastgenes = &read_fasta($ARGV[0],%blastgenes);
%orfs = &read_fasta($ARGV[1],%orfs);
%blastrnas = &read_fasta($ARGV[5],%blastrnas);
%blastrnas = &read_fasta($ARGV[6],%blastrnas);

#Read the plastome before the RNA blasts: rna_blasts needs its length to clamp
#the coordinates it extends.
my $plastome;
open $file, "<", $ARGV[4] or die "Cannot open plastome $ARGV[4]: $!"; #full plastome being annotated
while(<$file>){
	chomp;
	if(/^>/){
			next;
	}
	else{
		$plastome.=$_;
	}
}
close $file;
die "No sequence read from $ARGV[4]\n" unless defined $plastome && length $plastome;
$plastome =~ tr/atcg/ATCG/;   #v2.6

($torf_counter,$idorfs_f, $orf_pos) = &rna_blasts($torf_counter, $ARGV[7], \%idorfs_f, \%orf_pos, $ARGV[10], length($plastome));
%orf_pos = %$orf_pos;
%idorfs_f = %$idorfs_f;
($torf_counter,$idorfs_f,$orf_pos) = &rna_blasts($torf_counter, $ARGV[8], \%idorfs_f, \%orf_pos, $ARGV[11], length($plastome));
%orf_pos = %$orf_pos;
%idorfs_f = %$idorfs_f;

for my $gid (sort keys %blastgenes){
		if(length($blastgenes{$gid})<25){
			$short_exons{$gid}=1;
		}
		#25-34 nt is handled SEPARATELY, and deliberately not by adding it to
		#%short_exons. Membership there does two things at once: it enables the
		#exact-sequence search, and it routes the gene down an emit path that writes
		#raw ORF coordinates with no boundary refinement. rps12_exon3 (26 nt, 29 in
		#some guides) needs the first and must not have the second - tblastx does
		#reach it in the nine-genome set and refines it exactly, and bypassing that
		#pushed the boundary ~98 nt out in five genomes. So it gets its own list,
		#used only by the last-exon search below.
		if(length($blastgenes{$gid})>=25 && length($blastgenes{$gid})<35
		   && $gid =~ /_exon(\d+)XXX/ && $1 > 1){
			$short_last_exons{$gid}=1;
		}
}

for my $gid (sort keys %blastrnas){
		if(length($blastrnas{$gid})<26){
			$short_exons{$gid}=1;
		}
}

#Guide genomes disagree about exon boundaries, and the disagreement is often
#lineage specific rather than error: grasses carry different exon counts from
#tobacco for several genes. Whichever single guide wins the annotation score has
#its structure transferred verbatim, so picking a guide from the wrong lineage
#transfers the wrong splicing pattern.
#
#An earlier version of this took the majority exon length across all guides. That
#is exactly wrong when the target's own lineage is the minority in the guide set -
#annotating a grass against mostly eudicot guides would vote the grass structure
#away. Instead, prefer the reference from the most closely related guide, ranked
#by k-mer similarity to the genome being annotated (see rank_references_by_kmer.pl).
#
#ANNOBTD_REF_RANKING names a file of "guide_index <TAB> guide_name <TAB> rank",
#rank 1 being the closest. Without it, the score-chosen reference is used as-is.
my %guide_rank;
if($ENV{ANNOBTD_REF_RANKING} && -s $ENV{ANNOBTD_REF_RANKING}){
	open my $rh, "<", $ENV{ANNOBTD_REF_RANKING} or die "Cannot open ranking: $!";
	while(<$rh>){
		chomp;
		my @f = split /\t/;
		next unless defined $f[2] && $f[0] =~ /^\d+$/;
		$guide_rank{$f[0]} = $f[2];
	}
	close $rh;
	print STDERR "reference ranking: " . scalar(keys %guide_rank) . " guides ranked\n";
}

my %ref_by_gene;
for my $rid (sort keys %blastgenes){
	next unless $rid =~ /^(.+)XXX(.+)$/;
	push @{$ref_by_gene{$1}}, $rid;
}

#Swap the score-chosen reference for the same gene from a more closely related
#guide, when one exists. The ORF and the reported gene name are untouched; only
#the sequence used for boundary refinement changes.
sub closest_reference {
	my $rid = shift;
	return $rid unless %guide_rank;
	return $rid unless defined $rid && $rid =~ /^(.+)XXX(.+)$/;
	my ($base, $idx) = ($1, $2);
	$idx =~ s/_\d+$//;   #strip the duplicate-name disambiguator added upstream
	my $sibs = $ref_by_gene{$base} or return $rid;
	return $rid unless @$sibs > 1;

	my $best = $rid;
	my $best_rank = defined $guide_rank{$idx} ? $guide_rank{$idx} : 1e9;
	for my $s (sort @$sibs){
		next unless $s =~ /XXX(.+)$/;
		my $r = defined $guide_rank{$1} ? $guide_rank{$1} : 1e9;
		if($r < $best_rank){ ($best, $best_rank) = ($s, $r) }
	}
	return $best;
}

my %blast_overlaps;
my %best_orf_scores_id;   #v2.6: annotation score per (reference, orf)
open $file, "<", $ARGV[2] or die "Cannot open $ARGV[2]: $!"; #best id to orf
while(<$file>){
		chomp;
		if(/intron/){
			next;
		}
		my @tarray = split /\s+/;
		next unless defined $tarray[1];
		$tarray[1] = &closest_reference($tarray[1]);
		$tarray[1] =~ /(.*?)XXX.+/;
		my $temp_geneid = $1;

		$idorfs_f{$tarray[1]}{$tarray[0]}{"Start"}=$orf_pos{$tarray[0]}{"Start"};
		$idorfs_f{$tarray[1]}{$tarray[0]}{"End"}=$orf_pos{$tarray[0]}{"End"};		
		$blast_overlaps{$tarray[1]}{$tarray[0]} = abs($orf_pos{$tarray[0]}{"End"}-$orf_pos{$tarray[0]}{"Start"});
		#v2.6: keep the score too - resolving overlaps by which ORF matched the
		#reference better beats resolving them by which ORF is longer.
		$best_orf_scores_id{$tarray[1]}{$tarray[0]} = defined $tarray[2] ? $tarray[2] : 0;
		
}
close $file;

print STDERR "TRACE A: before hit_cleaner\n" if $ENV{ANNOBTD_TRACE};
%idorfs_f = &hit_cleaner(%idorfs_f); #cleans up positions of genes to reconcile with overlaps


print STDERR "TRACE B: after hit_cleaner, before short-exon loop\n" if $ENV{ANNOBTD_TRACE};
my $orf_counter=0;
for my $shid (sort keys %short_exons, sort keys %short_last_exons){
		print STDERR "TRACE   short_exon $shid\n" if $ENV{ANNOBTD_TRACE};

		my $tseq;
		my ($gene, $exonnum, $refid);
		if($shid =~ /trn/){
			$tseq=$blastrnas{$shid};
			$shid =~ /(trn\w+-\w\w\w)_exon(\d+)(XXX.+)/;
			$gene=$1;
			$exonnum=$2;
			$refid=$3;
		}
		else{
			$tseq = $blastgenes{$shid};
			$shid =~ /(\w+)_exon(\d+)(XXX.+)/;
			$gene = $1;
			$exonnum = $2;
			$refid = $3;
		}
		
		
		
		#A short reference whose name does not carry the _exonN pattern leaves these
		#undefined - trnG_exon1XXX1 parses, a bare short gene does not. Skip rather
		#than comparing undef.
		next unless defined $exonnum && defined $gene && defined $refid;
		if($exonnum == 1){
				my $intron_name = $gene . "_intron1";

				#if(exists $blastgenes{$intron_name}){
					$exonnum++;
					my $next_orf = $gene . "_exon" . $exonnum . $refid;
					
					if(scalar keys %{$idorfs_f{$next_orf}} == 1){
							for my $norf (sort keys %{$idorfs_f{$next_orf}}){
								if($norf =~ /\-r/){ #taking care of the reverse complement first
									my @tarray;
									#my $add =substr($blastgenes{$intron_name},-3); #adding some intron seq
									#$tseq = $add . $tseq;
									$tseq = reverse($tseq);
									$tseq =~ tr/ATCGatcg/TAGCtagc/;
									my $j=0;
									for (my $i = 0; $i < length($plastome);$i=$j){
										if(index($plastome,$tseq,$i) > 0){
											push(@tarray,index($plastome,$tseq,$i));
											$j=index($plastome,$tseq,$i)+1;
											
										}
										else{
												$j=length($plastome);
										}

									}
									
									my $close_end = $orf_pos{$norf}{"End"};
									print STDERR "TRACE     rev candidates=".scalar(@tarray)."\n" if $ENV{ANNOBTD_TRACE};
								my $good_start = &pick_short_exon(\@tarray, length($tseq), $close_end, "-");
									
									#$good_start+=5;
									my $elen;
									if($shid =~ /trn/){
										$elen = length($blastrnas{$shid});
									}
									else{
										$elen = length($blastgenes{$shid});
									}
									$orf_counter++;
									my $new_orf = "orfx" . $orf_counter. "-r";
									$idorfs_f{$shid}{$new_orf}{"Start"}=$good_start+1;
									$idorfs_f{$shid}{$new_orf}{"End"}=$good_start+$elen-1+1;

								}
								else{
									my @tarray;
									#my $add =substr($blastgenes{$intron_name},0,3); #adding some intron seq
									#$tseq = $add . $tseq;
									my $j=0;
									for (my $i = 0; $i < length($plastome);$i=$j){
										if(index($plastome,$tseq,$i) > 0){
											push(@tarray,index($plastome,$tseq,$i));
											$j=index($plastome,$tseq,$i)+1;
											
										}
										else{
												$j=length($plastome);
										}

									}
									
									my $close_start = $orf_pos{$norf}{"Start"};
									print STDERR "TRACE     fwd candidates=".scalar(@tarray)."\n" if $ENV{ANNOBTD_TRACE};
								my $good_start = &pick_short_exon(\@tarray, length($tseq), $close_start, "+");
									
									my $elen;
									if($shid =~ /trn/){
										$elen = length($blastrnas{$shid});
									}
									else{
										$elen = length($blastgenes{$shid});
									}
									$orf_counter++;
									my $new_orf = "orfx" . $orf_counter;
									$idorfs_f{$shid}{$new_orf}{"Start"}=$good_start+1;
									$idorfs_f{$shid}{$new_orf}{"End"}=$good_start+$elen-1+1;

								}
							}
					}
				#}
		}
		elsif($exonnum > 1 && $short_last_exons{$shid}){
			#A short LAST exon, anchored on the exon before it. rps12_exon3 is 26 nt and
			#tblastx never reports it, but its exon2 is placed correctly, so the exon can
			#be found by exact sequence search relative to that anchor - the same trick the
			#exon1 branch above uses in mirror image. 9 of 15 rps12_exon3 guide references
			#match this genome exactly, so the search has plenty to work with.
			#
			#Unlike the exon1 branch this does NOT require the anchor exon to be unique:
			#rps12_exon2 is IR-duplicated and both copies need their own exon3.
			#Skip when some guide has already located this exon. tblastx DOES reach
			#rps12_exon3 in the nine-genome set and refines it exactly; a synthetic exon
			#on top of that competed with the correct call. The Asparagaceae need this
			#branch only because tblastx reports nothing for their copy at all.
			my $already = 0;
			for my $k (keys %idorfs_f){
				next unless $k =~ /^\Q$gene\E_exon${exonnum}XXX/;
				if(scalar keys %{$idorfs_f{$k}}){ $already = 1; last }
			}
			next if $already;
			my $prev_orf = $gene . "_exon" . ($exonnum - 1) . $refid;
			next unless exists $idorfs_f{$prev_orf};
			for my $norf (sort keys %{$idorfs_f{$prev_orf}}){
				my $dir = ($norf =~ /\-r/) ? "-" : "+";
				my $probe = $tseq;
				if($dir eq "-"){ $probe = reverse($probe); $probe =~ tr/ATCGatcg/TAGCtagc/ }
				next unless defined $probe && length($probe) >= 12;
				my @cands;
				my $j = 0;
				for (my $i = 0; $i < length($plastome); $i = $j){
					my $hit = index($plastome, $probe, $i);
					if($hit > 0){ push @cands, $hit; $j = $hit + 1 }
					else { $j = length($plastome) }
				}
				next unless @cands;
				my $anchor = ($dir eq "-") ? $orf_pos{$norf}{"Start"} : $orf_pos{$norf}{"End"};
				next unless defined $anchor;
				print STDERR "TRACE     last-exon $shid candidates=".scalar(@cands)."\n" if $ENV{ANNOBTD_TRACE};
				my $good_start = &pick_last_exon(\@cands, length($probe), $anchor, $dir);
				next unless defined $good_start;
				$orf_counter++;
				my $new_orf = "orfx" . $orf_counter . ($dir eq "-" ? "-r" : "");
				$idorfs_f{$shid}{$new_orf}{"Start"} = $good_start + 1;
				$idorfs_f{$shid}{$new_orf}{"End"}   = $good_start + length($probe);
			}
		}
}
print STDERR "TRACE C: after short-exon loop\n" if $ENV{ANNOBTD_TRACE};
my %gene_features;
for my $gid (sort keys %idorfs_f){
	if($gid =~ /exon/){
			$gid =~ /(.+)_exon(\d+)/;
			my $full_name = $1. "_exon" . $2;
			$gene_features{$1}{$full_name}{$gid}=1;
	}
	else{
			
			$gene_features{$gid}{$gid}=1;
	}
}
#Identify start/stop for genes (not tRNA, rRNA, and short exons)
#Diagnostic hook: set ANNOBTD_DEBUG_COORDS=<file> to dump the intermediate values
#behind every emitted CDS coordinate, for reconciling against known-good annotation.
my $DEBUG_COORDS;
if($ENV{ANNOBTD_DEBUG_COORDS}){
	open($DEBUG_COORDS, ">", $ENV{ANNOBTD_DEBUG_COORDS})
		or die "Cannot open coord debug file: $!";
	print $DEBUG_COORDS join("\t", qw(sub gid eid branch orf_start orf_end frame
		frame_restore mod_shift nearstart nearend start_match start_alt end_match
		emit_start emit_end dir)), "\n";
}
sub dbg_coord {
	return unless $DEBUG_COORDS;
	print $DEBUG_COORDS join("\t", map { defined $_ ? $_ : "" } @_), "\n";
}
my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');
my %final_annotation;
for my $gid (sort keys %idorfs_f){
	print "$gid\n";
	if($ENV{ANNOBTD_TRACE}){ $| = 1; print STDERR "TRACE   gene $gid\n" }
	for my $eid (sort keys %{$idorfs_f{$gid}}){
		#$eid =~ /^orfx/ covers the synthetic exons built by exact sequence match
			#below. Those coordinates are exact by construction, so they take the direct
			#emit path like any other short exon - rps12_exon3 is not in %short_exons
			#(see the note there) and would otherwise be sent through boundary
			#refinement, which dropped one of its two IR copies.
		if(exists $short_exons{$gid} || $eid =~ /^orfx/ || $gid =~ /rrn/ || $gid =~ /trn/){
			my $dir = "+";
			if($eid =~ /\-r/){
				$dir = "-";
			}
			if($eid =~ /rrn45/){
				$eid =~ /rrn45(XXX.+)/;
				$eid = "rrn4.5" . $1;
			}
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}}{$idorfs_f{$gid}{$eid}{"End"}}{$dir}=$gid;
		}
		else{
			my $blastsize = length($blastgenes{$gid});
			my $blmod = $blastsize%3;
			if($gid =~ /exon/){
				$gid =~ /(.+)_exon(\d+)/;
				my $full_name = $1. "_exon" . $2;
				my $totalexons = scalar keys %{$gene_features{$1}};
				$gid =~ /exon(\d+)/;
				my $cur_exon = $1;
				if(($cur_exon-1) == 0){
					&exon_mods_start($blmod,$gid,$eid);
				}
				elsif(($totalexons - $cur_exon) == 0){
					&exon_mods_end($blmod,$gid,$eid);
					#&exon_mods_start($blmod,$gid,$eid);
				}
				else{
					#&exon_mods_start($blmod,$gid,$eid);
					&exon_mods_mid($blmod,$gid,$eid);					
				}
			}

			else{
				&exon_mods_start(0,$gid,$eid);
			}

		}
	}
}

# Splice-site refinement is opt-in: see the note on the sub. It moves correct
# boundaries when a genuine intron lacks a GT donor, which is common.
print STDERR "TRACE D: after main gene loop\n" if $ENV{ANNOBTD_TRACE};
&refine_intron_boundaries() if $ENV{ANNOBTD_SPLICE_REFINE};

my %temp_store;
for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (sort {$a <=> $b} keys %{$final_annotation{$start}}){
		if(%temp_store){
		for my $temp_start (sort {$a <=> $b} keys %temp_store){
			for my $temp_end (sort {$a <=> $b} keys %{$temp_store{$temp_start}}){
				if($start == $temp_start && $end == $temp_end){
					for my $dir (sort keys %{$final_annotation{$start}{$end}}){
							if($final_annotation{$start}{$end}{$dir} eq $temp_store{$temp_start}{$temp_end}){
								next;
							}
					}
				}
				elsif($start > $temp_start && $end <= $temp_end){
					delete $final_annotation{$start}{$end};

				}
				elsif($temp_start > $start && $temp_end <= $end){
					delete $final_annotation{$temp_start}{$temp_end};
					for my $dir (sort keys %{$final_annotation{$start}{$end}}){
						my $next_gene = $final_annotation{$start}{$end}{$dir};
						$next_gene =~ /(.*?)XXX.+/;
						$next_gene = $1;
						$temp_store{$start}{$end}=$next_gene;
					}
				}
				elsif($start >= $temp_start && $end < $temp_end){
					delete $final_annotation{$start}{$end};
				}
				elsif($temp_start >= $start && $temp_end < $end){
					delete $final_annotation{$temp_start}{$temp_end};
					for my $dir (sort keys %{$final_annotation{$start}{$end}}){
						my $next_gene = $final_annotation{$start}{$end}{$dir};
						$next_gene =~ /(.*?)XXX.+/;
						$next_gene = $1;
						$temp_store{$start}{$end}=$next_gene;
					}
				}
				elsif($start < $temp_start && $temp_end > $end){
					if(abs($end-$start) < abs($temp_end-$temp_start)){
						delete $final_annotation{$start}{$end};
					}
					if(abs($temp_end-$temp_start) < abs($end-$start)){
						delete $final_annotation{$temp_start}{$temp_end};
						for my $dir (sort keys %{$final_annotation{$start}{$end}}){
							my $next_gene = $final_annotation{$start}{$end}{$dir};
							$next_gene =~ /(.*?)XXX.+/;
							$next_gene = $1;
							$temp_store{$start}{$end}=$next_gene;
						}
					}
				} 
				elsif($temp_start > $start && $temp_end > $end){
					if(abs($end-$start) < abs($temp_end-$temp_start)){
						delete $final_annotation{$start}{$end};
					}
					if(abs($temp_end-$temp_start) < abs($end-$start)){
						delete $final_annotation{$temp_start}{$temp_end};
						for my $dir (sort keys %{$final_annotation{$start}{$end}}){
							my $next_gene = $final_annotation{$start}{$end}{$dir};
							$next_gene =~ /(.*?)XXX.+/;
							$next_gene = $1;
							$temp_store{$start}{$end}=$next_gene;
						}
					}
				}
				else{
					for my $dir (sort keys %{$final_annotation{$start}{$end}}){
						my $next_gene = $final_annotation{$start}{$end}{$dir};
						$next_gene =~ /(.*?)XXX.+/;
						$next_gene = $1;
							$temp_store{$start}{$end}=$next_gene;
						}
				}
			}
		}
		}
		else{
			for my $dir (sort keys %{$final_annotation{$start}{$end}}){
				my $next_gene = $final_annotation{$start}{$end}{$dir};
				$next_gene =~ /(.*?)XXX.+/;
				$next_gene = $1;
				$temp_store{$start}{$end}=$next_gene;
			}	
		}
	}
}


my $prev;
my $prev_end;
for my $start (sort {$a <=> $b} keys %final_annotation){
	
	for my $end (sort {$a <=> $b} keys %{$final_annotation{$start}}){
		
		for my $dir (sort keys %{$final_annotation{$start}{$end}}){
			my $tids = $final_annotation{$start}{$end}{$dir};
			$tids =~ /(.*?)XXX.+/;
			$tids = $1;
			if(!$prev){
					$final_annotation{1}{$start-1}{"+"}="Start\~" . $tids;
					$prev = $tids;
					$prev_end=$end;
			}
			else{
				if($start-$prev_end <=0){
					$prev = $tids;
					$prev_end=$end;
				}
				elsif ($prev =~ /exon/ && $tids =~ /exon/){
					$prev =~ /(.+)_exon/;
					my $pgene = $1;
					$tids =~ /(.+)_exon/;
					my $tgene = $1;
					if($tgene eq $pgene){
						$prev =~ /_exon(\d)/;
						my $pnum = $1;
						$tids =~ /_exon(\d)/;
						my $tnum = $1;
						if ($pnum == 3 && $tnum == 2){
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron2";
							$prev = $tids;
							$prev_end = $end;
						}
						elsif ($pnum ==2 && $tnum ==1) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
							$prev = $tids;
							$prev_end = $end;
						}
						elsif ($pnum ==1 && $tnum ==2) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
							$prev = $tids;
							$prev_end = $end;
						}
						elsif ($pnum ==2 && $tnum ==3) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron2";
							$prev = $tids;
							$prev_end = $end;
						}
					}
					else{
						my $next_gene = $final_annotation{$start}{$end}{$dir};
						$next_gene =~ /(.*?)XXX.+/;
						$next_gene = $1;
						$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "~" . $next_gene;
						$prev = $next_gene;
						$prev_end = $end;
					}
				}
				else{
					my $next_gene = $final_annotation{$start}{$end}{$dir};
					$next_gene =~ /(.*?)XXX.+/;
					$next_gene = $1;
					$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "~" . $next_gene;
					$prev = $next_gene;
					$prev_end = $end;
				}
			}
		}
	}
}

$final_annotation{$prev_end+1}{length($plastome)}{"+"}=$prev . "~End";




my $lsc;
my $ssc;
my $irb;
my $ira;

my $boundary1= substr(reverse($plastome), 0, 20);
$boundary1 =~ tr/ATCGatcg/TAGCtagc/;
my $start_irb = index($plastome, $boundary1);

$lsc = substr($plastome, 0, $start_irb);

my $i = $start_irb;

my $irb_end = substr($plastome, $i, 20);
$irb_end = reverse($irb_end);
$irb_end =~ tr/ATCGatcg/TAGCtagc/;

until($plastome !~ /$irb_end/){
	$i++;
	$irb_end = substr($plastome, $i, 20);
	$irb_end = reverse($irb_end);
	$irb_end =~ tr/ATCGatcg/TAGCtagc/;
}
$i--;
$irb_end = substr($plastome, $i, 20);
$irb = substr($plastome, $start_irb, ($i+20-$start_irb));
my $ira_seq = $irb_end;
$ira_seq = reverse($ira_seq);
$ira_seq =~ tr/ATCGatcg/TAGCtagc/;
my $ira_start = index($plastome, $ira_seq);
$ssc = substr($plastome, $i+20, $ira_start-($i+20));
$ira= substr($plastome, $ira_start);


my $lsc_range = (index($plastome, $lsc)+1) . "-" . (index($plastome, $lsc)+length($lsc));
my $ssc_range = (index($plastome, $ssc)+1) . "-" . (index($plastome, $ssc)+length($ssc));
my $ira_range = (index($plastome, $ira)+1) . "-" . (index($plastome, $ira)+length($ira));
my $irb_range = (index($plastome, $irb)+1) . "-" . (index($plastome, $irb)+length($irb));
my $full_range = 1 . "-" . length($plastome);



###########################
sub read_fasta{
my ($file, %temphash) = @_;
open my $tfile, "<", $file; #gene database fasta file
while(<$tfile>){
		chomp;
		if(/>/){
				$sid =substr($_, 1);
		}
		else{
				$temphash{$sid}.=$_;

		}
}
close $file;
return(%temphash);
}
###########################
sub hit_cleaner{
my %idorfs = @_;
for my $bid (sort keys %idorfs){
	
	if($bid =~ /intron/){
			next;
	}
	if(scalar keys %{$idorfs{$bid}}>1){
		my %temporfs;
		for my $oid (sort keys %{$idorfs{$bid}}){
			
				if(%temporfs){
						my $cur_start = $orf_pos{$oid}{"Start"};
						my $cur_end = $orf_pos{$oid}{"End"};
						for my $toid (sort keys %temporfs){
							if($temporfs{$toid}{"Start"} > $cur_start && $temporfs{$toid}{"End"} < $cur_end){
								if(($best_orf_scores_id{$bid}{$toid}||0) > ($best_orf_scores_id{$bid}{$oid}||0)){
									delete $idorfs{$bid}{$oid};
								}
								else{
									delete $temporfs{$toid};
									delete $idorfs{$bid}{$toid};
									$temporfs{$oid}{"Start"}=$cur_start;
									$temporfs{$oid}{"End"}=$cur_end;
								}
							}
							elsif($temporfs{$toid}{"Start"} < $cur_start && $temporfs{$toid}{"End"} > $cur_end){
								if(($best_orf_scores_id{$bid}{$toid}||0) > ($best_orf_scores_id{$bid}{$oid}||0)){
									delete $idorfs{$bid}{$oid};
								}
								else{
									delete $temporfs{$toid};
									delete $idorfs{$bid}{$toid};
									$temporfs{$oid}{"Start"}=$cur_start;
									$temporfs{$oid}{"End"}=$cur_end;
								}
							}
							elsif($temporfs{$toid}{"Start"} < $cur_start && $temporfs{$toid}{"End"} < $cur_end && $temporfs{$toid}{"End"} > $cur_start){
								if(($best_orf_scores_id{$bid}{$toid}||0) >= ($best_orf_scores_id{$bid}{$oid}||0)){
									delete $idorfs{$bid}{$oid};
								}
								else{
									delete $temporfs{$toid};
									delete $idorfs{$bid}{$toid};
									$temporfs{$oid}{"Start"}=$cur_start;
									$temporfs{$oid}{"End"}=$cur_end;
								}
							}
							elsif($temporfs{$toid}{"Start"} > $cur_start && $temporfs{$toid}{"End"} > $cur_end && $temporfs{$toid}{"Start"} < $cur_end){
								if(($best_orf_scores_id{$bid}{$toid}||0) > ($best_orf_scores_id{$bid}{$oid}||0)){
									delete $idorfs{$bid}{$oid};
								}
								else{
									delete $temporfs{$toid};
									delete $idorfs{$bid}{$toid};
									$temporfs{$oid}{"Start"}=$cur_start;
									$temporfs{$oid}{"End"}=$cur_end;
								}
							}
							else{
								$temporfs{$oid}{"Start"}=$cur_start;
								$temporfs{$oid}{"End"}=$cur_end;
							}
						}
				}
				else{
					my $cur_start = $orf_pos{$oid}{"Start"};
					my $cur_end = $orf_pos{$oid}{"End"};
					$temporfs{$oid}{"Start"}=$cur_start;
					$temporfs{$oid}{"End"}=$cur_end;
				}
				
				
		}

	}
}
return(%idorfs);
}
###########################
#Choose among candidate genomic positions for a very short exon.
#
#These exons are 6-9 nt, so their sequence occurs dozens of times by chance - the
#Spinacia petB exon 1 6-mer matches 42 places in that genome - and proximity to
#the neighbouring exon is not enough to pick between them: the wrong candidate
#there sat 47 nt CLOSER to exon 2 than the right one.
#
#The intron that follows a real short exon opens with GT in all 27 short-exon
#introns of the reference set (petB, petD, rpl16 across nine genomes), so prefer
#candidates carrying that donor. It is a preference, not a requirement: when no
#candidate has one, plain proximity still decides, so a gene whose donor differs
#is no worse off than before.
#
#  $cands  arrayref of 0-based candidate start positions for the short exon
#  $elen   exon length in nt
#  $anchor 0-based ORF coordinate of the neighbouring exon to measure from:
#          its Start for a + strand gene, its End for a - strand gene
#  $dir    "+" or "-", the transcription direction
#Choose a short LAST exon from exact-match candidates, anchored on the exon
#BEFORE it. pick_short_exon() handles the mirror image - a short FIRST exon
#anchored on the one after it - and its intron arithmetic assumes that geometry,
#so the two are kept separate rather than overloading one sub with a mode flag.
#
#$cands are 0-based plastome offsets; the exon occupies [c+1, c+$elen] 1-based.
sub pick_last_exon {
	my ($cands, $elen, $anchor, $dir) = @_;
	return(undef) unless $cands && @$cands;
	my @ok;
	for my $c (@$cands){
		my $ilen;
		if($dir eq '-'){
			#the preceding exon lies at HIGHER coordinates on this strand
			next unless $c + $elen < $anchor;
			$ilen = $anchor - ($c + $elen) - 1;
		}
		else{
			next unless $c + 1 > $anchor;
			$ilen = ($c + 1) - $anchor - 1;
		}
		next if $ilen < 200 || $ilen > 3000;      #a plausible plastid intron
		push @ok, [ $ilen, $c ];
	}
	return(undef) unless @ok;
	my ($best) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @ok;
	return($best->[1]);
}

sub pick_short_exon {
	my ($cands, $elen, $anchor, $dir) = @_;
	return(undef) unless $cands && @$cands;
	return($cands->[0]) if @$cands == 1;
	return($cands->[0]) unless defined $anchor && defined $elen && $elen > 0;

	my (@with_donor, @all);
	for my $c (@$cands){
		push @all, [ abs($c - $anchor), $c ];

		#Intron bounds in 1-based plastome coordinates. For a + strand gene the
		#intron follows the exon; for a - strand gene it precedes it genomically.
		my ($istart, $iend);
		if($dir eq '-'){ $istart = $anchor + 2;      $iend = $c; }
		else           { $istart = $c + $elen + 1;   $iend = $anchor; }
		next if $iend < $istart;
		my $ilen = $iend - $istart + 1;
		next if $ilen < 200 || $ilen > 3000;      #a plausible plastid intron
		push @with_donor, [ abs($c - $anchor), $c ]
			if &splice_donor_ok($istart, $iend, $dir);
	}

	my $pool = @with_donor ? \@with_donor : \@all;
	my ($best) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$pool;
	return($best->[1]);
}
###########################
###########################
#Locate the ORF-protein residue that aligns with the reference's C-terminal
#anchor. The mirror of locate_ref_start.
#
#The anchor is reference residue L-4, matching what the original code probed with
#(`substr(substr($gene_prot,-4),0,3)` is residues L-4..L-2) and what $nearend is
#used as downstream - a search start for match_end.
#
#Two problems with the original. Three residues are not specific, exactly as at
#the 5' end. And the miss path, `for (my $i=2; $nearend<0; $i++)`, has no bound:
#substr($gene_prot,-3-$i) does not run off the front and return undef, it returns
#the WHOLE string once $i passes the length, so the same 3-mer is retried forever.
#Once close guides put many references in play that hung the annotator outright.
#
#Here the probe is grown backwards from the anchor - longer, so more specific -
#and the index is shifted back to the anchor so the returned value keeps the
#original meaning. Returns -1 if nothing matches.
sub locate_ref_end {
	my ($orf_prot, $gene_prot, $from) = @_;
	return(-1) unless defined $orf_prot && defined $gene_prot;
	my $L = length($gene_prot);
	return(-1) if $L < 4 || !length($orf_prot);
	$from = 0 unless defined $from && $from >= 0;

	my $anchor = $L - 4;                      #the reference residue $nearend names
	my $fallback = -1;
	for my $back (9, 7, 5, 3, 1, 0){
		my $begin = $anchor - $back;
		next if $begin < 0;
		my $probe = substr($gene_prot, $begin, 3 + $back);
		next unless defined $probe && length($probe) == 3 + $back;
		next if $probe =~ /_/;                #never probe across a stop codon
		my $hit = index($orf_prot, $probe, $from);
		$hit = index($orf_prot, $probe) if $hit < 0;
		next if $hit < 0;
		my $idx = $hit + $back;               #shift forward to the anchor
		return($idx) if index($orf_prot, $probe, $hit + 1) < 0;
		$fallback = $idx if $fallback < 0;
	}
	return($fallback) if $fallback >= 0;

	#Nothing anchored at the tail: walk a 3-mer inwards, bounded by the reference.
	for (my $i = 1; $i <= 30 && $anchor - $i >= 0; $i++){
		my $probe = substr($gene_prot, $anchor - $i, 3);
		next unless defined $probe && length($probe) == 3;
		my $hit = index($orf_prot, $probe);
		return($hit + $i) if $hit >= 0;
	}
	return(-1);
}
###########################
#Locate where the reference's N-terminus sits in the ORF protein, as a residue
#index into $orf_prot.
#
#This used to probe with three residues of the reference starting at residue 1.
#Three amino acids are not specific: in a 150-residue ORF a given 3-mer turns up
#by chance, and a spurious early hit pushed $nearstart tens of codons upstream.
#Everything downstream is bounded by that value - match_start and alt_start only
#search at or before $nearstart*3 - so the real start became unreachable. That is
#how atpF exon 2 lost 57 nt in Marchantia and 22 nt in Pinus.
#
#Take the longest probe that matches, and prefer one that matches exactly once.
#Returns the residue index of the reference's residue 0, or -1.
sub locate_ref_start {
	my ($orf_prot, $gene_prot) = @_;
	return(-1) unless defined $orf_prot && defined $gene_prot;
	return(-1) if length($gene_prot) < 2 || length($orf_prot) < 1;

	my $fallback = -1;
	for my $len (12, 10, 8, 6, 5, 4, 3){
		last if length($gene_prot) < 1 + $len;
		my $probe = substr($gene_prot, 1, $len);
		next if $probe =~ /_/;                  #never probe across a stop codon
		my $first = index($orf_prot, $probe);
		next if $first < 0;
		#A probe that occurs once is trustworthy; take it immediately.
		return($first - 1) if index($orf_prot, $probe, $first + 1) < 0;
		$fallback = $first - 1 if $fallback < 0;
	}
	return($fallback) if $fallback >= 0;

	#Nothing anchored from residue 1: slide the probe further into the reference.
	my $limit = length($gene_prot) - 3;
	$limit = 30 if $limit > 30;
	for (my $i = 2; $i <= $limit; $i++){
		my $hit = index($orf_prot, substr($gene_prot, $i, 3));
		return($hit - $i) if $hit >= 0;
	}
	return(-1);
}
###########################
#Nudge cis-spliced exon junctions onto real splice sites.
#
#Exon boundaries are inherited from whichever guide genome won the score, so a
#junction can sit 1-2 nt off. The tell is that both sides move together: exon1's
#3' end and exon2's 5' start are wrong by the same amount, i.e. the intron is the
#right length but displaced. Shifting the whole junction therefore leaves the CDS
#length - and the reading frame - untouched.
#
#Plastid group II introns are GT...AY, NOT the spliceosomal GT...AG. Surveyed over
#the 7-genome truth set: 5' is GT in 90/135 introns, and the 3' dinucleotide is
#AC (51) or AT (46) against a single AG. The tRNA introns of trnI, trnA and trnL
#are a different class and are left alone.
#
#Conservative by construction: a junction that already reads GT...AY is never
#touched, and a shift is only taken when it produces an exact GT...AY.
sub refine_intron_boundaries {
	return unless defined $plastome && length $plastome;

	#Gather exons per (gene, strand) out of the coordinates already emitted.
	my %exons;
	for my $s (sort {$a <=> $b} keys %final_annotation){
		for my $e (sort {$a <=> $b} keys %{$final_annotation{$s}}){
			for my $d (sort keys %{$final_annotation{$s}{$e}}){
				my $name = $final_annotation{$s}{$e}{$d};
				next unless defined $name;
				next unless $name =~ /^(.+)_exon(\d+)(?:XXX\d+)?$/;
				my ($gene, $num) = ($1, $2);
				next if $gene =~ /^trn/;          #group I / tRNA introns, different signal
				push @{$exons{$gene}{$d}}, { s => $s, e => $e, num => $num };
			}
		}
	}

	my $moved = 0;
	for my $gene (sort keys %exons){
		for my $d (sort keys %{$exons{$gene}}){
			my @ex = sort { $a->{s} <=> $b->{s} } @{$exons{$gene}{$d}};
			next unless @ex >= 2;

			for my $i (0 .. $#ex - 1){
				my ($left, $right) = ($ex[$i], $ex[$i+1]);

				#Only consecutive exons of one copy bound a real intron. Without this
				#an IR-duplicated gene pairs copy 1's last exon with copy 2's first
				#across the whole genome, and trans-spliced rps12 pairs exons that are
				#not spliced to each other at all.
				my $DBG = $ENV{ANNOBTD_DEBUG_SPLICE};
				if($DBG && abs($left->{num} - $right->{num}) != 1){
					print STDERR "splice: $gene skip (exon nums $left->{num},$right->{num})\n";
				}
				next unless abs($left->{num} - $right->{num}) == 1;

				my $istart = $left->{e} + 1;         #1-based, inclusive
				my $iend   = $right->{s} - 1;
				my $ilen   = $iend - $istart + 1;
				#Real plastid introns in the reference set run 304-2559 nt.
				if($DBG && ($ilen < 200 || $ilen > 3000)){
					print STDERR "splice: $gene skip (intron $ilen nt)\n";
				}
				next if $ilen < 200 || $ilen > 3000;

				#Move only to a strictly better splice signal; a junction that is
				#already the best on offer stays where it is. Gating instead on
				#"donor must newly become GT" was tried and measured no better, so
				#the simpler rule stands.
				my $here = &splice_score($istart, $iend, $d);
				my ($best, $best_score);
				for my $delta (-1, 1, -2, 2, -3, 3){
					my ($ns, $ne) = ($istart + $delta, $iend + $delta);
					next if $ns < 1 || $ne > length($plastome);
					next if $left->{s} >= $ns || $ne >= $right->{e};   #don't invert an exon
					my $sc = &splice_score($ns, $ne, $d);
					if($sc > $here && (!defined $best_score || $sc > $best_score)){
						($best, $best_score) = ($delta, $sc);
					}
				}
				if($DBG){
					my @sc = map { my ($a,$b)=($istart+$_,$iend+$_);
						($a<1||$b>length($plastome)) ? "x" : &splice_score($a,$b,$d) } (-3..3);
					printf STDERR "splice: %-16s %s %d-%d (%d nt) scores[-3..3]=%s -> %s\n",
						$gene, $d, $istart, $iend, $ilen, join(",",@sc),
						(defined $best ? "shift $best" : "no move");
				}
				next unless defined $best;

				#Move both sides of the junction by the same amount.
				&move_annotation($left->{s},  $left->{e},  $d, $left->{s},           $left->{e} + $best);
				&move_annotation($right->{s}, $right->{e}, $d, $right->{s} + $best,  $right->{e});
				$left->{e}  += $best;
				$right->{s} += $best;
				$moved++;
			}
		}
	}
	print STDERR "refine_intron_boundaries: shifted $moved junction(s) onto GT...AY\n" if $moved;
}
###########################
#True when the span opens with the group II donor GT on the transcribed strand.
#This is the strong half of the signal and the only one allowed to move a boundary.
sub splice_donor_ok {
	my ($istart, $iend, $dir) = @_;
	return 0 if $iend - $istart + 1 < 10;
	my $seq = uc substr($plastome, $istart - 1, $iend - $istart + 1);
	if($dir eq '-'){
		$seq = reverse $seq;
		$seq =~ tr/ACGTacgt/TGCAtgca/;
	}
	return substr($seq, 0, 2) eq 'GT' ? 1 : 0;
}
###########################
#Score how well a genomic span reads as a group II intron on the transcribed
#strand. Higher is better; 0 means no splice signal at either end.
#
#A hard GT...AY test is too strict for real plastid introns - tobacco rpl2 is
#genuinely GT...AA, and rejecting it left the boundary uncorrected. Across the
#7-genome truth set the donor is GT in 90/135 introns, while the acceptor ends in
#A in 104/135 but in AY in only 97. So the donor is scored as the stronger signal
#and the acceptor is graded: AY best, any A still credited.
sub splice_score {
	my ($istart, $iend, $dir) = @_;
	return -1 if $iend - $istart + 1 < 10;
	my $seq = uc substr($plastome, $istart - 1, $iend - $istart + 1);
	if($dir eq '-'){
		$seq = reverse $seq;
		$seq =~ tr/ACGTacgt/TGCAtgca/;
	}
	my $score = 0;
	$score += 2 if substr($seq, 0, 2) eq 'GT';           #donor
	my $acceptor = substr($seq, -2);
	$score += 2 if $acceptor =~ /^A[CT]$/;               #canonical AY
	$score += 1 if $acceptor =~ /^A.$/ && $acceptor !~ /^A[CT]$/;
	return $score;
}
###########################
#Re-key one %final_annotation entry from (start,end) to (new_start,new_end).
sub move_annotation {
	my ($s, $e, $d, $ns, $ne) = @_;
	return if $s == $ns && $e == $ne;
	return unless exists $final_annotation{$s}{$e}{$d};
	my $name = $final_annotation{$s}{$e}{$d};
	delete $final_annotation{$s}{$e}{$d};
	delete $final_annotation{$s}{$e} unless %{$final_annotation{$s}{$e}};
	delete $final_annotation{$s}    unless %{$final_annotation{$s}};
	$final_annotation{$ns}{$ne}{$d} = $name;
}
###########################
#Walk a tRNA/rRNA HSP out to the full length of its reference gene.
#
#Hits are accepted at >=97% of the reference length, so up to 3% of the gene can
#lie outside the alignment. Storing the raw query coordinates truncates the
#annotation by exactly that unaligned remainder, which is where the persistent
#1-2 nt slips on tRNAs and rRNAs came from.
#
#  $qs,$qe  query (plastome) start/end - blastn always reports $qs < $qe
#  $ss,$se  subject (reference gene) start/end - $ss > $se means a minus-strand hit
#  $reflen  full length of the reference gene
#  $glen    plastome length, for clamping (optional)
#Returns the extended ($qs,$qe).
sub extend_rna_hit{
	my ($qs, $qe, $ss, $se, $reflen, $glen) = @_;
	return($qs, $qe) unless defined $reflen && $reflen > 0;

	my ($missing_low, $missing_high);
	if($se >= $ss){                          #plus/plus: qs<->ss, qe<->se
		$missing_low  = $ss - 1;             #reference 5' end sits at the qs side
		$missing_high = $reflen - $se;       #reference 3' end sits at the qe side
	}
	else{                                    #plus/minus: qs<->ss (high), qe<->se (low)
		$missing_low  = $reflen - $ss;       #reference 3' end sits at the qs side
		$missing_high = $se - 1;             #reference 5' end sits at the qe side
	}
	$missing_low  = 0 if $missing_low  < 0;
	$missing_high = 0 if $missing_high < 0;

	$qs -= $missing_low;
	$qe += $missing_high;
	$qs = 1 if $qs < 1;
	$qe = $glen if defined $glen && $glen && $qe > $glen;
	return($qs, $qe);
}
###########################
sub rna_blasts{

my ($counter,$file, $idorfs, $torf_pos, $bestmatchfile, $genome_len) = @_;
my %idorfs = %$idorfs;
my %torf_pos = %$torf_pos;
my %bestmatch;
open my $tfile, "<", $bestmatchfile;
while(<$tfile>){
		chomp;
		my @tarray = split/\s+/;
		$bestmatch{$tarray[1]}=1;
}

open $tfile, "<", $file; #blast output
while(<$tfile>){
		chomp;
		my $cur_orfname = "orfRNA" . $counter;
		my @tarray = split /\s+/;
		unless(exists $bestmatch{$tarray[1]}){
			next;
		}
		my $hit_truelen=length($blastrnas{$tarray[1]});
		if (($tarray[9]-$tarray[8])<0){
			$cur_orfname = $cur_orfname . "-r";
		}
		my $hit_curlen=abs(($tarray[9]-$tarray[8]))+1;
		if($hit_curlen/$hit_truelen >= 0.96 && $tarray[2] >= 85){   #v2.6: add an identity floor
				my ($qs, $qe) = &extend_rna_hit($tarray[6], $tarray[7],
					$tarray[8], $tarray[9], $hit_truelen, $genome_len);

				$idorfs{$tarray[1]}{$cur_orfname}{"Start"}=$qs;
				$idorfs{$tarray[1]}{$cur_orfname}{"End"}=$qe;
				$torf_pos{$cur_orfname}{"Start"}=$qs;
				$torf_pos{$cur_orfname}{"End"}=$qe;
				$counter++;


		}
}
return($counter,\%idorfs,\%torf_pos);
close $file;
}
###########################
sub translate {
	my @seqs;
	my $count=0;
	my $codon;
	for my $seq (@_){
        my $protein;

        for(my $i=0;$i<(length($seq)-2);$i+=3){
                $codon=substr($seq,$i,3);
                #$codon= uc $codon;
                if (exists $codons{$codon}){
                    $protein .= $codons{$codon};
                }
                
                unless(exists $codons{$codon}){

                  	$protein .= "X";
                }
                
               
        }

        $seqs[$count]=$protein;
        $count++;
	}
	return @seqs;
}
###########################
sub best_match {
	my $total=0;
	my $forward=0;
	my $reverse=0;
	for (my $i=0; $i<length($_[0])-2; $i+=3){
		my $test = substr($_[0],$i,3);
		$total++;
		if(index($_[1], $test) >= 0){
			$forward++
		}
		if(index($_[2], $test) >= 0){
			$reverse++;
		}
	}
	my $fortot = $forward/$total;
	my $revtot = $reverse/$total;
	return($fortot,$revtot);
}
###########################
sub best_match_frames {
	my $total=0;
	my $forward_0=0;
	my $reverse_0=0;
	my $forward_1=0;
	my $reverse_1=0;
	my $forward_2=0;
	my $reverse_2=0;
	for (my $i=0; $i<length($_[0])-2; $i+=3){
		my $test = substr($_[0],$i,3);
		$total++;
		if(index($_[1], $test) >= 0){
			$forward_0++
		}
		if(index($_[2], $test) >= 0){
			$reverse_0++;
		}
		if(index($_[3], $test) >= 0){
			$forward_1++
		}
		if(index($_[4], $test) >= 0){
			$reverse_1++;
		}
		if(index($_[5], $test) >= 0){
			$forward_2++
		}
		if(index($_[6], $test) >= 0){
			$reverse_2++;
		}
	}
	my $for_frame=0;
	my $rev_frame=0;
	my $fortot = $forward_0/$total;
	if(($forward_1/$total)>$fortot){
		$fortot = $forward_1/$total;
		$for_frame=1;
	}
	if(($forward_2/$total)>$fortot){
		$fortot = $forward_2/$total;
		$for_frame=2;
	}
	my $revtot = $reverse_0/$total;
	if(($reverse_1/$total)>$revtot){
		$revtot = $reverse_1/$total;
		$rev_frame=1;
	}
	if(($reverse_2/$total)>$revtot){
		$revtot = $reverse_2/$total;
		$rev_frame=2;
	}

	return($fortot,$for_frame,$revtot,$rev_frame);
}
###########################
#Locate the ATG start codon of a CDS inside an ORF.
# $_[0] ORF nucleotides, already shifted into the matching reading frame
# $_[1] $nearstart - offset in AMINO ACIDS where the reference N-terminus aligns,
#       so the reference-implied start sits at about $nearstart*3 nt
# $_[2] reference gene nucleotides (unused here; kept for the call signature)
#
#Only offsets that are a multiple of 3 are legal starts: the previous version
#scanned every offset and kept the LAST hit, so it happily returned an
#out-of-frame ATG, and biased starts downstream when several were in range.
sub match_start{
	my ($seq, $nearstart, $gene_seq) = @_;
	return(-1) unless defined $seq;
	my $len = length($seq);
	return(-1) if $len < 3;

	my $target = (defined $nearstart ? $nearstart : 0) * 3;
	$target = 0    if $target < 0;
	$target = $len if $target > $len;

	#Closest in-frame ATG at or before the reference-implied start.
	my $best = -1;
	for (my $i = 0; $i + 3 <= $len && $i <= $target; $i += 3){
		$best = $i if substr($seq, $i, 3) eq "ATG";
	}
	return($best) if $best >= 0;

	#Nothing upstream: take the first in-frame ATG just downstream instead.
	my $limit = $target + 30;
	$limit = $len - 3 if $limit > $len - 3;
	for (my $i = $target - ($target % 3); $i <= $limit; $i += 3){
		next if $i < 0;
		return($i) if substr($seq, $i, 3) eq "ATG";
	}
	return(-1);
}
###########################
#As match_start, but looks for the reference's own start codon rather than ATG,
#so non-ATG plastid starts (ACG, GTG) are recovered. Same in-frame constraint.
sub alt_start{
	my ($seq, $nearstart, $gene_seq) = @_;
	return(-1) unless defined $seq && defined $gene_seq;
	my $len = length($seq);
	my $codon = substr($gene_seq, 0, 3);
	return(-1) if $len < 3 || length($codon) < 3;

	my $target = (defined $nearstart ? $nearstart : 0) * 3;
	$target = 0    if $target < 0;
	$target = $len if $target > $len;

	my $best = -1;
	for (my $i = 0; $i + 3 <= $len && $i <= $target; $i += 3){
		$best = $i if substr($seq, $i, 3) eq $codon;
	}
	return($best);
}
###########################
#Locate a CDS 5' boundary when BOTH match_start and alt_start have failed.
#
#The old fallback was pure arithmetic on $nearstart, and with $nearstart 0 it
#simply kept the ORF's own end. That put Agave virginica ndhK's start on a TAA -
#a CDS beginning with a stop codon - 39 nt (13 codons) outside the annotated
#boundary, giving an internal stop and a HIGH flag.
#
#Scan in frame from the ORF's 5' end inward instead and take the first ATG.
#ATG is preferred over the other table-11 initiators rather than taking whichever
#comes first, because it is 96.9% of starts measured over 2,248 curated plastid
#CDS, and because the first initiator encountered here is the wrong one: ndhK's
#ORF meets ATT at 714 nt before the true ATG at 702. Other initiators are used
#only when no ATG is found. Returns undef when nothing qualifies, leaving the
#caller's original arithmetic in place.
#Should scan_start_5p be consulted at all? Only when the boundary in hand is
#already suspect - either it sits on a stop codon, or the length it implies
#disagrees with the reference by more than 5%.
#
#Without this gate the scan corrupts a start that is RIGHT but unrecognisable
#from sequence. Sorghum bicolor rpl23 begins TAC and the record declares
#/transl_except=(pos:59411..59413,aa:Met); the arithmetic already placed it
#exactly, and hunting for an ATG moved it 12 nt inward. Its implied length
#matches the reference, so the gate leaves it alone. Arabidopsis psbK (273 nt
#against a ~186 nt reference) and Agave ndhK (a TAA start) both still qualify.
#Third scoring signal: a per-gene EXPECTATION mined from thousands of GenBank
#records (build_gene_profiles.pl), consulted only where the first two signals -
#the guide reference and the splice site - have already failed or tied. A guide
#can be wrong in a way that is copied faithfully; against hundreds of records of
#the same gene that error is a minority. Stratified by family, then order, then
#everything, so lineage-specific lengths (infA 77 aa in eudicots, 107 in Poaceae)
#are simply the expectation for that lineage rather than an anomaly.
#
#It is an annotation consensus, not expression evidence, and a lineage whose
#records share a propagated error will carry it (Poaceae ycf3 sits at 172 aa on
#the strength of 27 related records; every other family reads 168-170). The
#cross-lineage fallback is what exposes that; the family cell alone cannot.
sub expected_nt {
	my ($gid, $cur) = @_;   # $cur: current CDS length, enables weak mode
	(my $name = $gid) =~ s/XXX.*//;
	print STDERR "TRACE expected_nt $gid lineage=$EXP_FAM/$EXP_ORD loaded=".scalar(keys %EXPECT)."\n" if $ENV{ANNOBTD_TRACE};
	for my $key ($name, ($name =~ /^(.+)_\d+$/ && $1 !~ /_exon$/ ? $1 : ())) {
		next unless $EXPECT{$key};
		# Family first; the order only when the family is thin (n < MIN_N). A populated
		# but unsettled family cell means the gene is variable or the annotation
		# convention is split inside the family, and a sibling family's convention
		# must not be imposed on it (Nicotiana ndhD: Solanaceae 1503/1383 split,
		# Convolvulaceae-only order cell said 1530). No ALL fallback.
		# The genus refines only where the family is AMBIGUOUS (thin or unsettled): a
		# settled family cell is kept as is. Measured genus-first at 5 species: +4/-2 on
		# LOO, the harms being Daucus ycf1 (family settled at 5475, genus 5454 over a gene
		# with codon-scale indels) and Zea ndhK (5 species, all the modern convention).
		my $fam = $EXPECT{$key}{"F:$EXP_FAM"}; my $fam_settled = $fam && $fam->[1] >= $EXP_MIN_N && $fam->[2] >= $EXP_MIN_CONS;
		if (!$fam_settled && $EXP_GEN ne '' && $EXP_GEN ne 'NA') { my $g = $EXPECT{$key}{"G:$EXP_GEN"};
			if ($g && $g->[1] >= $EXP_MIN_N_GENUS) { return $g->[0] if $g->[2] >= $EXP_MIN_CONS; return undef } }   # populated but unsettled genus: no expectation
		for my $lin ("F:$EXP_FAM", "O:$EXP_ORD") {
			next if $lin =~ /:NA$/;
			my $e = $EXPECT{$key}{$lin}; next unless $e && $e->[1] >= $EXP_MIN_N;
			return $e->[0] if $e->[2] >= $EXP_MIN_CONS;
			if ($cur && $e->[2] >= $EXP_WEAK_CONS) {
				my $near = sub { my $m = shift; $m > 0 && abs($cur - $m) / $m <= 0.02 };
				unless ($near->($e->[0]) || $near->($e->[3])) {
					print STDERR "PROFILE weak $gid cur $cur modes $e->[0]/$e->[3] cons $e->[2]\n" if $ENV{ANNOBTD_TRACE};
					return $e->[0];
				}
			}
			last;   # populated but unsettled, and the call sits on one of its modes
		}
	}
	return undef;
}

#Profile as a third signal at the sites where a locator SUCCEEDED. Measured on 17
#genomes, 40 emitted boundaries were wrong AND >5% off a settled expectation, while
#29 were exact-but-unusual; so this only ever moves a boundary when every check
#agrees: the cell is settled (expected_nt's gate), the current length is >5% off,
#and an in-frame candidate exists within 2% of the expectation that begins on an
#initiator and reaches the 3' end without an internal stop. ATG is preferred and a
#current ATG is never traded for a non-ATG. Returns the new 5' coordinate or undef.
sub profile_adjust_5p {
	my ($gid, $lo, $hi, $dir) = @_;
	my $cur = $hi - $lo + 1; return undef if $cur <= 0;
	my $exp = &expected_nt($gid, $cur); return undef unless $exp && $exp > 0;
	return undef if abs($cur - $exp) / $exp <= 0.05;
	my $five = $dir eq '+' ? $lo : $hi; my $three = $dir eq '+' ? $hi : $lo;
	my $codon_at = sub { my $p = shift; return undef if $p < 3 || $p + 2 > length $plastome;
		$dir eq '+' ? uc substr($plastome, $p - 1, 3) : do { my $z = substr($plastome, $p - 3, 3); $z = reverse $z; $z =~ tr/ACGTacgt/TGCAtgca/; uc $z } };
	my $clean = sub { my $p = shift; my ($a, $b) = $dir eq '+' ? ($p, $three) : ($three, $p); return 0 if $b < $a;
		my $seq = substr($plastome, $a - 1, $b - $a + 1); if ($dir eq '-') { $seq = reverse $seq; $seq =~ tr/ACGTacgt/TGCAtgca/ }
		$seq = uc $seq; for (my $i = 0; $i + 5 < length $seq; $i += 3) { return 0 if substr($seq, $i, 3) =~ /^(TAA|TAG|TGA)$/ } 1 };
	my $curc = $codon_at->($five) // '';
	my $reach = abs($exp - $cur) + 30; my @cand;
	for (my $step = -$reach; $step <= $reach; $step += 3) { next unless $step;
		my $p = $dir eq '+' ? $five + $step : $five - $step; my $len = abs($three - $p) + 1; next if $len < 30 || $len % 3;
		next if abs($len - $exp) / $exp > 0.02;
		my $c = $codon_at->($p) // next; next unless $c =~ /^(ATG|GTG|TTG|ATT|ATC|ATA|CTG|ACG)$/;
		next if $curc eq 'ATG' && $c ne 'ATG';
		next unless $clean->($p);
		push @cand, [ abs($len - $exp), ($c eq 'ATG' ? 0 : 1), abs($step), $p, $c, $len ] }
	return undef unless @cand;
	my ($b) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @cand;
	print STDERR "PROFILE adjust $gid $dir 5' $five->$b->[3] len $cur->$b->[5] exp $exp start $curc->$b->[4]\n" if $ENV{ANNOBTD_TRACE};
	return $b->[3];
}

sub start_is_suspect {
	my ($codon, $curlen, $reflen) = @_;
	return 1 if defined $codon && $codon =~ /^(TAA|TAG|TGA)$/;
	return 0 unless $reflen && $curlen;
	return abs($curlen - $reflen) / $reflen > 0.05 ? 1 : 0;
}

sub scan_start_5p {
	my ($plast, $three_prime, $dir, $orf_edge, $reflen, $expect) = @_;
	return undef unless defined $orf_edge && defined $three_prime;
	#With a profile expectation, the candidate whose length best matches it wins,
	#ATG breaking ties. Without one, the first ATG inward (the earlier rule). This
	#is what separates Schoenolirion ndhK (882 nt, an ATG at 702 also available)
	#from Agave ndhK (702 nt): the family expectation decides, not scan order.
	my @cand;
	#Do not shrink the CDS below 60% of the reference; past that the scan is
	#guessing rather than correcting.
	my $floor = defined $reflen && $reflen > 0 ? int(0.6 * $reflen) : 60;
	my @init;
	for (my $step = 0; $step <= 400; $step += 3){
		my ($pos, $len);
		if($dir eq '+'){ $pos = $orf_edge + $step; $len = $three_prime - $pos + 1 }
		else           { $pos = $orf_edge - $step; $len = $pos - $three_prime + 1 }
		last if $len < $floor;
		next if $len % 3;
		my $codon = ($dir eq '+')
			? substr($plast, $pos - 1, 3)
			: do { my $z = substr($plast, $pos - 3, 3); $z = reverse $z; $z =~ tr/ACGTacgt/TGCAtgca/; $z };
		$codon = uc $codon;
		next unless length($codon) == 3;
		next if $codon =~ /^(TAA|TAG|TGA)$/;      #never begin a CDS on a stop
		my $is_init = ($codon eq 'ATG' || $codon =~ /^(GTG|TTG|ATT|ATC|ATA|CTG|ACG)$/);
		next unless $is_init;
		if (!defined $expect) { return $pos if $codon eq 'ATG'; push @init, $pos; next }
		push @cand, [ abs($len - $expect), ($codon eq 'ATG' ? 0 : 1), $step, $pos ];
	}
	if (@cand) { my ($b) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @cand; return $b->[3] }
	return $init[0] if @init;
	return undef;
}

sub match_end{
	#my $nearend = length($orf_prot)-index(reverse($orf_prot), substr($gene_prot,-1))-1;;
	#my $endpos=length($orf_prot)-1;
	my $tempend;
	my $seq = $_[0];
	if(index($seq, substr($_[1],-4)) > 0){
		$tempend= index($seq, substr($_[1],-4),$_[2])+3;
	}
	else{
	$tempend = index($seq,substr($_[1],-1),$_[2]);
	my $boundary;
	if(index($seq, "_", $_[2]) > 0){
		if(index($seq, "_", $_[2]) > $_[2]+5){
			$boundary = $_[2]+5;
		}
		else{
			$boundary=index($seq, "_", $_[2]);
		}
	}
	else{
		$boundary=length($seq);
	}

	for (my $i=$tempend; $i<=$boundary;$i++){
		if(index($seq,substr($_[1],-1),$i)>=0){#second parameter is gene_seq
			if(index($seq,substr($_[1],-1),$i) < $boundary){
				$tempend =index($seq,substr($_[1],-1),$i);

				
			}
		}
	}
}
	#$seq is the PROTEIN, so this offset is measured on 3*residues - it does not
	#include any trailing partial codon. The stop-codon path in the callers
	#measures on length($orf_seq) instead, so the callers normalise the two onto
	#one scale rather than this sub guessing.
	$tempend = length($seq)*3-($tempend*3 + 2)-1;
	return($tempend);

}
###########################
sub exon_mods_start {
	my $frame_restore = $_[0];
	my $gid = $_[1]; #pass $gid to function
	my $eid = $_[2]; #pass $eid to function
	my $gene_seq=substr($blastgenes{$gid},0, length($blastgenes{$gid})-$frame_restore);
	my $orf_seq = substr($plastome,$idorfs_f{$gid}{$eid}{"Start"},($idorfs_f{$gid}{$eid}{"End"}-$idorfs_f{$gid}{$eid}{"Start"}+1));
	my $rev_orf_seq = reverse($orf_seq);
	$rev_orf_seq =~ tr/ATCGatcg/TAGCtagc/;
	my $orf_seq_1=substr($orf_seq,1);
	my $rev_orf_seq_1 = substr(reverse($orf_seq),1);
	$rev_orf_seq_1 =~ tr/ATCGatcg/TAGCtagc/;
	my $orf_seq_2=substr($orf_seq,2);
	my $rev_orf_seq_2 = substr(reverse($orf_seq),2);
	$rev_orf_seq_2 =~ tr/ATCGatcg/TAGCtagc/;
	my ($orf_prot, $rev_orf_prot, $orf_prot_1, $rev_orf_prot_1, $orf_prot_2, $rev_orf_prot_2, $gene_prot)=&translate($orf_seq,$rev_orf_seq,$orf_seq_1,$rev_orf_seq_1,$orf_seq_2,$rev_orf_seq_2,$gene_seq);
	my ($formatch, $for_frame,$revmatch,$rev_frame)=&best_match_frames($gene_prot, $orf_prot, $rev_orf_prot,$orf_prot_1,$rev_orf_prot_1,$orf_prot_2,$rev_orf_prot_2);
	if($formatch < 0.5 && $revmatch < 0.5){
		return;
	}
	if($formatch>$revmatch){
		if($for_frame == 1){
			$orf_seq = $orf_seq_1;
			$orf_prot = $orf_prot_1;
		}
		if($for_frame == 2){
			$orf_seq = $orf_seq_2;
			$orf_prot = $orf_prot_2;
		}
		my $mod_shift = length($orf_seq)%3;
		my $nearstart = &locate_ref_start($orf_prot, $gene_prot);
		$nearstart = 0 if $nearstart < 0;   #a residue index, never negative
		my $nearend = &locate_ref_end($orf_prot, $gene_prot, $nearstart);
		$nearend = 0 if $nearend < 0;   #a residue index, never negative
=item		for (my $i=0; $i<(length($orf_prot)-1);$i++){
		if(index($orf_prot,substr(substr($gene_prot,-4),0,3),$i)>=0){#third parameter is gene_seq
			if(index($orf_prot,substr(substr($gene_prot,-4),0,3),$i) < (length($orf_prot))){
				$nearend =index($orf_prot,substr(substr($gene_prot,-4),0,3),$i);
				$i=$nearend;
			}
			
			}
		}
=cut	
		my $end_match;
		my $stop_at = (substr($gene_prot, -1) eq "_") ? index($orf_prot,"_", $nearstart) : -1;
		if($stop_at >= 0){
			$end_match = length($orf_seq)-($stop_at*3 + 2)-1;
		}
		else{
			#The reference carries a stop but this ORF has none at or after
			#$nearstart - the ORF finder stopped just short of it. Guarding this
			#matters: index() returns -1 and the old arithmetic turned that into
			#length($orf_seq)-(-3+2)-1, i.e. the whole ORF length, which put the 3'
			#end back on the ORF's own start. Marchantia petD_exon2 came out as
			#"73216 73204 +" that way - a 13 nt reversed span for a 475 nt exon.
			#Falling through to match_end aligns the reference's own tail instead.
			$end_match = &match_end($orf_prot,$gene_prot,$nearend);
			#match_end measures on 3*residues; the stop path above measures on
			#length($orf_seq). Add the trailing partial codon so both are on the
			#same scale and one coordinate formula serves both.
			$end_match += length($orf_seq) % 3;
		}
		my $start_match = &match_start($orf_seq,$nearstart,$gene_seq);
		my $start_alt = &alt_start($orf_seq,$nearstart,$gene_seq);
		if(substr($gene_seq,0,3) eq substr($orf_seq,$start_match,3)){
				$start_alt = -1;
		}


		my $ref_target = $nearstart * 3;
		my $use_alt = 0;
		if($start_alt >= 0){
			$use_alt = ($start_match < 0)
				|| (abs($start_alt - $ref_target) < abs($start_match - $ref_target));
		}
		if($start_match == -1 && $start_alt == -1){
			&dbg_coord("exon_mods_start","$gid","$eid",1,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","+");
{
				my $lo = $idorfs_f{$gid}{$eid}{"Start"}+$for_frame+($nearstart*3)-3+1;
				my $hi = $idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1;
				my $cc = uc substr($plastome, $lo - 1, 3);
				my $exp_nt = &expected_nt($gid);
				if(&start_is_suspect($cc, $hi - $lo + 1, $exp_nt // length($blastgenes{$gid}))){
					my $scan = &scan_start_5p($plastome, $hi, "+", $lo, length($blastgenes{$gid}), $exp_nt);
					$lo = $scan if defined $scan;
				}
				$final_annotation{$lo}{$hi}{"+"}=$gid;
			}
		}
		elsif(!$use_alt){
			&dbg_coord("exon_mods_start","$gid","$eid",2,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","+");
			{
			my $lo = $idorfs_f{$gid}{$eid}{"Start"}+$for_frame+$start_match+1;
			my $hi = $idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1;
			my $adj = &profile_adjust_5p($gid, $lo, $hi, "+");
			$lo = $adj if defined $adj;
			$final_annotation{$lo}{$hi}{"+"}=$gid;
		}
		}
		else{
			&dbg_coord("exon_mods_start","$gid","$eid",3,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","+");
			{
			my $lo = $idorfs_f{$gid}{$eid}{"Start"}+$for_frame+$start_alt+1;
			my $hi = $idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1;
			my $adj = &profile_adjust_5p($gid, $lo, $hi, "+");
			$lo = $adj if defined $adj;
			$final_annotation{$lo}{$hi}{"+"}=$gid;
		}
		}
	}
							
	else{
		if($rev_frame == 1){
			$rev_orf_seq = $rev_orf_seq_1;
			$rev_orf_prot = $rev_orf_prot_1;
		}
		if($rev_frame == 2){
			$rev_orf_seq = $rev_orf_seq_2;
			$rev_orf_prot = $rev_orf_prot_2;

		}
		my $mod_shift = length($rev_orf_seq)%3;
		my $nearstart = &locate_ref_start($rev_orf_prot, $gene_prot);
		$nearstart = 0 if $nearstart < 0;   #a residue index, never negative
		my $nearend = &locate_ref_end($rev_orf_prot, $gene_prot, $nearstart);
		$nearend = 0 if $nearend < 0;   #a residue index, never negative
		
		my $end_match;
		my $stop_at = (substr($gene_prot, -1) eq "_") ? index($rev_orf_prot,"_", $nearstart) : -1;
		if($stop_at >= 0){
			$end_match = length($rev_orf_seq)-($stop_at*3 + 2)-1;
			$mod_shift=0;
		}
		else{
			#See the forward branch: a missing stop must not reach the arithmetic
			#as -1. $mod_shift is left alone here, as it is on the no-stop path.
			$end_match = &match_end($rev_orf_prot,$gene_prot, $nearend);
			#match_end measures on 3*residues; the stop path above measures on
			#length($rev_orf_seq). Add the trailing partial codon so both are on the
			#same scale and one coordinate formula serves both.
			$end_match += length($rev_orf_seq) % 3;
		}
		my $start_match = &match_start($rev_orf_seq,$nearstart,$gene_seq);
		my $start_alt = &alt_start($rev_orf_seq,$nearstart,$gene_seq);
		if(substr($gene_seq,0,3) eq substr($rev_orf_seq,$start_match,3)){
				$start_alt = -1;
		}
		my $ref_target = $nearstart * 3;
		my $use_alt = 0;
		if($start_alt >= 0){
			$use_alt = ($start_match < 0)
				|| (abs($start_alt - $ref_target) < abs($start_match - $ref_target));
		}
		if($start_match == -1 && $start_alt == -1){
			&dbg_coord("exon_mods_start","$gid","$eid",4,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","-");
			{
				my $lo = $idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1;
				my $hi = $idorfs_f{$gid}{$eid}{"End"}-($nearstart*3)+3+1;
				my $cc = do { my $z = substr($plastome, $hi - 3, 3); $z = reverse $z; $z =~ tr/ACGTacgt/TGCAtgca/; uc $z };
				my $exp_nt = &expected_nt($gid);
				if(&start_is_suspect($cc, $hi - $lo + 1, $exp_nt // length($blastgenes{$gid}))){
					my $scan = &scan_start_5p($plastome, $lo, "-", $hi, length($blastgenes{$gid}), $exp_nt);
					$hi = $scan if defined $scan;
				}
				$final_annotation{$lo}{$hi}{"-"}=$gid;
			}
		}
		elsif(!$use_alt){
			&dbg_coord("exon_mods_start","$gid","$eid",5,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","-");
			{
			my $lo = $idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1;
			my $hi = $idorfs_f{$gid}{$eid}{"End"}-$start_match-$rev_frame+1;
			my $adj = &profile_adjust_5p($gid, $lo, $hi, "-");
			$hi = $adj if defined $adj;
			$final_annotation{$lo}{$hi}{"-"}=$gid;
		}
		}
		else{
			&dbg_coord("exon_mods_start","$gid","$eid",6,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","-");
			{
			my $lo = $idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1;
			my $hi = $idorfs_f{$gid}{$eid}{"End"}-$start_alt-$rev_frame+1;
			my $adj = &profile_adjust_5p($gid, $lo, $hi, "-");
			$hi = $adj if defined $adj;
			$final_annotation{$lo}{$hi}{"-"}=$gid;
		}
		}
	}

}
###########################
sub exon_mods_end {
	my $frame_restore = $_[0];
	my $gid = $_[1]; #pass $gid to function
	my $eid = $_[2]; #pass $eid to function
	my $gene_seq=substr($blastgenes{$gid},$frame_restore, length($blastgenes{$gid}));
	my $orf_seq = substr($plastome,$idorfs_f{$gid}{$eid}{"Start"},($idorfs_f{$gid}{$eid}{"End"}-$idorfs_f{$gid}{$eid}{"Start"}+1));
	my $rev_orf_seq = reverse($orf_seq);
	$rev_orf_seq =~ tr/ATCGatcg/TAGCtagc/;
	my $orf_seq_1=substr($orf_seq,1);
	my $rev_orf_seq_1 = substr(reverse($orf_seq),1);
	$rev_orf_seq_1 =~ tr/ATCGatcg/TAGCtagc/;
	my $orf_seq_2=substr($orf_seq,2);
	my $rev_orf_seq_2 = substr(reverse($orf_seq),2);
	$rev_orf_seq_2 =~ tr/ATCGatcg/TAGCtagc/;
	my ($orf_prot, $rev_orf_prot, $orf_prot_1, $rev_orf_prot_1, $orf_prot_2, $rev_orf_prot_2, $gene_prot)=&translate($orf_seq,$rev_orf_seq,$orf_seq_1,$rev_orf_seq_1,$orf_seq_2,$rev_orf_seq_2,$gene_seq);
	my ($formatch, $for_frame,$revmatch,$rev_frame)=&best_match_frames($gene_prot, $orf_prot, $rev_orf_prot,$orf_prot_1,$rev_orf_prot_1,$orf_prot_2,$rev_orf_prot_2);

	if($formatch>$revmatch){
		if($for_frame == 1){
			$orf_seq = $orf_seq_1;
			$orf_prot = $orf_prot_1;
		}
		if($for_frame == 2){
			$orf_seq = $orf_seq_2;
			$orf_prot = $orf_prot_2;
		}
		my $nearend = &locate_ref_end($orf_prot, $gene_prot, 0);
		$nearend = 0 if $nearend < 0;   #a residue index, never negative
		my $end_match;
		#index() returns -1 when this ORF carries no stop, and feeding that to the
		#arithmetic below yields length($orf_seq)-(-3+2)-1 = the whole ORF length,
		#which puts the 3' end back on the ORF's own start. See exon_mods_start.
		my $stop_at = (substr($gene_prot, -1) eq "_") ? index($orf_prot,"_", $nearend) : -1;
		if($stop_at >= 0){
			$end_match = length($orf_seq)-($stop_at*3 + 2)-1;
		}
		else{
			$end_match = &match_end($orf_prot,$gene_prot,$nearend);
			#match_end measures on 3*residues; the stop path above measures on
			#length($orf_seq). Add the trailing partial codon so both are on the
			#same scale and one coordinate formula serves both.
			$end_match += length($orf_seq) % 3;
		}
		my $nearstart = &locate_ref_start($orf_prot, $gene_prot);
		$nearstart = 0 if $nearstart < 0;   #a residue index, never negative
		#my $start_match = &match_start($orf_seq);
		my $start_alt = &alt_start($orf_seq,$nearstart,$gene_seq);
		#my $end_match = &match_end($orf_seq);

		
		&dbg_coord("exon_mods_end","$gid","$eid",7,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,"",$nearstart,$nearend,"",$start_alt,$end_match,"","","+");
		#alt_start returns -1 when it cannot place the 5' boundary. Letting that reach
		#the formula shifts the boundary by one instead of trimming it, so the exon
		#keeps whatever length the ORF happened to have. The reference's own length is
		#a far better estimate: the 3' end is already known, so measure back from it.
		#See the reverse branch for the case this was found on.
		{
			my $hi = $idorfs_f{$gid}{$eid}{"End"} - $end_match + 1;
			my $lo = ($start_alt >= 0)
			       ? $idorfs_f{$gid}{$eid}{"Start"} + $start_alt - $frame_restore + 1 + $for_frame
			       : $hi - length($gene_seq) + 1;
			#Deliberately NOT clamped to the ORF span. An exon boundary can legitimately
			#fall outside it - the ORF finder's ends are stop-to-stop, not exon ends -
			#and clamping cost rps16_exon2, whose 5' boundary sits 3 nt past its ORF.
			$final_annotation{$lo}{$hi}{"+"}=$gid;
		}
		
	}
							
	else{
		if($rev_frame == 1){
			$rev_orf_seq = $rev_orf_seq_1;
			$rev_orf_prot = $rev_orf_prot_1;
		}
		if($rev_frame == 2){
			$rev_orf_seq = $rev_orf_seq_2;
			$rev_orf_prot = $rev_orf_prot_2;

		}
		#my $start_match = &match_start($rev_orf_seq);
		my $nearstart = &locate_ref_start($rev_orf_prot, $gene_prot);
		$nearstart = 0 if $nearstart < 0;   #a residue index, never negative
		my $nearend = &locate_ref_end($rev_orf_prot, $gene_prot, $nearstart);
		$nearend = 0 if $nearend < 0;   #a residue index, never negative
		my $end_match;
		#See the forward branch: a missing stop must not reach the arithmetic as -1.
		my $stop_at = (substr($gene_prot, -1) eq "_") ? index($rev_orf_prot,"_",$nearstart) : -1;
		if($stop_at >= 0){
			$end_match = length($rev_orf_seq)-($stop_at*3 + 2)-1;
		}
		else{
			$end_match = &match_end($rev_orf_prot,$gene_prot, $nearstart);
			#match_end measures on 3*residues; the stop path above measures on
			#length($rev_orf_seq). Add the trailing partial codon so both are on the
			#same scale and one coordinate formula serves both.
			$end_match += length($rev_orf_seq) % 3;
		}
		my $start_alt = &alt_start($rev_orf_seq,$nearstart,$gene_seq);
		#my $end_match = &match_end($rev_orf_seq);

		&dbg_coord("exon_mods_end","$gid","$eid",8,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,"",$nearstart,$nearend,"",$start_alt,$end_match,"","","-");
		#As the forward branch. This is the case it was found on: Agave virginica
		#clpP_exon3 sits on orf4710-r, 71917-72333, against a 252 nt reference.
		#alt_start returned -1, so the 5' boundary came out as End+1 = 72335 instead
		#of being trimmed, giving a 418 nt exon whose protein carried five internal
		#stops. Measuring 252 nt back from the known 3' end gives 72169 - the
		#annotated boundary exactly. Schoenolirion croceum fails and recovers the same
		#way.
		{
			my $lo = $idorfs_f{$gid}{$eid}{"Start"} + 1 + $end_match;
			my $hi = ($start_alt >= 0)
			       ? $idorfs_f{$gid}{$eid}{"End"} - $start_alt + $frame_restore + 1 - $rev_frame
			       : $lo + length($gene_seq) - 1;
			#Deliberately NOT clamped to the ORF span. An exon boundary can legitimately
			#fall outside it - the ORF finder's ends are stop-to-stop, not exon ends -
			#and clamping cost rps16_exon2, whose 5' boundary sits 3 nt past its ORF.
			$final_annotation{$lo}{$hi}{"-"}=$gid;
		}
		
	}

}
###########################
sub exon_mods_mid{
	my $frame_restore = $_[0];
	my $gid = $_[1]; #pass $gid to function
	my $eid = $_[2]; #pass $eid to function
	my $gene_seq=$blastgenes{$gid};
	my $frame=0;
	#v2.6: search all three frames and take the one that reads longest before a
	#stop. The old loop only searched 0..$frame_restore and mis-scored a frame
	#with no stop at all as position -1.
	my $max_length_frame=0;
	for (my $j=0; $j <= 2; $j++){
		my @cur_frame = &translate(substr($gene_seq, $j));
		if($cur_frame[0] =~ /_/){
			if(index($cur_frame[0], "_") > $max_length_frame){
				$max_length_frame = index($cur_frame[0], "_");
				$frame = $j;
			}
		}
		else{
			if(length($cur_frame[0]) > $max_length_frame){
				$max_length_frame = length($cur_frame[0]);
				$frame = $j;
			}
		}
	}
	my $end_remove = length(substr($gene_seq,$frame))%3;
	$gene_seq = substr($blastgenes{$gid}, $frame,length(substr($gene_seq,$frame))-$end_remove);
	my $orf_seq = substr($plastome,$idorfs_f{$gid}{$eid}{"Start"},($idorfs_f{$gid}{$eid}{"End"}-$idorfs_f{$gid}{$eid}{"Start"}+1));
	my $rev_orf_seq = reverse($orf_seq);
	$rev_orf_seq =~ tr/ATCGatcg/TAGCtagc/;
	my $orf_seq_1=substr($orf_seq,1);
	my $rev_orf_seq_1 = substr(reverse($orf_seq),1);
	$rev_orf_seq_1 =~ tr/ATCGatcg/TAGCtagc/;
	my $orf_seq_2=substr($orf_seq,2);
	my $rev_orf_seq_2 = substr(reverse($orf_seq),2);
	$rev_orf_seq_2 =~ tr/ATCGatcg/TAGCtagc/;
	my ($orf_prot, $rev_orf_prot, $orf_prot_1, $rev_orf_prot_1, $orf_prot_2, $rev_orf_prot_2, $gene_prot)=&translate($orf_seq,$rev_orf_seq,$orf_seq_1,$rev_orf_seq_1,$orf_seq_2,$rev_orf_seq_2,$gene_seq);
	my ($formatch, $for_frame,$revmatch,$rev_frame)=&best_match_frames($gene_prot, $orf_prot, $rev_orf_prot,$orf_prot_1,$rev_orf_prot_1,$orf_prot_2,$rev_orf_prot_2);

	
	if($formatch>$revmatch){
		if($for_frame == 1){
			$orf_seq = $orf_seq_1;
			$orf_prot = $orf_prot_1;
		}
		if($for_frame == 2){
			$orf_seq = $orf_seq_2;
			$orf_prot = $orf_prot_2;
		}
		my $nearstart = &locate_ref_start($orf_prot, $gene_prot);
		$nearstart = 0 if $nearstart < 0;   #a residue index, never negative
		my $nearend = &locate_ref_end($orf_prot, $gene_prot, $nearstart);
		$nearend = 0 if $nearend < 0;   #a residue index, never negative
		my $start_alt = &alt_start($orf_seq,$nearstart,$gene_seq);
		my $end_match = &match_end($orf_prot,$gene_prot,$nearend);
		&dbg_coord("exon_mods_mid","$gid","$eid",9,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame,$end_remove,$nearstart,$nearend,"",$start_alt,$end_match,"","","+");
		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_alt-$frame+1}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$end_remove+1}{"+"}=$gid;
	}
	else{
		if($rev_frame == 1){
			$rev_orf_seq = $rev_orf_seq_1;
			$rev_orf_prot = $rev_orf_prot_1;
		}
		if($rev_frame == 2){
			$rev_orf_seq = $rev_orf_seq_2;
			$rev_orf_prot = $rev_orf_prot_2;

		}
		my $nearstart = &locate_ref_start($rev_orf_prot, $gene_prot);
		$nearstart = 0 if $nearstart < 0;   #a residue index, never negative

		my $nearend = &locate_ref_end($rev_orf_prot, $gene_prot, $nearstart);
		$nearend = 0 if $nearend < 0;   #a residue index, never negative
		my $start_alt = &alt_start($rev_orf_seq,$nearstart,$gene_seq);
		my $end_match = &match_end($rev_orf_prot,$gene_prot,$nearend);
		&dbg_coord("exon_mods_mid","$gid","$eid",10,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame,$end_remove,$nearstart,$nearend,"",$start_alt,$end_match,"","","-");
		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}-$end_remove+1+$end_match}{$idorfs_f{$gid}{$eid}{"End"}-$start_alt+$frame+1}{"-"}=$gid;
	}
}
###########################
# Species-prefixed: unprefixed names collide when two annotations run at once.
open my $out, ">", $ARGV[9] . "_blast_orfid.txt";
for my $bid (sort keys %idorfs_f){
	for my $oid (sort keys %{$idorfs_f{$bid}}){
		print $out "$bid\t$oid\n";
	}
}

open my $outfile2, ">", $ARGV[9] . "_VERDANT_cleaned_annotation.txt";

for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (sort {$a <=> $b} keys %{$final_annotation{$start}}){
		for my $dir (sort keys %{$final_annotation{$start}{$end}}){
			my $out_gene;
			if($final_annotation{$start}{$end}{$dir} !~ /XXX/){
				$out_gene = $final_annotation{$start}{$end}{$dir};
			}
			else{
				$out_gene = $final_annotation{$start}{$end}{$dir};
				$out_gene =~ /(.*?)XXX.+/;
				$out_gene = $1;
			}
			#Drop the _2 / _3 disambiguator that get_annotated_regions_fromverdant.pl
			#adds to tell same-named loci apart. It is an internal reference id, not
			#a gene name: all three trnS loci are reported as trnS.
			#Guarded so _exon2 and similar are untouched.
			$out_gene =~ s/_\d+$// if defined $out_gene;
			next unless defined $out_gene && length $out_gene;

			# A record whose end precedes its start is malformed and must not reach
			# the output: coordinates are always ascending here, on both strands,
			# with $dir carrying the orientation. Downstream anything computing
			# end-start gets a negative length from these.
			#
			# Four occurred across the nine-genome regression, and none was a
			# formatting slip - each was a call that had already failed:
			#   Spinacia  "Start~trnH 1 0"      trnH begins at position 1, so the
			#                                   leading spacer is empty
			#   Pinus     "psbA 1028 976"       spurious; the real psbA spans the
			#                                   origin at 1-976 / 119622-119707
			#   Pinus     "rpl2_exon2 94944 94711"  spurious duplicate; the real one
			#                                   is called correctly at 63799-64245
			#   Marchantia "petD_exon2 73216 73204"  a 13 nt reversed span where the
			#                                   truth is 73216-73690
			# So dropping them removes two spurious features and one 474 nt-wrong
			# boundary, and costs nothing that was right. They are reported on
			# stderr rather than silently discarded, because a malformed record
			# means an upstream boundary computation went wrong and that is worth
			# seeing even though the output is now clean.
			if($end < $start){
				my $kind = ($out_gene =~ /~/) ? "empty spacer" : "reversed span";
				warn sprintf("%s: dropping malformed record (%s): %s %d-%d %s\n",
				             $ARGV[9], $kind, $out_gene, $start, $end, $dir);
				next;
			}
			print $outfile2 "$out_gene\t$start\t$end\t$dir\n";
		}
	}
}
$lsc_range =~ /(\d+)-(\d+)/;
print $outfile2 "LSC\t$1\t$2\t+\n";
$irb_range =~ /(\d+)-(\d+)/;
print $outfile2 "IRB\t$1\t$2\t+\n";
$ssc_range =~ /(\d+)-(\d+)/;
print $outfile2 "SSC\t$1\t$2\t+\n";
$ira_range =~ /(\d+)-(\d+)/;
print $outfile2 "IRA\t$1\t$2\t+\n";
$full_range =~ /(\d+)-(\d+)/;
print $outfile2 "FULL\t$1\t$2\t+\n";
