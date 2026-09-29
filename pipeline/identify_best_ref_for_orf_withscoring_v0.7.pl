#!/usr/bin/perl -w
use strict;

open(STDOUT, '>', $ARGV[4] . '_results.log') or die "Can't open log";
open(STDERR, '>', $ARGV[4] . '_results_error.log') or die "Can't open log";

my %blastgenes;
my %orfseqs;
my %refs2orfs;
my %orfs2refs;
#Clean up blast files for annoBTD
# 1: blastfile 2: blastgeneseqs 3:query orfs/seqs 4:Feature type
#Read in ORFs sequences and reference tRNA, rRNA, and gene sequences
%blastgenes = &read_fasta($ARGV[1],%blastgenes); #These are the REFERENCE genes for annotation
%orfseqs   = &read_fasta($ARGV[2],%orfseqs);  #query seqs (ORFs, or the plastome for RNA runs)
my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');
my $feature=$ARGV[3];

# Scoring knobs, read once. See annotation_score() for why k defaults to 2.
my $SCORE_K   = defined $ENV{ANNOBTD_SCORE_K}   ? $ENV{ANNOBTD_SCORE_K}   : 2;
my $SCORE_MIN = defined $ENV{ANNOBTD_SCORE_MIN} ? $ENV{ANNOBTD_SCORE_MIN} : 0.70;
my $RANK_K    = defined $ENV{ANNOBTD_RANK_K}    ? $ENV{ANNOBTD_RANK_K}    : 3;


my %blast_gene_hsps;
if($feature eq "protein"){
	my ($current_hsp_query,$current_hsp_subject,%best_orf_ref,%orf_ref_scores,$dir);

	open my $file, "<", $ARGV[0] or die "Cannot open BLAST output $ARGV[0]: $!"; #blast output
	while(<$file>){
		chomp;
		if(/intron/){
			next;
		}
		my @tarray = split /\s+/;
		if($tarray[2]< 0.5){
				next;
		}
		if($tarray[1] =~ /rpl16/ && $tarray[0] !~ /\-r/){
			next;
		}
		if(exists $blast_gene_hsps{$tarray[0]}){
			if(exists $blast_gene_hsps{$tarray[0]}{$tarray[1]}){
				if($tarray[9] > $tarray[8]){
					for(my $i=$tarray[8]-1; $i<=$tarray[9]-1; $i++){
						$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=1;
					}
				}
				else{
					for(my $i=$tarray[9]-1; $i<=$tarray[8]-1; $i++){
						$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=1;
					}
				}

			}
			else{
				for (my $i=0; $i<=length($blastgenes{$tarray[1]})-1; $i++){
					$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=0;
				}
				if($tarray[9] > $tarray[8]){
					for(my $i=$tarray[8]-1; $i<=$tarray[9]-1; $i++){
						$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=1;
					}
				}
				else{
					for(my $i=$tarray[9]-1; $i<=$tarray[8]-1; $i++){
						$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=1;
					}
				}

			}
		}
		else{
			if($current_hsp_query){
				
				for my $orfid (sort keys %blast_gene_hsps){
					my $best_match;
					my $best_score=0;
					for my $refgene (sort keys %{$blast_gene_hsps{$orfid}}){
						if($refgene =~ /rpl16/){
								if($orfid !~ /-r/){
										next;
								}
						}
					####Get sequence for current reference and orf######
						my $orf_seq = $orfseqs{$orfid};
						my $ref_seq = $blastgenes{$refgene};

					####Identify percent match to reference for current orf-ref pair####
					
					####Translate orf in all direction to identify best match####
						my $orf_seq2 = substr($orf_seq,1);
						my $orf_seq3 = substr($orf_seq,2);
						my $rev_orf_seq = reverse($orf_seq);
						$rev_orf_seq =~ tr/ATCGatcg/TAGCtagc/;
						my $rev_orf_seq2 = substr($rev_orf_seq,1);
						my $rev_orf_seq3 = substr($rev_orf_seq,2);
						
						my $trans_orf = translate($orf_seq);
						my $trans_orf2 = translate($orf_seq2);
						my $trans_orf3 = translate($orf_seq3);
						my $trans_rev_orf = translate($rev_orf_seq);
						my $trans_rev_orf2 = translate($rev_orf_seq2);
						my $trans_rev_orf3 = translate($rev_orf_seq3);
						
						my $trans_ref = translate($ref_seq);
						#Pick the reference reading frame that runs furthest before a stop.
						#Exon references do not necessarily begin in frame 0 - an exon's frame
						#depends on the exons before it - so without this, petD_exon2,
						#rps12_exon3 and rps16_exon2 score too low to pass the 0.70 cutoff and
						#drop out of the annotation entirely.
						#The original wrote $max_transreflen from $trans_ref in every branch
						#rather than from the frame under test; this compares them directly.
						{
							my $best_run = -1;
							for my $cand ($trans_ref, translate(substr($ref_seq,1)), translate(substr($ref_seq,2))){
								next unless defined $cand && length $cand;
								my $stop = index($cand, "_");
								my $run  = ($stop == -1) ? length($cand) : $stop;
								if($run > $best_run){ $best_run = $run; $trans_ref = $cand }
							}
						}

					####Calculate annotation score for all directions####
						my @scores;
						push(@scores, annotation_score($trans_orf, $trans_ref));
						push(@scores, annotation_score($trans_orf2, $trans_ref));
						push(@scores, annotation_score($trans_orf3, $trans_ref));
						push(@scores, annotation_score($trans_rev_orf, $trans_ref));
						push(@scores, annotation_score($trans_rev_orf2, $trans_ref));
						push(@scores, annotation_score($trans_rev_orf3, $trans_ref));

						my $max_score = -1;
						for my $tscore (@scores){
								if($tscore > $max_score){
										$max_score = $tscore;
								}
						}
						#Admission uses $max_score (k=2). The choice between competing
						#references uses the more specific k=$RANK_K over the same frames.
						my $rank_score = -1;
						for my $fr ($trans_orf,$trans_orf2,$trans_orf3,$trans_rev_orf,$trans_rev_orf2,$trans_rev_orf3){
								my $r = annotation_score($fr, $trans_ref, $RANK_K);
								$rank_score = $r if $r > $rank_score;
						}
						print "$refgene\t$max_score\n";
						if($max_score < $SCORE_MIN){
							#next;
							delete $blast_gene_hsps{$orfid}{$refgene};
							delete $refs2orfs{$refgene}{$orfid};
							delete $orfs2refs{$orfid}{$refgene};
							next;
				#delete $orfseqs{$orfid};
						}
						$orfid =~ /(.*?)XXX.+/;
						my $true_gene = $1;
		####NEED TO KEEP THE RIGHT REFERENCE ID BUT LOOK AT THEM ALL SEPARATELY. 
						if(exists $orf_ref_scores{$orfid}){

							if ($rank_score > $orf_ref_scores{$orfid}{Rank}
							    || ($rank_score == $orf_ref_scores{$orfid}{Rank}
							        && $max_score > $orf_ref_scores{$orfid}{Score})){
								#Drop the previous claim. refs2orfs has to stay a true
								#inverse of the current best assignment: overlap groups are
								#built per reference and the loser of a group is deleted, so
								#a stale entry lets an ORF be judged - and destroyed - under
								#a gene it no longer claims.
								my $prev_ref = $orf_ref_scores{$orfid}{Real_Ref};
								if(defined $prev_ref && $prev_ref ne $refgene){
									delete $refs2orfs{$prev_ref}{$orfid};
									delete $orfs2refs{$orfid}{$prev_ref};
								}
								$orf_ref_scores{$orfid}{ORF}=$orfid;
								$orf_ref_scores{$orfid}{Score}=$max_score;
								$orf_ref_scores{$orfid}{Rank}=$rank_score;
								$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
								$refs2orfs{$refgene}{$orfid}=1;
								$orfs2refs{$orfid}{$refgene}=1;
							}
						}
						else{
							$orf_ref_scores{$orfid}{ORF}=$orfid;
							$orf_ref_scores{$orfid}{Score}=$max_score;
							$orf_ref_scores{$orfid}{Rank}=$rank_score;
							$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
							$refs2orfs{$refgene}{$orfid}=1;
							$orfs2refs{$orfid}{$refgene}=1;
						}

					}
				  }	
				}
				%blast_gene_hsps = ();
				for (my $i=0; $i<=length($blastgenes{$tarray[1]})-1; $i++){
					$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=0;
				}
				if($tarray[9] > $tarray[8]){
					for(my $i=$tarray[8]-1; $i<=$tarray[9]-1; $i++){
						$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=1;
					}
				}
				else{
					for(my $i=$tarray[9]-1; $i<=$tarray[8]-1; $i++){
						$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$i}=1;
					}
				}
				$current_hsp_query = $tarray[0];
				$current_hsp_subject = $tarray[1];
			}
		}
#close $file;

	for my $orfid (sort keys %blast_gene_hsps){
					my $best_match;
					my $best_score=0;
					for my $refgene (sort keys %{$blast_gene_hsps{$orfid}}){
					####Get sequence for current reference and orf######
						my $orf_seq = $orfseqs{$orfid};
						my $ref_seq = $blastgenes{$refgene};
					####Translate orf in all direction to identify best match####
						my $orf_seq2 = substr($orf_seq,1);
						my $orf_seq3 = substr($orf_seq,2);
						my $rev_orf_seq = reverse($orf_seq);
						$rev_orf_seq =~ tr/ATCGatcg/TAGCtagc/;
						my $rev_orf_seq2 = substr($rev_orf_seq,1);
						my $rev_orf_seq3 = substr($rev_orf_seq,2);
						
						my $trans_orf = translate($orf_seq);
						my $trans_orf2 = translate($orf_seq2);
						my $trans_orf3 = translate($orf_seq3);
						my $trans_rev_orf = translate($rev_orf_seq);
						my $trans_rev_orf2 = translate($rev_orf_seq2);
						my $trans_rev_orf3 = translate($rev_orf_seq3);
						
						my $trans_ref = translate($ref_seq);
						#Pick the reference reading frame that runs furthest before a stop.
						#Exon references do not necessarily begin in frame 0 - an exon's frame
						#depends on the exons before it - so without this, petD_exon2,
						#rps12_exon3 and rps16_exon2 score too low to pass the 0.70 cutoff and
						#drop out of the annotation entirely.
						#The original wrote $max_transreflen from $trans_ref in every branch
						#rather than from the frame under test; this compares them directly.
						{
							my $best_run = -1;
							for my $cand ($trans_ref, translate(substr($ref_seq,1)), translate(substr($ref_seq,2))){
								next unless defined $cand && length $cand;
								my $stop = index($cand, "_");
								my $run  = ($stop == -1) ? length($cand) : $stop;
								if($run > $best_run){ $best_run = $run; $trans_ref = $cand }
							}
						}

					####Calculate annotation score for all directions####
						my @scores;
						push(@scores, annotation_score($trans_orf, $trans_ref));
						push(@scores, annotation_score($trans_orf2, $trans_ref));
						push(@scores, annotation_score($trans_orf3, $trans_ref));
						push(@scores, annotation_score($trans_rev_orf, $trans_ref));
						push(@scores, annotation_score($trans_rev_orf2, $trans_ref));
						push(@scores, annotation_score($trans_rev_orf3, $trans_ref));

						my $max_score = -1;
						for my $tscore (@scores){
								if($tscore > $max_score){
										$max_score = $tscore;
								}
						}
						#Admission uses $max_score (k=2). The choice between competing
						#references uses the more specific k=$RANK_K over the same frames.
						my $rank_score = -1;
						for my $fr ($trans_orf,$trans_orf2,$trans_orf3,$trans_rev_orf,$trans_rev_orf2,$trans_rev_orf3){
								my $r = annotation_score($fr, $trans_ref, $RANK_K);
								$rank_score = $r if $r > $rank_score;
						}
						#print "$max_score\n";
						if($max_score < $SCORE_MIN){
							#next;
							delete $blast_gene_hsps{$orfid}{$refgene};
							delete $orfs2refs{$orfid}{$refgene};
							#was: delete $orfs2refs{$refgene}{$orfid} - wrong hash, so a
							#reference the ORF had just FAILED kept the ORF in refs2orfs
							#and could still drag it into an overlap group and delete it.
							delete $refs2orfs{$refgene}{$orfid};
							next;
				#delete $orfseqs{$orfid};
						}
						$orfid =~ /(.*?)XXX.+/;
						my $true_gene = $1;
		####NEED TO KEEP THE RIGHT REFERENCE ID BUT LOOK AT THEM ALL SEPARATELY. 
						if(exists $orf_ref_scores{$orfid}){

							if ($rank_score > $orf_ref_scores{$orfid}{Rank}
							    || ($rank_score == $orf_ref_scores{$orfid}{Rank}
							        && $max_score > $orf_ref_scores{$orfid}{Score})){
								#Drop the previous claim. refs2orfs has to stay a true
								#inverse of the current best assignment: overlap groups are
								#built per reference and the loser of a group is deleted, so
								#a stale entry lets an ORF be judged - and destroyed - under
								#a gene it no longer claims.
								my $prev_ref = $orf_ref_scores{$orfid}{Real_Ref};
								if(defined $prev_ref && $prev_ref ne $refgene){
									delete $refs2orfs{$prev_ref}{$orfid};
									delete $orfs2refs{$orfid}{$prev_ref};
								}
								$orf_ref_scores{$orfid}{ORF}=$orfid;
								$orf_ref_scores{$orfid}{Score}=$max_score;
								$orf_ref_scores{$orfid}{Rank}=$rank_score;
								$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
								$refs2orfs{$refgene}{$orfid}=1;
								$orfs2refs{$orfid}{$refgene}=1;
							}
						}
						else{
							$orf_ref_scores{$orfid}{ORF}=$orfid;
							$orf_ref_scores{$orfid}{Score}=$max_score;
							$orf_ref_scores{$orfid}{Rank}=$rank_score;
							$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
							#was: keyed on $true_gene, from "$orfid =~ /(.*?)XXX.+/".
							#$orfid is an ORF id and never contains XXX, so that match
							#always fails and $1 silently holds a leftover from an earlier
							#regex. The other three call sites all key on $refgene.
							$refs2orfs{$refgene}{$orfid}=1;
							$orfs2refs{$orfid}{$refgene}=1;
						}

					}
				  }	

	my %coordinates;
	open my $coords, "<", $ARGV[5] or die "Cannot open coords $ARGV[5]: $!"; #Read in coordinates for ORFs
	while(<$coords>){
			chomp;
			my @tarray = split/\s+/;
			$coordinates{$tarray[0]}=$tarray[1];
	}
	
	# Redundant-ORF removal, among ORFs claiming the SAME reference gene.
	#
	# Only a COMPLETELY NESTED ORF is redundant. Two ORFs that merely share some
	# bases are not competing for the same thing: an ORF runs to the next stop
	# codon, so it routinely overruns its gene and laps into a neighbour that the
	# annotation does not overlap at all - rpoB and rpoC1 sit 37 nt apart in Zea
	# mays, psaA and psaB 25 nt. Real annotated genes do overlap, but by very
	# little: measured over 44 plastomes, 95 of 174 overlapping pairs share no more
	# than 5% of the shorter gene (atpE/atpB 4 nt, rps3/rpl22 16 nt, ndhK/ndhC a
	# mean of 20 nt). The exception is matK, which sits wholly inside the trnK
	# intron - but those are different references and never compete here anyway.
	#
	# The previous code grouped on ANY shared base, transitively, so a single
	# nucleotide of contact merged ORFs into one winner-take-all group and every
	# loser was deleted outright. Set ANNOBTD_OVERLAP_MODE=any to restore it.
	# Redundancy is judged by RECIPROCAL overlap, not containment alone. Two ORFs
	# can each stick out past the other and still be the same locus read in
	# different frames or on opposite strands: Arabidopsis psaI drew orf3217-r
	# (59116-59331) against orf4826 (59234-59359), 78% of the shorter, and
	# Marchantia psaJ drew a pair overlapping 97%. Neither is nested, both tie the
	# correct ORF exactly on score, and downstream picked the wrong one - psaI came
	# out as "59247 59117 +", start past end.
	#
	# This is safe for adjacent genes because the comparison is per reference:
	# atpE and atpB claim different references and never meet here. The only
	# legitimate distinct loci for ONE gene are the IR copies, tens of kb apart.
	my $OVMIN = defined $ENV{ANNOBTD_OVERLAP_MIN} ? $ENV{ANNOBTD_OVERLAP_MIN} : 0.5;
	if(($ENV{ANNOBTD_OVERLAP_MODE} || 'nested') eq 'nested'){
	  for my $tempref (sort keys %refs2orfs){
		my @orfs = grep { exists $orf_ref_scores{$_} && exists $coordinates{$_} }
		           sort keys %{$refs2orfs{$tempref}};
		next unless @orfs > 1;
		my %span;
		for my $o (@orfs){ my ($s,$e) = split /-/, $coordinates{$o}; $span{$o} = [$s,$e] }
		for my $a (@orfs){
			for my $b (@orfs){
				next if $a eq $b;
				next unless exists $orf_ref_scores{$a} && exists $orf_ref_scores{$b};
				# nested, or overlapping by at least $OVMIN of the shorter ORF
				my $lo = $span{$a}[0] > $span{$b}[0] ? $span{$a}[0] : $span{$b}[0];
				my $hi = $span{$a}[1] < $span{$b}[1] ? $span{$a}[1] : $span{$b}[1];
				my $ov = $hi - $lo + 1;
				next if $ov <= 0;
				my $la = $span{$a}[1] - $span{$a}[0] + 1;
				my $lb = $span{$b}[1] - $span{$b}[0] + 1;
				my $shorter = $la < $lb ? $la : $lb;
				next unless $ov >= $OVMIN * $shorter;
				# identical spans: settle it one way round only, so the pair does
				# not delete each other depending on iteration order
				next if $span{$a}[0] == $span{$b}[0] && $span{$a}[1] == $span{$b}[1] && $a ge $b;
				my $ra = defined $orf_ref_scores{$a}{Rank} ? $orf_ref_scores{$a}{Rank} : $orf_ref_scores{$a}{Score};
				my $rb = defined $orf_ref_scores{$b}{Rank} ? $orf_ref_scores{$b}{Rank} : $orf_ref_scores{$b}{Score};
				# Frames of one locus tie exactly on k-mer containment, so the tie
				# needs a real signal: prefer the ORF whose length sits closest to
				# the reference it claims.
				my $fita = abs($la - length($blastgenes{ $orf_ref_scores{$a}{Real_Ref} } // ''));
				my $fitb = abs($lb - length($blastgenes{ $orf_ref_scores{$b}{Real_Ref} } // ''));
				if($ra > $rb
				   || ($ra == $rb && $orf_ref_scores{$a}{Score} > $orf_ref_scores{$b}{Score})
				   || ($ra == $rb && $orf_ref_scores{$a}{Score} == $orf_ref_scores{$b}{Score} && $fita < $fitb)){
					delete $orf_ref_scores{$b};
				}
				else{
					delete $orf_ref_scores{$a};
					last;                 # $a is gone; stop comparing it
				}
			}
		}
	  }
	}
	else{
		for my $tempref (sort keys %refs2orfs){
			my %tempposition;
			my %orf_groups;
			###NEED TO add a way to group orfs into regions of the genome
			if(scalar (keys %{$refs2orfs{$tempref}}) > 1){
				for my $temporf (sort keys %{$refs2orfs{$tempref}}){
					my @tempcoord = split(/-/, $coordinates{$temporf});
					for (my $i=$tempcoord[0]; $i<=$tempcoord[1]; $i++){
							$tempposition{$i}++;
							push(@{$orf_groups{$i}}, $temporf);
					}
				}
			}
					
		my %temp_overlap;
		my $group_count=0;
	
		for my $ctemp (sort {$a <=> $b} keys %orf_groups){
			if(scalar (@{$orf_groups{$ctemp}}) >= 1){
				my $fill;
				for my $idt (@{$orf_groups{$ctemp}}){
					if(exists $temp_overlap{$group_count}){
					if(exists $temp_overlap{$group_count}{$idt}){
						$fill=1;
					}
				}
			}
	
				if($fill){
						for my $idt (@{$orf_groups{$ctemp}}){
							$temp_overlap{$group_count}{$idt}=1;
						}
				}
				else{
						$group_count++;
						for my $idt (@{$orf_groups{$ctemp}}){
							$temp_overlap{$group_count}{$idt}=1;
						}
				}
			
			}

		}
		for my $group_temp (sort keys %temp_overlap){
			for my $group_mem (sort keys %{$temp_overlap{$group_temp}}){
				for my $group_temp2 (sort keys %temp_overlap){
					if($group_temp2 eq $group_temp){
							next;
					}
					for my $group_mem2 (sort keys %{$temp_overlap{$group_temp}}){
						if ($group_mem2 eq $group_mem){
							$temp_overlap{$group_temp}=$temp_overlap{$group_temp2};
							delete $temp_overlap{$group_temp2};
						}
					}
				}
			}
		}
		#Overlapping ORFs compete and only one survives the group. This is a
		#DISCRIMINATION step, so it compares the specific k=$RANK_K score, not the
		#permissive admission score. Ranking here on the k=2 score threw away ndhK in
		#Zea mays and Sorghum bicolor: k=2 compresses everything towards 1.0, so a
		#neighbouring ORF edged out orf4460-r even though that ORF matched an ndhK
		#reference at 1.000. The admission score breaks ties.
		for my $group_temp (sort keys %temp_overlap){
			my $temp_max=0;
			my $temp_max_admit=0;
			my $prev_temp;
			for my $group_mem (sort keys %{$temp_overlap{$group_temp}}){
		   
					if(!exists $orf_ref_scores{$group_mem}){
						next;
					}
			    
					my $r = defined $orf_ref_scores{$group_mem}{Rank}
					      ? $orf_ref_scores{$group_mem}{Rank} : $orf_ref_scores{$group_mem}{Score};
					if($r > $temp_max
					   || ($r == $temp_max && $orf_ref_scores{$group_mem}{Score} > $temp_max_admit)){
						$temp_max = $r;
						$temp_max_admit = $orf_ref_scores{$group_mem}{Score};
						if($prev_temp){
							delete $orf_ref_scores{$prev_temp};
						}
						$prev_temp=$group_mem;
					}
					else{
						delete $orf_ref_scores{$group_mem};	

					}
				}
			}
		}
	}


	open my $out, ">", $ARGV[4] . "_best_orfs_for_refs_SCORE.txt";
	for my $ref_gene (sort keys %orf_ref_scores){
		#if($orf_ref_scores{$ref_gene}{Real_Ref} =~ /petDXXX/){
		#	$orf_ref_scores{$ref_gene}{Real_Ref} =~ /petD(XXX.*?)/;	
		#	$orf_ref_scores{$ref_gene}{Real_Ref} = "petD_exon2$1";
		#}
		print $out "$orf_ref_scores{$ref_gene}{ORF}\t$orf_ref_scores{$ref_gene}{Real_Ref}\t$orf_ref_scores{$ref_gene}{Score}\n";
	}
 }

else{
	my %best_orf_ref;
	my %orf_ref_scores;
	open my $file, "<", $ARGV[0] or die "Cannot open BLAST output $ARGV[0]: $!"; #blast output
	while(<$file>){
		chomp;
		if(/intron/){
			next;
		}
		my @tarray = split /\s+/;
		if($tarray[2]< 0.5){
				next;
		}
		my $dir;
		if($tarray[7]-$tarray[6] > 0 && $tarray[9]-$tarray[8] > 0){
			$dir="+";
		}
		elsif($tarray[7]-$tarray[6] < 0 && $tarray[9]-$tarray[8] < 0){
			$dir="+";
		}
		else{
			$dir = "-";
		}
		$blast_gene_hsps{$tarray[0]}{$tarray[1]}{$tarray[6]}{$tarray[7]}=$dir;
		my $ref_seq = $blastgenes{$tarray[1]};
		my $unknown = substr($orfseqs{$tarray[0]},$tarray[6]-1, abs($tarray[7]-$tarray[6]));
		if($dir eq "-"){
			$unknown = reverse($unknown);
			$unknown =~ tr/ATCGatcg/TAGCtagc/;
		}
		my $score = annotation_score_rna($unknown, $ref_seq);
		#Isoacceptors sharing a name used to collapse here, but the cause was
		#upstream: get_annotated_regions_fromverdant.pl overwrote same-named loci,
		#so only one trnS reference ever existed. With that fixed the references
		#arrive as trnS, trnS_2, trnS_3 and this non-greedy match already keeps
		#them apart. Keying by the full reference id instead was tried and was
		#worse: it keeps one entry per guide as well, so the same gene from every
		#guide is accepted and rna_blasts over-generates placements.
		$tarray[1] =~ /(.*?)XXX\d/;
		my $short_name=$1;
	    if(exists $orf_ref_scores{$short_name}){

			if ($score > $orf_ref_scores{$short_name}{Score}){
				$orf_ref_scores{$short_name}{ORF}=$tarray[0];
				$orf_ref_scores{$short_name}{Score}=$score;
				$orf_ref_scores{$short_name}{Real_Ref}=$tarray[1];
			}
		}
		else{
			$orf_ref_scores{$short_name}{ORF}=$tarray[0];
			$orf_ref_scores{$short_name}{Score}=$score;
			$orf_ref_scores{$short_name}{Real_Ref}=$tarray[1];
		}
	}


close $file;
if($feature eq "tRNA"){
	open my $out, ">", $ARGV[4] . "_best_tRNA_ref_SCORE.txt";
	for my $ref_gene (sort keys %orf_ref_scores){
		print $out "$orf_ref_scores{$ref_gene}{ORF}\t$orf_ref_scores{$ref_gene}{Real_Ref}\n";
	}
}
else{
	open my $out, ">", $ARGV[4] . "_best_rRNA_ref_SCORE.txt";
	for my $ref_gene (sort keys %orf_ref_scores){
		print $out "$orf_ref_scores{$ref_gene}{ORF}\t$orf_ref_scores{$ref_gene}{Real_Ref}\n";
	}
}

}





###########################
sub get_seq{
my ($file, $seqid) = @_;
open my $tempfasta, "<", $file;
my $tsid;
my $tempseq;
while(<$tempfasta>){
	chomp;
	if(/>/){
		my $tempid = substr($_,1);
		if($tempid eq $seqid){
			$tsid = $seqid;
		}
		else{
			$tsid = ();
		}
	}
	else{
		if($tsid){
		$tempseq.=$_;
	}}
}
#close $tempfasta;
return($tempseq);
}
###########################
sub read_fasta{
my ($file, %temphash) = @_;
my $sid;
open my $tfile, "<", $file or die "Cannot open FASTA $file: $!";
while(<$tfile>){
		chomp;
		if(/^>/){
				$sid = substr($_, 1);
				$sid =~ s/\s.*$//;   #BLAST truncates subject ids at the first space
		}
		else{
				$temphash{$sid}.=$_;

		}
}
close $tfile;
return(%temphash);
}
###########################
sub translate {
	my $seqs;
	my $count=0;
	my $codon;
	my $seq = $_[0];
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

        $seqs=$protein;
        $count++;
	
	return $seqs;
}
###########################
sub annotation_score{
	#USAGE: 1: translated ORF from the plastome being annotated  2: translated reference
	#
	#Fraction of the REFERENCE's k-peptides that occur somewhere in the ORF. This
	#asks the question the caller actually needs answered - does this ORF contain
	#this gene - and it is normalised by what was tested, so it cannot exceed 1.
	#
	#Two things changed here, both measured on the test set:
	#
	#  denominator. The original counted the ORF's tripeptides that appear in the
	#  reference but divided by the REFERENCE's length in residues. That mixes two
	#  different scales, so the score drifts with the ORF/reference length ratio
	#  and can exceed 1 (rps18 in Zea mays scored 1.012). Small in effect, but it
	#  made the number hard to reason about and to threshold.
	#
	#  k=2, not 3. An exact tripeptide match needs three consecutive identical
	#  residues, so the score falls off as roughly p^3 in the identity p. Against a
	#  close reference that is fine; against a distant one it collapses. Measured
	#  per gene, the fraction clearing the 0.70 cutoff was 100% for Zea mays, which
	#  has a congener in the set, but only 30% for Marchantia and 32% for Pinus,
	#  which have no close relative - and those two are the worst genomes in the
	#  regression at 42% and 46% recall. At k=2 they clear 86% and 82% while Zea
	#  stays at 100%. k=2 admits more weak pairs, so the cutoff still has to reject
	#  them; that is what ANNOBTD_SCORE_MIN is for.
	#
	#The two jobs this score does are separated at the call sites. ADMITTING a pair
	#("is there a gene here at all") uses k=2, which survives divergence. CHOOSING
	#between competing references ("which gene is it") uses ANNOBTD_RANK_K (3),
	#because k=2 cannot tell paralogues apart: ranking on k=2 alone let ndhJ and
	#ndhC outscore ndhK on ndhK's own ORF in both Zea mays and Sorghum bicolor, and
	#ndhK vanished from the annotation - its span came back as the ndhJ~ndhC spacer.
	#
	#Set ANNOBTD_SCORE_K=3 to restore the old k. The old denominator is not
	#restorable by a switch: it was wrong rather than a tunable choice.
	my ($orf, $ref, $k) = @_;
	return(0) unless defined $ref && defined $orf;
	$k = $SCORE_K unless defined $k;
	return(0) if length($ref) < $k;
	my ($total, $seq_score) = (0, 0);
	for (my $i=0; $i+$k <= length($ref); $i++){
		$total++;
		$seq_score++ if index($orf, substr($ref,$i,$k)) >= 0;
	}
	return($total ? $seq_score/$total : 0);
}
###########################
sub annotation_score_rna{
	#USAGE: 1: Sequence from plastome being annotated 2:Reference for sequence 3: Direction of sequence
	my $total=0;
	my $seq_score=0;
	for (my $i=0; $i<length($_[1])-5; $i+=6){
		my $test = substr($_[1],$i,6);
		$total++;
		
		if(index($_[0], $test) >= 0){
			$seq_score++;
		}
	}
	
	return(0) unless $total;
	my $final_score = $seq_score/$total;

	return($final_score);
}
