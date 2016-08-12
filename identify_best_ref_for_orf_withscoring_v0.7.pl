#!/usr/bin/perl -w
use strict;
use diagnostics;

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
my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');
my $feature=$ARGV[3];


my %blast_gene_hsps;
if($feature eq "protein"){
	my ($current_hsp_query,$current_hsp_subject,%best_orf_ref,%orf_ref_scores,$dir);

	open my $file, "<", $ARGV[0]; #blast output
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
				
				for my $orfid (keys %blast_gene_hsps){
					my $best_match;
					my $best_score=0;
					for my $refgene (keys %{$blast_gene_hsps{$orfid}}){
						if($refgene =~ /rpl16/){
								if($orfid !~ /-r/){
										next;
								}
						}
					####Get sequence for current reference and orf######
						my $orf_seq = &get_seq($ARGV[2], $orfid);
						my $ref_seq = &get_seq($ARGV[1], $refgene);

					####Identify percent match to reference for current orf-ref pair####
					
						my $overlap = 0;
						for my $hit (keys %{$blast_gene_hsps{$orfid}{$refgene}}){

						}	
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
						my $trans_ref2 = translate(substr($ref_seq,1));
						my $trans_ref3 = translate(substr($ref_seq,2));	

						my $max_transreflen = 0;
						my $temp_index = index($trans_ref,"_");
						if($temp_index == "-1"){
								$max_transreflen=100000000;
						}
						
						if(index($trans_ref,"_") > $max_transreflen){
								$max_transreflen = index($trans_ref,"_");
								$trans_ref = $trans_ref;
						}
						if(index($trans_ref2,"_") > $max_transreflen){
								$max_transreflen = index($trans_ref,"_");
								$trans_ref = $trans_ref2;
						}
						if(index($trans_ref3,"_") > $max_transreflen){
								$max_transreflen = index($trans_ref,"_");
								$trans_ref = $trans_ref3;
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
						print "$refgene\t$max_score\n";
						if($max_score < 0.70){
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

							if ($max_score > $orf_ref_scores{$orfid}{Score}){
								$orf_ref_scores{$orfid}{ORF}=$orfid;
								$orf_ref_scores{$orfid}{Score}=$max_score;
								$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
								$refs2orfs{$refgene}{$orfid}=1;
								$orfs2refs{$orfid}{$refgene}=1;
							}
						}
						else{
							$orf_ref_scores{$orfid}{ORF}=$orfid;
							$orf_ref_scores{$orfid}{Score}=$max_score;
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

	for my $orfid (keys %blast_gene_hsps){
					my $best_match;
					my $best_score=0;
					for my $refgene (keys %{$blast_gene_hsps{$orfid}}){
					####Get sequence for current reference and orf######
						my $orf_seq = &get_seq($ARGV[2], $orfid);
						my $ref_seq = &get_seq($ARGV[1], $refgene);
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
						my $trans_ref2 = translate(substr($ref_seq,1));
						my $trans_ref3 = translate(substr($ref_seq,2));	

						my $max_transreflen = 0;
						my $temp_index = index($trans_ref,"_");
						if($temp_index == "-1"){
								$max_transreflen=100000000;
						}
						if(index($trans_ref,"_") > $max_transreflen){
								$max_transreflen = index($trans_ref,"_");
								$trans_ref = $trans_ref;
						}
						if(index($trans_ref2,"_") > $max_transreflen){
								$max_transreflen = index($trans_ref,"_");
								$trans_ref = $trans_ref2;
						}
						if(index($trans_ref3,"_") > $max_transreflen){
								$max_transreflen = index($trans_ref,"_");
								$trans_ref = $trans_ref3;
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
						#print "$max_score\n";
						if($max_score < 0.70){
							#next;
							delete $blast_gene_hsps{$orfid}{$refgene};
							delete $orfs2refs{$orfid}{$refgene};
							delete $orfs2refs{$refgene}{$orfid};
							next;
				#delete $orfseqs{$orfid};
						}
						$orfid =~ /(.*?)XXX.+/;
						my $true_gene = $1;
		####NEED TO KEEP THE RIGHT REFERENCE ID BUT LOOK AT THEM ALL SEPARATELY. 
						if(exists $orf_ref_scores{$orfid}){

							if ($max_score > $orf_ref_scores{$orfid}{Score}){
								$orf_ref_scores{$orfid}{ORF}=$orfid;
								$orf_ref_scores{$orfid}{Score}=$max_score;
								$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
								$refs2orfs{$refgene}{$orfid}=1;
								$orfs2refs{$orfid}{$refgene}=1;
							}
						}
						else{
							$orf_ref_scores{$orfid}{ORF}=$orfid;
							$orf_ref_scores{$orfid}{Score}=$max_score;
							$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
							$refs2orfs{$true_gene}{$orfid}=1;
							$orfs2refs{$orfid}{$true_gene}=1;
						}

					}
				  }	

	my %coordinates;
	open my $coords, "<", $ARGV[5]; #Read in coordinates for ORFs
	while(<$coords>){
			chomp;
			my @tarray = split/\s+/;
			$coordinates{$tarray[0]}=$tarray[1];
	}
	
	for my $tempref (keys %refs2orfs){
		my %tempposition;
		my %orf_groups;
		###NEED TO add a way to group orfs into regions of the genome
		if(scalar (keys %{$refs2orfs{$tempref}}) > 1){
			for my $temporf (keys %{$refs2orfs{$tempref}}){
				my @tempcoord = split(/-/, $coordinates{$temporf});
				for (my $i=$tempcoord[0]; $i<=$tempcoord[1]; $i++){
						$tempposition{$i}++;
						push(@{$orf_groups{$i}}, $temporf);
				}
			}
		}
					
	my %temp_overlap;
	my $group_count=0;
	
	for my $ctemp (keys %orf_groups){
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
	for my $group_temp (keys %temp_overlap){
		for my $group_mem (keys %{$temp_overlap{$group_temp}}){
			for my $group_temp2 (keys %temp_overlap){
				if($group_temp2 eq $group_temp){
						next;
				}
				for my $group_mem2 (keys %{$temp_overlap{$group_temp}}){
					if ($group_mem2 eq $group_mem){
						$temp_overlap{$group_temp}=$temp_overlap{$group_temp2};
						delete $temp_overlap{$group_temp2};
					}
				}
			}
		}
	}
	for my $group_temp (keys %temp_overlap){
		my $temp_max=0;
		my $prev_temp;
		for my $group_mem (keys %{$temp_overlap{$group_temp}}){
		   
				if(!exists $orf_ref_scores{$group_mem}){
					next;
				}
			    
				if($orf_ref_scores{$group_mem}{Score} > $temp_max){
					$temp_max = $orf_ref_scores{$group_mem}{Score};
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

	open my $out, ">", $ARGV[4] . "_best_orfs_for_refs_SCORE.txt";
	for my $ref_gene (sort keys %orf_ref_scores){
		#if($orf_ref_scores{$ref_gene}{Real_Ref} =~ /petDXXX/){
		#	$orf_ref_scores{$ref_gene}{Real_Ref} =~ /petD(XXX.*?)/;	
		#	$orf_ref_scores{$ref_gene}{Real_Ref} = "petD_exon2$1";
		#}
		print $out "$orf_ref_scores{$ref_gene}{ORF}\t$orf_ref_scores{$ref_gene}{Real_Ref}\t$orf_ref_scores{$ref_gene}{Score}q\n";
	}
 }

else{
	my %best_orf_ref;
	my %orf_ref_scores;
	open my $file, "<", $ARGV[0]; #blast output
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
	#USAGE: 1: Sequence from plastome being annotated 2:Reference for sequence 3: Direction of sequence
	my $total=0;
	my $seq_score=0;
	for (my $i=0; $i<length($_[0])-2; $i++){
		my $test = substr($_[0],$i,3);
		$total++;
		
		if($_[1] =~ /$test/){
			$seq_score++;
		}
	}
	
	my $final_score= $seq_score/length($_[1]); ##($seq_score/length($_[1]));
	#print "$final_score\n";
	return($final_score);
}
###########################
sub annotation_score_rna{
	#USAGE: 1: Sequence from plastome being annotated 2:Reference for sequence 3: Direction of sequence
	my $total=0;
	my $seq_score=0;
	for (my $i=0; $i<length($_[1])-5; $i+=6){
		my $test = substr($_[1],$i,6);
		$total++;
		
		if($_[0] =~ /$test/){
			$seq_score++;
		}
	}
	
	my $final_score= $seq_score/$total;

	return($final_score);
}
