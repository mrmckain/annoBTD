#!/usr/bin/perl -w
use strict;
use diagnostics;

open(STDOUT, '>', $ARGV[4] . '_results.log') or die "Can't open log";
open(STDERR, '>', $ARGV[4] . '_results_error.log') or die "Can't open log";

my %blastgenes;
my %orfseqs;
my %orfgenes;
#Clean up blast files for annoBTD
# 1: blastfile 2: blastgeneseqs 3:query orfs/seqs 4:Feature type
#Read in ORFs sequences and reference tRNA, rRNA, and gene sequences
%blastgenes = &read_fasta($ARGV[1],%blastgenes); #These are the REFERENCE genes for annotation
%orfgenes = &read_fasta($ARGV[2],%orfgenes);
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
		#if($tarray[2]< 0.5){
		#		next;
		#}
		
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
					####Get sequence for current reference and orf######
						my $orf_seq = $orfgenes{$orfid};
						my $ref_seq = $blastgenes{$refgene};
						#if(length($orf_seq)/length($ref_seq) < 0.75){
						#	delete $blast_gene_hsps{$orfid}{$refgene};
						#	next;
						#}
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

						if($max_score < 0.75){
							#next;
							delete $blast_gene_hsps{$orfid}{$refgene};
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
							}
						}
						else{
							$orf_ref_scores{$orfid}{ORF}=$orfid;
							$orf_ref_scores{$orfid}{Score}=$max_score;
							$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
						}

					}
				  }	
				}
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

						if($max_score < 0.75){
							#next;
							delete $blast_gene_hsps{$orfid}{$refgene};
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
							}
						}
						else{
							$orf_ref_scores{$orfid}{ORF}=$orfid;
							$orf_ref_scores{$orfid}{Score}=$max_score;
							$orf_ref_scores{$orfid}{Real_Ref}=$refgene;
						}

					}
				  }	

	open my $out, ">", $ARGV[4] . "_best_orfs_for_refs_SCORE.txt";
	for my $ref_gene (sort keys %orf_ref_scores){
		print $out "$orf_ref_scores{$ref_gene}{ORF}\t$orf_ref_scores{$ref_gene}{Real_Ref}\n";
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
		my $score = &annotation_score_rna($unknown, $ref_seq);
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
	elsif($tsid){
		$seqid.=$_;
	}
}
#close $tempfasta;
return($seqid);
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
	for (my $i=0; $i<length($_[0])-2; $i+=3){
		my $test = substr($_[0],$i,3);
		$total++;
		
		if($_[1] =~ /$test/){
			$seq_score++;
		}
	}
	
	my $final_score= ($seq_score/$total)*(length($_[0])/length($_[1]));

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
