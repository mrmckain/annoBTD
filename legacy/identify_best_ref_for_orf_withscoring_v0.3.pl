#!/usr/bin/perl -w
use strict;
use diagnostics;

open(STDOUT, '>', $ARGV[4] . '_results.log') or die "Can't open log";
open(STDERR, '>', $ARGV[4] . '_results_error.log') or die "Can't open log";

#system("pmpath ExtUtils::MakeMaker");
#system("pmpath Acme::AandB");

my %blastgenes;
my %orfseqs;
#Clean up blast files for annoBTD
# 1: blastfile 2: blastgeneseqs 3:query orfs/seqs 4:Feature type
#Read in ORFs sequences and reference tRNA, rRNA, and gene sequences
%blastgenes = &read_fasta($ARGV[1],%blastgenes);
%orfseqs = &read_fasta($ARGV[2], %orfseqs);
my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');
my $feature=$ARGV[3];
#%blastrnas = &read_fasta($ARGV[5],%blastrnas);
#%blastrnas = &read_fasta($ARGV[6],%blastrnas);

my %blast_gene_hsps;
if($feature eq "protein"){
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
close $file;
my %best_orf_ref;
my %orf_ref_scores;

for my $ref_gene (keys %blast_gene_hsps){
	my $best_match;
	my $best_score=0;

	for my $orfid (keys %{$blast_gene_hsps{$ref_gene}}){
		my $orf_seq = $orfseqs{$ref_gene};
		if($ref_gene =~ /\-r/){
			$orf_seq = reverse($orf_seq);
			$orf_seq =~ tr/ATCGatcg/TAGCtagc/;
		}

		my $trans_orf = translate($orf_seq);

		my $ref_seq = $blastgenes{$orfid};
		my $longest;
		my $longest_length=0;
		if(length($blastgenes{$orfid})%3 > 0){
			my $frame1_ref = translate($blastgenes{$orfid});
			my $frame2_ref = translate(substr($blastgenes{$orfid}, 1));
			my $frame3_ref = translate(substr($blastgenes{$orfid}, 2));
			
			my $temp_longest_length=0;

			for (my $i=0; $i<=length($frame1_ref)-1; $i++){
				if(substr($frame1_ref,$i,1) ne "X" || substr($frame1_ref,$i,1) ne "_"){
					$temp_longest_length++;
				}
				else{
					if($temp_longest_length > $longest_length){
						$ref_seq = $frame1_ref;
						$longest_length=$temp_longest_length;
					}
				}
			}
			for (my $i=0; $i<=length($frame2_ref)-1; $i++){
				if(substr($frame2_ref,$i,1) ne "X" || substr($frame2_ref,$i,1) ne "_"){
					$temp_longest_length++;
				}
				else{
					if($temp_longest_length > $longest_length){
						$ref_seq = $frame2_ref;
						$longest_length=$temp_longest_length;
					}
				}
			}
			for (my $i=0; $i<=length($frame3_ref)-1; $i++){
				if(substr($frame3_ref,$i,1) ne "X" || substr($frame3_ref,$i,1) ne "_"){
					$temp_longest_length++;
				}
				else{
					if($temp_longest_length > $longest_length){
						$ref_seq = $frame3_ref;
						$longest_length=$temp_longest_length;
					}
				}
			}
		}
		my $trans_ref = translate($ref_seq);
		my $score = annotation_score($trans_orf, $trans_ref);
		if($score < 0.85){
				next;
				delete $blast_gene_hsps{$ref_gene}{$orfid};
				#delete $orfseqs{$orfid};
				
		}
		$orfid =~ /(.*?)XXX.+/;
		my $true_gene = $1;
		####NEED TO KEEP THE RIGHT REFERENCE ID BUT LOOK AT THEM ALL SEPARATELY. 
		if(exists $orf_ref_scores{$ref_gene}){

			if ($score > $orf_ref_scores{$ref_gene}{Score}){
				$orf_ref_scores{$ref_gene}{ORF}=$ref_gene;
				$orf_ref_scores{$ref_gene}{Score}=$score;
				$orf_ref_scores{$ref_gene}{Real_Ref}=$orfid;
			}
		}
		else{
			$orf_ref_scores{$ref_gene}{ORF}=$ref_gene;
			$orf_ref_scores{$ref_gene}{Score}=$score;
			$orf_ref_scores{$ref_gene}{Real_Ref}=$orfid;
		}

	}

=item	
	for my $orfid (keys %{$blast_gene_hsps{$ref_gene}}){
		my $total = 0;
		my $overlap = 0;
		for my $pos (keys %{$blast_gene_hsps{$ref_gene}{$orfid}}){
			if($blast_gene_hsps{$ref_gene}{$orfid}{$pos} == 1){
				$total++;
				$overlap++;
			}
			else{
				$total++;
			}
		}
		my $temp_score = $overlap/$total;
		if($temp_score > $best_score){
			$best_match = $orfid;
			$best_score = $temp_score;
		}
	}
	$best_orf_ref{$ref_gene}=$best_match;
=cut	
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
