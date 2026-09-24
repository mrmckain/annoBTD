#!/usr/bin/perl -w
use strict;
#USAGE 1-plastome of guide species 2-verdant annotation file of guide species 3-Species to be annotated
my %elements;
my $plastome;

my $sid;
open my $pfile, "<", $ARGV[0]; #plastome seq
while(<$pfile>){
		chomp;
		if(/^>/){
			$sid=$_;
		}
		else{
			$plastome.=$_;
		}
}

open my $file, "<", $ARGV[1];  #verdant annotation file
while(<$file>){
		chomp;
		my @tarray = split /\s+/;
		if($tarray[0] !~ /\-/){
			if($tarray[0] !~ /\~/){
				if(/IRA/ || /IRB/ || /LSC/ || /SSC/ || /FULL/ || /intron/){
					next;
				}
				my $seq = substr($plastome, $tarray[1]-1, ($tarray[2]-$tarray[1]+1));
				if($tarray[3] eq "-"){
					$seq = reverse($seq);
					$seq =~ tr/ATCGatcg/TAGCtagc/;
				}
				$elements{$tarray[0]}=$seq;
			}
		}
		elsif($tarray[0] =~ /^trn\w+\-\w\w\w$/){
			my $seq = substr($plastome, $tarray[1]-1, ($tarray[2]-$tarray[1]+1));
			if($tarray[3] eq "-"){
				$seq = reverse($seq);
				$seq =~ tr/ATCGatcg/TAGCtagc/;
			}
			$elements{$tarray[0]}=$seq;
		}
		elsif($tarray[0] =~ /^trn\w+\-\w\w\w_exon\d$/){
			my $seq = substr($plastome, $tarray[1]-1, ($tarray[2]-$tarray[1]+1));
			if($tarray[3] eq "-"){
				$seq = reverse($seq);
				$seq =~ tr/ATCGatcg/TAGCtagc/;
			}
			$elements{$tarray[0]}=$seq;
		}
}

open my $outfile, ">>", $ARGV[2] . "_annotated_regions_fromverdant_genes.fsa";
open my $outfile2, ">>", $ARGV[2] . "_annotated_regions_fromverdant_trnas.fsa";
open my $outfile3, ">>", $ARGV[2] . "_annotated_regions_fromverdant_rrnas.fsa";
for my $gid (sort keys %elements){
		if($gid =~ /trn/){
				print $outfile2 ">$gid" . "XXX$ARGV[3]\n$elements{$gid}\n";
		}
		elsif($gid =~ /rrn/){
			print $outfile3 ">$gid" . "XXX$ARGV[3]\n$elements{$gid}\n";
		}
		else{
			print $outfile ">$gid" . "XXX$ARGV[3]\n$elements{$gid}\n";
		}
		
}
