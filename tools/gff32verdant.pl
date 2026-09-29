#!/usr/bin/perl -w
use strict;

#GFF3 to Verdant format parser
#Usage: 1-GFF3 file
open my $out, ">", $ARGV[0] . "_genbank_annotation.vt";
my %final_annotation;
my %temp_annotation;
my %direction;
my $plastome;
open my $file, "<", $ARGV[0];
while(<$file>){
		chomp;
		my @tarray = split/\s+/;
		if(@tarray <2){
				next;
		}
		if(/^\#/){
			next;
		}
		if($tarray[2] eq "region"){
			$plastome=$tarray[4];
		}
		
		if($tarray[2] eq "CDS"){
				$tarray[8] =~ /Parent\=(\w+)/;
				my $id = $1;
				$temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
				$direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];
		}
		if($tarray[2] eq "tRNA" ){
				$tarray[8] =~ /ID\=(\w+\-\w+)/;
				my $id = $1;
				$temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
				$direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];

		}
		if($tarray[2] eq "rRNA"){
			$tarray[8] =~ /ID\=(\w+)/;
				my $id = $1;
				$temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
				$direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];
		}
}

for my $geneid (keys %temp_annotation){
	my %temp_pos;
	my $posmin=10000000;
	my $posmax=-1;
	my $minmin=10000000;
	my $minmax=-1;
	
	for my $start (keys %{$temp_annotation{$geneid}}){
		for my $stop (keys %{$temp_annotation{$geneid}{$start}}){
			if(keys %{$direction{$geneid}{$temp_annotation{$geneid}{$start}{$stop}}} > 1){
				$temp_pos{$start}=$stop;
				if($temp_annotation{$geneid}{$start}{$stop} eq "+"){
					if($start>$posmax){
						$posmax=$start;
					}
					if($start<$posmin){
						$posmin=$start;
					}
					if($stop>$posmax){
						$posmax=$stop;
					}
					if($stop<$posmin){
						$posmin=$stop;
					}
				}
				if($temp_annotation{$geneid}{$start}{$stop} eq "-"){
					if($start>$minmax){
						$minmax=$start;
					}
					if($start<$minmin){
						$minmin=$start;
					}
					if($stop>$minmax){
						$minmax=$stop;
					}
					if($stop<$minmin){
						$minmin=$stop;
					}
				}			
			
		
		
				for my $tstart (keys %temp_pos){
					if($tstart != $posmin || $temp_pos{$tstart} != $posmax || $tstart != $minmin || $temp_pos{$tstart} != $minmax){

						if($temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}} eq "+"){
							my $exonnumber = keys %{$direction{$geneid}{"+"}};
							if($tstart == $posmin){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon1";
							}
							elsif($temp_pos{$tstart} == $posmax){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon" . $exonnumber;
							}
							else{
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon2";
							}
						
						}
						else{
							my $exonnumber = keys %{$direction{$geneid}{"-"}};
							if($tstart == $minmin){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon" . $exonnumber;
							}
							elsif($temp_pos{$tstart} == $minmax){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon1";
							}
							else{
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon2";
							}
						}
					
					}
				}
			}
			else{
		
				$final_annotation{$start}{$stop}{$temp_annotation{$geneid}{$start}{$stop}}=$geneid;
			}
		}		
		
	}
}
my $prev;
my $prev_end;
for my $start (sort {$a <=> $b} keys %final_annotation){
	
	for my $end (keys %{$final_annotation{$start}}){
		
		for my $dir (keys %{$final_annotation{$start}{$end}}){
			my $tids = $final_annotation{$start}{$end}{$dir};
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
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
						elsif ($pnum ==2 && $tnum ==1) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
						elsif ($pnum ==1 && $tnum ==2) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
						elsif ($pnum ==2 && $tnum ==3) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron2";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
					}
					else{
						$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "~" . $final_annotation{$start}{$end}{$dir};
						$prev = $final_annotation{$start}{$end}{$dir};
						$prev_end = $end;
					}
				}
				else{
					$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "~" . $final_annotation{$start}{$end}{$dir};
					$prev = $final_annotation{$start}{$end}{$dir};
					$prev_end = $end;
				}
			}
		}
	}
}

$final_annotation{$prev_end+1}{$plastome}{"+"}=$prev . "~End";


for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (keys %{$final_annotation{$start}}){
		for my $dir (keys %{$final_annotation{$start}{$end}}){
			print $out "$final_annotation{$start}{$end}{$dir}\t$start\t$end\t$dir\n";
		}
	}
}