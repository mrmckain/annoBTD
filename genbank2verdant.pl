#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0]; #genbank annotation file
open my $out, ">", $ARGV[1] . "_genbank_annotation.vt";

my %final_annotation;
my $plastome;
while(<$file>){
		chomp;
		if(/LOCUS/){
			/(\d+)\sbp/;
			$plastome=$1;
		}
		if(/gene/ || /tRNA/ || /rRNA/){
			my @exons = $_ =~ /\d+/g;
			my $nextline = readline($file);
			unless($nextline =~ /gene/){
				next;
			}
			$nextline =~ /\"(.+)\"/;
			my $gene=$1;
			my $count=1;
			my $size = @exons;

			for (my $i=0; $i<$size; $i+=2){
				my $start = shift(@exons);
				my $end = shift(@exons);
				$final_annotation{$start}{$end}{"+"}=$gene;
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
					$final_annotation{1}{$start-1}{"+"}="Start\-" . $tids;
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
						$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "-" . $final_annotation{$start}{$end}{$dir};
						$prev = $final_annotation{$start}{$end}{$dir};
						$prev_end = $end;
					}
				}
				else{
					$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "-" . $final_annotation{$start}{$end}{$dir};
					$prev = $final_annotation{$start}{$end}{$dir};
					$prev_end = $end;
				}
			}
		}
	}
}

$final_annotation{$prev_end+1}{$plastome}{"+"}=$prev . "-End";


for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (keys %{$final_annotation{$start}}){
		for my $dir (keys %{$final_annotation{$start}{$end}}){
			print $out "$final_annotation{$start}{$end}{$dir}\t$start\t$end\t$dir\n";
		}
	}
}