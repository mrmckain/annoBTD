#!/usr/bin/perl -w
use strict;

my %blastgenes;
my %blastrnas;
my %orfs;
my %idorfs_f;
my %short_exons;
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
($torf_counter,$idorfs_f, $orf_pos) = &rna_blasts($torf_counter, $ARGV[7], \%idorfs_f, \%orf_pos);
%orf_pos = %$orf_pos;
%idorfs_f = %$idorfs_f;
($torf_counter,$idorfs_f,$orf_pos) = &rna_blasts($torf_counter, $ARGV[8], \%idorfs_f, \%orf_pos);
%orf_pos = %$orf_pos;
%idorfs_f = %$idorfs_f;

for my $gid (keys %blastgenes){
		if(length($blastgenes{$gid})<25){
			$short_exons{$gid}=1;
		}
}

for my $gid (keys %blastrnas){
		if(length($blastrnas{$gid})<26){
			$short_exons{$gid}=1;
		}
}

my %blast_overlaps;
open $file, "<", $ARGV[2]; #blast output
while(<$file>){
		chomp;
		if(/intron/){
			next;
		}
		my @tarray = split /\s+/;
		my $hit_truelen=length($blastgenes{$tarray[1]});
		my $hit_curlen=abs(($tarray[9]-$tarray[8]));
		if($hit_truelen > 4500){
			if($hit_curlen/$hit_truelen >= 0.05){
				$idorfs_f{$tarray[1]}{$tarray[0]}{"Start"}=$orf_pos{$tarray[0]}{"Start"};
				$idorfs_f{$tarray[1]}{$tarray[0]}{"End"}=$orf_pos{$tarray[0]}{"End"};
				if(exists $blast_overlaps{$tarray[1]}){
					if(exists $blast_overlaps{$tarray[1]}{$tarray[0]}){
						if(($hit_curlen/$hit_truelen) > $blast_overlaps{$tarray[1]}{$tarray[0]}){
							$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
						}
					}
					else{
						$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
					}
				}
				else{
					$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
				}
			}	

		}
		elsif($hit_truelen > 1400){
			if($hit_curlen/$hit_truelen >= 0.30){
				$idorfs_f{$tarray[1]}{$tarray[0]}{"Start"}=$orf_pos{$tarray[0]}{"Start"};
				$idorfs_f{$tarray[1]}{$tarray[0]}{"End"}=$orf_pos{$tarray[0]}{"End"};
				if(exists $blast_overlaps{$tarray[1]}){
					if(exists $blast_overlaps{$tarray[1]}{$tarray[0]}){
						if(($hit_curlen/$hit_truelen) > $blast_overlaps{$tarray[1]}{$tarray[0]}){
							$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
						}
					}
					else{
						$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
					}
				}
				else{
					$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
				}
			}	

		}
		elsif($hit_curlen/$hit_truelen >= 0.65){
				$idorfs_f{$tarray[1]}{$tarray[0]}{"Start"}=$orf_pos{$tarray[0]}{"Start"};
				$idorfs_f{$tarray[1]}{$tarray[0]}{"End"}=$orf_pos{$tarray[0]}{"End"};

				if(exists $blast_overlaps{$tarray[1]}){
					if(exists $blast_overlaps{$tarray[1]}{$tarray[0]}){
						if(($hit_curlen/$hit_truelen) > $blast_overlaps{$tarray[1]}{$tarray[0]}){
							$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
						}
					}
					else{
						$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
					}
				}
				else{
					$blast_overlaps{$tarray[1]}{$tarray[0]}= ($hit_curlen/$hit_truelen);
				}
		}
}
close $file;





%idorfs_f = &hit_cleaner(%idorfs_f);


my $plastome;
open $file, "<", $ARGV[4]; #full plastome being annotated
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
my $orf_counter=0;
for my $shid (keys %short_exons){

		my $tseq;
		my ($gene, $exonnum);
		if($shid =~ /trn/){
			$tseq=$blastrnas{$shid};
			$shid =~ /(trn\w+-\w\w\w)_exon(\d+)/;
			$gene=$1;
			$exonnum=$2;
		}
		else{
			$tseq = $blastgenes{$shid};
			$shid =~ /(\w+)_exon(\d+)/;
			$gene = $1;
			$exonnum = $2;
		}
		
		
		
		if($exonnum == 1){
				my $intron_name = $gene . "_intron1";

				#if(exists $blastgenes{$intron_name}){
					$exonnum++;
					my $next_orf = $gene . "_exon" . $exonnum;
					
					if(scalar keys %{$idorfs_f{$next_orf}} == 1){
							for my $norf (keys %{$idorfs_f{$next_orf}}){
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
									my $distance = 1000000;
									my $good_start;
									if(scalar @tarray > 1){ #checking for multiple matches for closest proximity match
										for my $pot_exon (@tarray){
											if(abs($pot_exon-$close_end) < $distance){
												$good_start = $pot_exon;
											}

										}
									}
									else{
										
											$good_start = shift(@tarray);
										}
									
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
									my $distance = 1000000;
									my $good_start;
									if(scalar @tarray > 1){ #checking for multiple matches for closest proximity match
										for my $pot_exon (@tarray){
											if(abs($close_start-$pot_exon) < $distance){
												$good_start = $pot_exon;
												$distance = abs($close_start-$pot_exon);
											}

										}
									}
									else{
										$good_start = shift(@tarray);
									}
									
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
}
my %gene_features;
for my $gid (keys %idorfs_f){
	if($gid =~ /exon/){
			$gid =~ /(.+)_exon(\d+)/;

			$gene_features{$1}{$gid}=1;
	}
	else{
			
			$gene_features{$gid}{$gid}=1;
	}
}
#Identify start/stop for genes (not tRNA, rRNA, and short exons)
my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');
my %final_annotation;
for my $gid (sort keys %idorfs_f){
	print "$gid\n";
	for my $eid (keys %{$idorfs_f{$gid}}){
		if(exists $short_exons{$gid} || $gid =~ /rrn/ || $gid =~ /trn/){
			my $dir = "+";
			if($eid =~ /\-r/){
				$dir = "-";
			}
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}}{$idorfs_f{$gid}{$eid}{"End"}}{$dir}=$gid;
		}
		else{
			my $blastsize = length($blastgenes{$gid});
			my $blmod = $blastsize%3;
			if($gid =~ /exon/){
				$gid =~ /(.+)_exon(\d+)/;

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

my %temp_store;
for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (keys %{$final_annotation{$start}}){
		if(%temp_store){
		for my $temp_start (sort {$a <=> $b} keys %temp_store){
			for my $temp_end (sort {$a <=> $b} keys %{$temp_store{$temp_start}}){
				if($start == $temp_start && $end == $temp_end){
					for my $dir (keys %{$final_annotation{$start}{$end}}){
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
					for my $dir (keys %{$final_annotation{$start}{$end}}){
						$temp_store{$start}{$end}=$final_annotation{$start}{$end}{$dir};
					}
				}
				elsif($start >= $temp_start && $end < $temp_end){
					delete $final_annotation{$start}{$end};
				}
				elsif($temp_start >= $start && $temp_end < $end){
					delete $final_annotation{$temp_start}{$temp_end};
					for my $dir (keys %{$final_annotation{$start}{$end}}){
						$temp_store{$start}{$end}=$final_annotation{$start}{$end}{$dir};
					}
				}
				elsif($start < $temp_start && $temp_end > $end){
					if(abs($end-$start) < abs($temp_end-$temp_start)){
						delete $final_annotation{$start}{$end};
					}
					if(abs($temp_end-$temp_start) < abs($end-$start)){
						delete $final_annotation{$temp_start}{$temp_end};
						for my $dir (keys %{$final_annotation{$start}{$end}}){
							$temp_store{$start}{$end}=$final_annotation{$start}{$end}{$dir};
						}
					}
				} 
				elsif($temp_start > $start && $temp_end > $end){
					if(abs($end-$start) < abs($temp_end-$temp_start)){
						delete $final_annotation{$start}{$end};
					}
					if(abs($temp_end-$temp_start) < abs($end-$start)){
						delete $final_annotation{$temp_start}{$temp_end};
						for my $dir (keys %{$final_annotation{$start}{$end}}){
							$temp_store{$start}{$end}=$final_annotation{$start}{$end}{$dir};
						}
					}
				}
				else{
					for my $dir (keys %{$final_annotation{$start}{$end}}){
							$temp_store{$start}{$end}=$final_annotation{$start}{$end}{$dir};
						}
				}
			}
		}
		}
		else{
			for my $dir (keys %{$final_annotation{$start}{$end}}){
				$temp_store{$start}{$end}=$final_annotation{$start}{$end}{$dir};
			}	
		}
	}
}


my $prev;
my $prev_end;
for my $start (sort {$a <=> $b} keys %final_annotation){
	
	for my $end (sort {$a <=> $b} keys %{$final_annotation{$start}}){
		
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
						for my $toid (keys %temporfs){
							if($temporfs{$toid}{"Start"} > $cur_start && $temporfs{$toid}{"End"} < $cur_end){
								if($blast_overlaps{$bid}{$toid} > $blast_overlaps{$bid}{$oid}){
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
								if($blast_overlaps{$bid}{$toid} > $blast_overlaps{$bid}{$oid}){
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
								if($blast_overlaps{$bid}{$toid} >= $blast_overlaps{$bid}{$oid}){
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
								if($blast_overlaps{$bid}{$toid} > $blast_overlaps{$bid}{$oid}){
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
sub rna_blasts{

my ($counter,$file, $idorfs, $torf_pos) = @_;
my %idorfs = %$idorfs;
my %torf_pos = %$torf_pos;

open my $tfile, "<", $file; #blast output
while(<$tfile>){
		chomp;
		my $cur_orfname = "orfRNA" . $counter;
		my @tarray = split /\s+/;
		my $hit_truelen=length($blastrnas{$tarray[1]});
		if (($tarray[9]-$tarray[8])<0){
			$cur_orfname = $cur_orfname . "-r";
		}
		my $hit_curlen=abs(($tarray[9]-$tarray[8]))+1;
		if($hit_curlen/$hit_truelen >= 0.98){
				$idorfs{$tarray[1]}{$cur_orfname}{"Start"}=$tarray[6];
				$idorfs{$tarray[1]}{$cur_orfname}{"End"}=$tarray[7];
				$torf_pos{$cur_orfname}{"Start"}=$tarray[6];
				$torf_pos{$cur_orfname}{"End"}=$tarray[7];
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
		if($_[1] =~ /$test/){
			$forward++
		}
		if($_[2] =~ /$test/){
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
		if($_[1] =~ /$test/){
			$forward_0++
		}
		if($_[2] =~ /$test/){
			$reverse_0++;
		}
		if($_[3] =~ /$test/){
			$forward_1++
		}
		if($_[4] =~ /$test/){
			$reverse_1++;
		}
		if($_[5] =~ /$test/){
			$forward_2++
		}
		if($_[6] =~ /$test/){
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
sub match_start{
	my $altstart=-1;
	my $seq = $_[0];
	for (my $i=0; $i<($_[1]*3);$i++){#second parameter is $nearstart
		if(index($seq,"ATG",$i)>=0){
			if(index($seq,"ATG",$i) < ($_[1]*3)){
				$altstart =index($seq,"ATG",$i);
				$i=$altstart;
			}

		}
	}
	return($altstart);
}
###########################
sub alt_start{
	my $matchstart=-1;
	my $seq = $_[0];
	for (my $i=0; $i<($_[1]*3+1);$i++){
		if(index($seq,substr($_[2],0,3),$i)>=0){#third parameter is gene_seq
			if(index($seq,substr($_[2],0,3),$i) < ($_[1]*3+1)){
				$matchstart =index($seq,substr($_[2],0,3),$i);
				$i=$matchstart;
			}
		}
	}
	return($matchstart);
}
###########################
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
		my $nearstart = index($orf_prot, substr($gene_prot,1,3));
		my $fix_nearstart=0;
		if($nearstart == -1){
			for (my $i =2; $nearstart<0; $i++){
				$nearstart = index($orf_prot, substr($gene_prot,$i,3));
				$fix_nearstart = $i;
			}
		}
		$nearstart = $nearstart-$fix_nearstart;
		my $nearend;
		my $fix_nearend=0;
		$nearend = index($orf_prot,substr(substr($gene_prot,-4),0,3), $nearstart);
		while($nearend == -1){
			for (my $i =2; $nearend<0; $i++){
				$nearend = index($orf_prot, substr(substr($gene_prot,(-3-$i)),0,3));
				$fix_nearend = $i;
			}
		}
=item		for (my $i=0; $i<(length($orf_prot)-1);$i++){
		if(index($orf_prot,substr(substr($gene_prot,-4),0,3),$i)>=0){#third parameter is gene_seq
			if(index($orf_prot,substr(substr($gene_prot,-4),0,3),$i) < (length($orf_prot))){
				$nearend =index($orf_prot,substr(substr($gene_prot,-4),0,3),$i);
				$i=$nearend;
			}
			
			}
		}
=cut	
		$nearend = $nearend + $fix_nearend;	
		my $end_match;
		if(substr($gene_prot, -1) eq "_"){
			$end_match = index($orf_prot,"_", $nearstart);

			$end_match = length($orf_seq)-($end_match*3 + 2)-1;
		}
		else{
			$end_match = &match_end($orf_prot,$gene_prot,$nearend);
		}
		my $start_match = &match_start($orf_seq,$nearstart,$gene_seq);
		my $start_alt = &alt_start($orf_seq,$nearstart,$gene_seq);
		if(substr($gene_seq,0,3) eq substr($orf_seq,$start_match,3)){
				$start_alt = -1;
		}


		if($start_match == -1 && $start_alt == -1){
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+($nearstart*3)-3+1}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1-$mod_shift}{"+"}=$gid;		}
		elsif($start_match>=$start_alt){
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_match-$for_frame+1}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1-$mod_shift}{"+"}=$gid;
		}
		else{
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_alt-$for_frame+1}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1-$mod_shift}{"+"}=$gid;
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
		my $nearstart = index($rev_orf_prot, substr($gene_prot,1,3));
		my $fix_nearstart=0;
		if($nearstart == -1){
			for (my $i =2; $nearstart<0; $i++){
				$nearstart = index($rev_orf_prot, substr($gene_prot,$i,3));
				$fix_nearstart = $i;
			}
		}
		$nearstart = $nearstart-$fix_nearstart;
		my $nearend;
		my $fix_nearend=0;
		$nearend = index($rev_orf_prot,substr(substr($gene_prot,-4),0,3),$nearstart);
		while($nearend == -1){
			for (my $i =2; $nearend<0; $i++){
				$nearend = index($rev_orf_prot, substr(substr($gene_prot,(-3-$i)),0,3));
				$fix_nearend = $i;
			}
		}
		$nearend = $nearend + $fix_nearend;	
		
		my $end_match;
		if(substr($gene_prot, -1) eq "_"){
			$end_match = index($rev_orf_prot,"_", $nearstart);
			$end_match = length($rev_orf_seq)-($end_match*3 + 2)-1;
			$mod_shift=0;
		}
		else{
			$end_match = &match_end($rev_orf_prot,$gene_prot, $nearend);
		}
		my $start_match = &match_start($rev_orf_seq,$nearstart,$gene_seq);
		my $start_alt = &alt_start($rev_orf_seq,$nearstart,$gene_seq);
		if(substr($gene_seq,0,3) eq substr($rev_orf_seq,$start_match,3)){
				$start_alt = -1;
		}
		if($start_match == -1 && $start_alt == -1){
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-($nearstart*3)+3+1}{"-"}=$gid;
		}
		elsif($start_match>$start_alt){
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-$start_match-$rev_frame+1}{"-"}=$gid;
		}
		else{
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-$start_alt-$rev_frame+1}{"-"}=$gid;
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
		my $nearend;
		my $fix_nearend=0;
		$nearend = index($orf_prot,substr(substr($gene_prot,-4),0,3));
		while($nearend == -1){
			for (my $i =2; $nearend<0; $i++){
				$nearend = index($orf_prot, substr(substr($gene_prot,(-3-$i)),0,3));
				$fix_nearend = $i;
			}
		}
		$nearend = $nearend + $fix_nearend;
		my $end_match;
		if(substr($gene_prot, -1) eq "_"){
			$end_match = index($orf_prot,"_", $nearend);

			$end_match = length($orf_seq)-($end_match*3 + 2)-1;
			$end_match = length($orf_seq) - $end_match;
		}
		else{
			$end_match = &match_end($orf_prot,$gene_prot,$nearend);
		}
		my $nearstart = index($orf_prot, substr($gene_prot,1,3));
		my $fix_nearstart=0;
		if($nearstart == -1){
			for (my $i =2; $nearstart<0; $i++){
				$nearstart = index($orf_prot, substr($gene_prot,$i,3));
				$fix_nearstart = $i;
			}
		}
		$nearstart = $nearstart-$fix_nearstart;
		#my $start_match = &match_start($orf_seq);
		my $start_alt = &alt_start($orf_seq,$nearstart,$gene_seq);
		#my $end_match = &match_end($orf_seq);

		
		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_alt-$frame_restore+1+$for_frame}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+1}{"+"}=$gid;
		
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
		my $nearstart = index($rev_orf_prot, substr($gene_prot,1,3));
		my $fix_nearstart=0;
		if($nearstart == -1){
			for (my $i =2; $nearstart<0; $i++){
				$nearstart = index($rev_orf_prot, substr($gene_prot,$i,3));
				$fix_nearstart = $i;
			}
		}
		$nearstart = $nearstart-$fix_nearstart;
		my $nearend;
		my $fix_nearend=0;
		$nearend = index($rev_orf_prot,substr(substr($gene_prot,-4),0,3));
		while($nearend == -1){
			for (my $i =2; $nearend<0; $i++){
				$nearend = index($rev_orf_prot, substr(substr($gene_prot,(-3-$i)),0,3));
				$fix_nearend = $i;
			}
		}
		$nearend = $nearend + $fix_nearend;
		my $end_match;
		if(substr($gene_prot, -1) eq "_"){
			$end_match = index($rev_orf_prot,"_",$nearstart);
			$end_match = length($rev_orf_seq)-($end_match*3 + 2)-1;
		}
		else{
			$end_match = &match_end($rev_orf_prot,$gene_prot, $nearstart);
		}
		my $start_alt = &alt_start($rev_orf_seq,$nearstart,$gene_seq);
		#my $end_match = &match_end($rev_orf_seq);

		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+1+$end_match}{$idorfs_f{$gid}{$eid}{"End"}-$start_alt+$frame_restore+1-$rev_frame}{"-"}=$gid;
		
	}

}
###########################
sub exon_mods_mid{
	my $frame_restore = $_[0];
	my $gid = $_[1]; #pass $gid to function
	my $eid = $_[2]; #pass $eid to function
	my $gene_seq=$blastgenes{$gid};
	my @frames;
	my $stop_pos=0;
	my $frame=0;
	for (my $j=0; $j <= $frame_restore; $j++){
		my @cur_frame = &translate(substr($gene_seq, $j));
		if(index($cur_frame[0], "_") > $stop_pos || $cur_frame[0] !~ /_/){
			$frame = $j;
			$stop_pos = index($cur_frame[0], "_");
			
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
		my $nearstart = index($orf_prot, substr($gene_prot,1,3));
		my $nearend;
		my $fix_nearend=0;
		$nearend = index($orf_prot,substr(substr($gene_prot,-4),0,3));
		while($nearend == -1){
			for (my $i =2; $nearend<0; $i++){
				$nearend = index($orf_prot, substr(substr($gene_prot,(-3-$i)),0,3));
				$fix_nearend = $i;
			}
		}
		$nearend = $nearend + $fix_nearend;
		my $start_alt = &alt_start($orf_seq,$nearstart,$gene_seq);
		my $end_match = &match_end($orf_prot,$gene_prot,$nearend);
		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_alt-$frame+1}{$idorfs_f{$gid}{$eid}{"End"}+-$end_remove+$end_remove+1}{"+"}=$gid;
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
		my $nearstart = index($rev_orf_prot, substr($gene_prot,1,3));

		my $nearend;
		my $fix_nearend=0;
		$nearend = index($rev_orf_prot,substr(substr($gene_prot,-4),0,3));
		while($nearend == -1){
			for (my $i =2; $nearend<0; $i++){
				$nearend = index($rev_orf_prot, substr(substr($gene_prot,(-3-$i)),0,3));
				$fix_nearend = $i;
			}
		}
		$nearend = $nearend + $fix_nearend;
		my $start_alt = &alt_start($rev_orf_seq,$nearstart,$gene_seq);
		my $end_match = &match_end($rev_orf_prot,$gene_prot,$nearend);
		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}-$end_remove+$end_match+1}{$idorfs_f{$gid}{$eid}{"End"}-$start_alt+$frame+1}{"-"}=$gid;
	}
}
###########################
open my $out, ">", "testing_blast_orfid.txt";
for my $bid (sort keys %idorfs_f){
	for my $oid (sort keys %{$idorfs_f{$bid}}){
		print $out "$bid\t$oid\n";
	}
}

open my $outfile2, ">", $ARGV[9] . "_VERDANT_cleaned_annotation.txt";

for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (keys %{$final_annotation{$start}}){
		for my $dir (keys %{$final_annotation{$start}{$end}}){
			print $outfile2 "$final_annotation{$start}{$end}{$dir}\t$start\t$end\t$dir\n";
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
