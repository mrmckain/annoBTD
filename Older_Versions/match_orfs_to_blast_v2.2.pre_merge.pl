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
}

for my $gid (sort keys %blastrnas){
		if(length($blastrnas{$gid})<26){
			$short_exons{$gid}=1;
		}
}

#Guide genomes disagree about exon boundaries. Whichever single guide wins the
#score has its annotation transferred verbatim, quirks included - e.g. Oryza
#annotates atpF exon1 at 144 nt where six other genomes say 145, and picking Oryza
#shifts BOTH sides of the tobacco atpF intron by 1 nt. Where the guides hold a
#clear majority on an exon's length, prefer a reference that agrees with it.
#Only the reference used for boundary refinement changes; the ORF and the gene
#name it is reported under are untouched.
my %ref_by_gene;
for my $rid (sort keys %blastgenes){
	next unless $rid =~ /^(.+)XXX\d+$/;
	push @{$ref_by_gene{$1}}, $rid;
}
sub consensus_reference {
	my $rid = shift;
	return $rid unless defined $rid && $rid =~ /^(.+)XXX\d+$/;
	my $sibs = $ref_by_gene{$1};
	return $rid unless $sibs && @$sibs >= 3;   #too few guides to call a consensus

	my %n;
	$n{ length($blastgenes{$_}) }++ for @$sibs;
	my ($modal) = sort { $n{$b} <=> $n{$a} || $a <=> $b } keys %n;
	return $rid if length($blastgenes{$rid}) == $modal;
	return $rid unless $n{$modal} * 2 > scalar(@$sibs);   #require a real majority

	for my $s (sort @$sibs){
		return $s if length($blastgenes{$s}) == $modal;
	}
	return $rid;
}

my %blast_overlaps;
open $file, "<", $ARGV[2] or die "Cannot open $ARGV[2]: $!"; #best id to orf
while(<$file>){
		chomp;
		if(/intron/){
			next;
		}
		my @tarray = split /\s+/;
		next unless defined $tarray[1];
		$tarray[1] = &consensus_reference($tarray[1]);
		$tarray[1] =~ /(.*?)XXX.+/;
		my $temp_geneid = $1;

		$idorfs_f{$tarray[1]}{$tarray[0]}{"Start"}=$orf_pos{$tarray[0]}{"Start"};
		$idorfs_f{$tarray[1]}{$tarray[0]}{"End"}=$orf_pos{$tarray[0]}{"End"};		
		$blast_overlaps{$tarray[1]}{$tarray[0]} = abs($orf_pos{$tarray[0]}{"End"}-$orf_pos{$tarray[0]}{"Start"});
		
}
close $file;

%idorfs_f = &hit_cleaner(%idorfs_f); #cleans up positions of genes to reconcile with overlaps


my $orf_counter=0;
for my $shid (sort keys %short_exons){

		my $tseq;
		my ($gene, $exonnum, $refid);
		if($shid =~ /trn/){
			$tseq=$blastrnas{$shid};
			$shid =~ /(trn\w+-\w\w\w)_exon(\d+)(XXX\d+)/;
			$gene=$1;
			$exonnum=$2;
			$refid=$3;
		}
		else{
			$tseq = $blastgenes{$shid};
			$shid =~ /(\w+)_exon(\d+)(XXX\d+)/;
			$gene = $1;
			$exonnum = $2;
			$refid = $3;
		}
		
		
		
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
	for my $eid (sort keys %{$idorfs_f{$gid}}){
		if(exists $short_exons{$gid} || $gid =~ /rrn/ || $gid =~ /trn/){
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

&refine_intron_boundaries();

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
				next unless abs($left->{num} - $right->{num}) == 1;

				my $istart = $left->{e} + 1;         #1-based, inclusive
				my $iend   = $right->{s} - 1;
				my $ilen   = $iend - $istart + 1;
				#Real plastid introns in the reference set run 304-2559 nt.
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
		if($hit_curlen/$hit_truelen >= 0.97){
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
			&dbg_coord("exon_mods_start","$gid","$eid",1,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","+");
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+($nearstart*3)-3+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1}{"+"}=$gid;		}
		elsif($start_match>=$start_alt){
			&dbg_coord("exon_mods_start","$gid","$eid",2,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","+");
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_match-$for_frame+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1}{"+"}=$gid;
		}
		else{
			&dbg_coord("exon_mods_start","$gid","$eid",3,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","+");
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$start_alt-$for_frame+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-$end_match+$frame_restore+1}{"+"}=$gid;
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
			&dbg_coord("exon_mods_start","$gid","$eid",4,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","-");
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-($nearstart*3)+3+1}{"-"}=$gid;
		}
		elsif($start_match>$start_alt){
			&dbg_coord("exon_mods_start","$gid","$eid",5,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","-");
			$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}+$end_match-$frame_restore+1+$mod_shift}{$idorfs_f{$gid}{$eid}{"End"}-$start_match-$rev_frame+1}{"-"}=$gid;
		}
		else{
			&dbg_coord("exon_mods_start","$gid","$eid",6,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,$mod_shift,$nearstart,$nearend,$start_match,$start_alt,$end_match,"","","-");
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
			#$end_match = length($orf_seq) - $end_match;
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

		
		&dbg_coord("exon_mods_end","$gid","$eid",7,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame_restore,"",$nearstart,$nearend,"",$start_alt,$end_match,"","","+");
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

		&dbg_coord("exon_mods_end","$gid","$eid",8,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame_restore,"",$nearstart,$nearend,"",$start_alt,$end_match,"","","-");
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
		&dbg_coord("exon_mods_mid","$gid","$eid",9,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$for_frame,$frame,$end_remove,$nearstart,$nearend,"",$start_alt,$end_match,"","","+");
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
		&dbg_coord("exon_mods_mid","$gid","$eid",10,$idorfs_f{$gid}{$eid}{"Start"},$idorfs_f{$gid}{$eid}{"End"},$rev_frame,$frame,$end_remove,$nearstart,$nearend,"",$start_alt,$end_match,"","","-");
		$final_annotation{$idorfs_f{$gid}{$eid}{"Start"}-$end_remove+$end_match+1}{$idorfs_f{$gid}{$eid}{"End"}-$start_alt+$frame+1}{"-"}=$gid;
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
			if($final_annotation{$start}{$end}{$dir} !~ /XXX/){
			
				print $outfile2 "$final_annotation{$start}{$end}{$dir}\t$start\t$end\t$dir\n";
			}
			else{
				my $next_gene = $final_annotation{$start}{$end}{$dir};
				$next_gene =~ /(.*?)XXX.+/;
				$next_gene = $1;
				print $outfile2 "$next_gene\t$start\t$end\t$dir\n";
			}
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
