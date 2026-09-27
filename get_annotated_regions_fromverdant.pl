#!/usr/bin/perl -w
use strict;
#USAGE 1-plastome of guide species 2-verdant annotation file of guide species 3-Species to be annotated
my %elements;
my %elem_order;      # preserve first-seen order for stable output
my %seen_seq;        # name -> [sequences], so IR copies do not duplicate a reference
my $elem_n = 0;

# Duplicate gene names are legitimate: isoacceptors (three distinct trnS loci) and
# IR copies both share a name. Keying by name alone dropped all but the last.
# Identical sequences (the IR case) collapse to one reference; genuinely different
# sequences under the same name each get their own, suffixed _2, _3, ...
sub add_element {
    my ($name, $seq) = @_;
    return unless defined $seq && length $seq;

    # IR copies of a gene are usually byte identical, but their two annotations
    # sometimes differ by a base or two, which exact matching treats as a second
    # distinct gene. That injects a near-duplicate reference for every such pair
    # and measurably degrades boundary choice, so near-identical counts as same.
    for my $prev (@{ $seen_seq{$name} || [] }) {
        return if &near_identical($prev, $seq);
    }
    push @{ $seen_seq{$name} }, $seq;

    my $n = scalar @{ $seen_seq{$name} };
    my $key = $n > 1 ? "${name}_${n}" : $name;
    $elements{$key} = $seq;
    $elem_order{$key} = $elem_n++;
}

# Same gene if one is contained in the other, or they are the same length and
# differ at no more than 5% of positions. Genuine isoacceptors (trnS-GCU vs
# trnS-UGA) are far more divergent than this and stay separate.
sub near_identical {
    my ($a, $b) = @_;
    return 1 if $a eq $b;
    my ($short, $long) = length($a) <= length($b) ? ($a, $b) : ($b, $a);
    return 0 if length($long) - length($short) > 6;
    return 1 if index($long, $short) >= 0;
    my $mm = 0;
    my $n  = length($short);
    return 0 unless $n;
    for (my $i = 0; $i < $n; $i++) {
        $mm++ if substr($short, $i, 1) ne substr($long, $i, 1);
        return 0 if $mm > $n * 0.05;
    }
    return 1;
}
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
				&add_element($tarray[0], $seq);
			}
		}
		elsif($tarray[0] =~ /^trn\w+\-\w\w\w$/){
			my $seq = substr($plastome, $tarray[1]-1, ($tarray[2]-$tarray[1]+1));
			if($tarray[3] eq "-"){
				$seq = reverse($seq);
				$seq =~ tr/ATCGatcg/TAGCtagc/;
			}
			&add_element($tarray[0], $seq);
		}
		elsif($tarray[0] =~ /^trn\w+\-\w\w\w_exon\d$/){
			my $seq = substr($plastome, $tarray[1]-1, ($tarray[2]-$tarray[1]+1));
			if($tarray[3] eq "-"){
				$seq = reverse($seq);
				$seq =~ tr/ATCGatcg/TAGCtagc/;
			}
			&add_element($tarray[0], $seq);
		}
}

open my $outfile, ">>", $ARGV[2] . "_annotated_regions_fromverdant_genes.fsa";
open my $outfile2, ">>", $ARGV[2] . "_annotated_regions_fromverdant_trnas.fsa";
open my $outfile3, ">>", $ARGV[2] . "_annotated_regions_fromverdant_rrnas.fsa";
for my $gid (sort { $elem_order{$a} <=> $elem_order{$b} } keys %elements){
		my $outname = $gid;
		$outname =~ s/_(\d+)$/_$1/;   # kept distinct in the FASTA id
		if($gid =~ /trn/){
				print $outfile2 ">$outname" . "XXX$ARGV[3]\n$elements{$gid}\n";
		}
		elsif($gid =~ /rrn/){
			print $outfile3 ">$outname" . "XXX$ARGV[3]\n$elements{$gid}\n";
		}
		else{
			print $outfile ">$outname" . "XXX$ARGV[3]\n$elements{$gid}\n";
		}
		
}
