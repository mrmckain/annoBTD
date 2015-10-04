#!/usr/bin/perl -w
use strict;

my %tks;

open my $file, "<", $ARGV[0];
while(<$file>){
		chomp;

		my @tarray= split /\s+/;
		$tarray[1] =~ /(TK\d+)/;
		$tks{$1}=1;
		
}

open my $out1, ">", $ARGV[1] . "_taxagb.txt";
open my $out2, ">", $ARGV[1] . "_nottaxagb.txt";

open my $file1, "<", $ARGV[1];

while(<$file1>){
	
	chomp;
	/(TK\d+)/;
	if( exists $tks{$1}){
		print $out1 "$_\n";

	}	
	else{
		print $out2 "$_\n"
	}	
}