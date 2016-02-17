#!/usr/bin/perl

@file=`grep Energy lowest`;
if (scalar(@ARGV)!=2){&usage};
foreach $line (@file){
	@words = split(/=|\s+/,$line);
	$i = $words[3];
	$energy = $words[5];
	$line=`TubeHBond lowest$i.1.pdb $ARGV[0] $ARGV[1]`;
	print "$i $energy $line";
}

sub usage{
	die "Usage: count_hb.pl <int ringsize> <int num_rings>\n";
}
