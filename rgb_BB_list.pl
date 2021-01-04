#!/usr/bin/perl


########################################################
#Example: perl rgb_BB_list.pl 2000 10000 1             #
########################################################

use strict;
use warnings;

my $begin = 2000.0 ;
my $end = 10000.0;
my $increment = 1.0;

my $T_flag = 0;

my @args = @ARGV;
if (@ARGV) {
	($begin) = shift @args;
  ($end) = (@args > 0) ? shift @args: "none";
	($increment) = (@args > 0) ? shift @args: "none";
	($T_flag) = (@args > 0) ? shift @args: "none";
}

if ($T_flag eq "-t"){
	$T_flag = 1;
}

for(my $i = $begin; $i <= $end; $i+=$increment ) {
	
	if ($T_flag){
		print "$i ";
	}
	
	my $cmd = "./StarColours -b $i";
	system($cmd);

}
