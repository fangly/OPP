#! /usr/bin/perl
$infile=$ARGV[0] or die "No OTU table input";
$outfile=$ARGV[1] or die "Please specify output file name";
open IN, "$infile";
open OUT, ">$outfile";

my $currentphylum="No";
@colors=("black","red","blue","green","grey");
$colorindex=0;
$colormax=scalar @colors;


while(<IN>){
	if(/No blast/){
		print OUT "$colors[$colorindex]\n";
	}
	if(/p__([A-Za-z]*)/){
		$phylum=$1;
		if($currentphylum ne $phylum){	
			$colorindex++;
			if($colorindex>=$colormax){
				$colorindex=0;
			}
			$currentphylum=$phylum;
		}
		print OUT "$colors[$colorindex]\n";
	}
}
