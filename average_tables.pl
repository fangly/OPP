#! /usr/bin/env perl
 $dir=$ARGV[0] or die "Please provide the folder where the rarefied OTU tables are located\n";
 $outfile=$ARGV[1] or die "Please specify output file name\n";

 #glob
  @files=<$dir*>;
  #how many files with otu tables
  $numrarefactions=scalar @files;
  #open first file as a test to count how many samples
  open IN, "$files[0]";
  #skip header line
  <IN>;
  $header=<IN>;
  @samplenames=split/\t/,$header;
  shift(@samplenames);
  pop(@samplenames);
  $numsamples=(scalar @samplenames);
  %samplecounts;
  %otuuniquenumbers;
  close IN;
  foreach $file(@files){
open OTUTABLE, "$file";
while(<OTUTABLE>){
	#skip headers
	if(/#/){
		next;
	}	
              chomp;
	my @info=split/\t/;
	my $otunumber=shift @info;
	my $otuname=pop @info;
	$otuuniquenumbers{$otunumber}=$otuname;
	#Create entry for OTU if it doesn't exist
	if(!exists($samplecounts{$otunumber})){
		$samplecounts{$otunumber}=[@info];
		
	}
	#Add to existing entries
	else{
	 	my @values=@{$samplecounts{$otunumber}};
		if(scalar(@values)!=scalar(@info)||scalar@values!=$numsamples){die;}
		else{
			for($i=0;$i<$numsamples;$i++){
				$values[$i]+=$info[$i];
			}
			$samplecounts{$otunumber}=[@values];
		}
	}
}
close OTUTABLE;			
  }
  open OUT, ">$outfile";
  print OUT "OTUID\t";
  foreach $name(@samplenames){
print OUT "$name\t";
  }
  print OUT "Consensus Lineage\n";
  @lineages=keys %samplecounts;
  foreach $key(@lineages){
my @counts=@{$samplecounts{$key}};
print OUT "$key\t";
foreach $count(@counts){
	$average=int($count/$numrarefactions);
	print OUT "$average\t";
}
print OUT "$otuuniquenumbers{$key}\n";
  }
  close OUT;
