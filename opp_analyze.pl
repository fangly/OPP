#! /usr/bin/env perl
###############################################################################
#
#    opp_make_results.pl
#    
#    Normalise and complete the QIIME pieline
#
#    Copyright (C) 2011 Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#CPAN modules
use File::Basename;

#locally-written modules
#load the pretty names for the fields
use FindBin qw($Bin);
use lib "$Bin";
use OPPConfig;

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
# GLOBALS
my $global_working_dir = "~";
my $global_conf_base = basename($options->{'config'});
my %global_samp_ID_list = ();
my $global_norm = 0;
my $global_rare_M = 10;
my $global_rare_X = 100;
my $global_rare_S = 10;
my $global_rare_N = 10;

print "Checking if all the config checks out...\t\t";

# Get the conf base
my @cb_1 = split /_/, $global_conf_base;
my @cb_2 = split /\./, $cb_1[1];
$global_conf_base = $cb_2[0];

# GET THE WORKING DIR
# working dir is the dir of the config file
# get the present dir
my $pwd = `pwd`;
chomp $pwd;
$global_working_dir = dirname($pwd."/".$options->{'config'});

# Override values from config
parse_config_results();

# Start the results pipeline!
print "All good!\n";

if (not -e "$global_working_dir/processing") {
   mkdir "$global_working_dir/processing" or die "Error: Could not write output directory $global_working_dir/processing\n$!\n";
}

if (not -e "$global_working_dir/results/") {
   mkdir "$global_working_dir/results/" or die "Error: Could not write output directory $global_working_dir/results\n$!\n";
}

# NORMALISE
`cp $global_working_dir/QA/denoised_acacia/acacia_out__all_tags.seqOut $global_working_dir/processing/normalised.fa`; 

#### Don't chdir
chdir "$global_working_dir/processing" or die "Error: Could not cd to directory $global_working_dir/processing\n$!\n";

print "Clustering OTUs...\n";
`pick_otus.py -i normalised.fa`;

print "Picking OTU representative sequences...\n";
`pick_rep_set.py -i uclust_picked_otus/normalised_otus.txt -f normalised.fa`;

print "Assigning taxonomy...\n";
`assign_taxonomy.py -i normalised.fa_rep_set.fasta -t /srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/taxonomies/otu_id_to_greengenes.txt -r /srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010.fasta -m blast -e 1e-50`;

#print "Treeing...\n";
#`align_seqs.py -i normalised.fa_rep_set.fasta -t /srv/whitlam/bio/db/gg/qiime_default/core_set_aligned.fasta.imputed -p 0.6`;
#`filter_alignment.py -i pynast_aligned/normalised.fa_rep_set_aligned.fasta`;
#`make_phylogeny.py -i normalised.fa_rep_set_aligned_pfiltered.fasta -r midpoint`;

print "Making OTU table...\n";
`make_otu_table.py -i uclust_picked_otus/normalised_otus.txt -t blast_assigned_taxonomy/normalised.fa_rep_set_tax_assignments.txt -o $global_working_dir/results/otu_table.txt`;


print "Rarefaction and diversity...\n";
`per_library_stats.py -i $global_working_dir/results/otu_table.txt > otu_table_stats.txt`;

#need to grep for params
my $min_seqs  = 1; #will not compute values for zero
my $max_seqs  = get_max_lib_size( 'otu_table_stats.txt' );
my $num_reps = 10;
my $num_steps = 30;
my $step_size = int(($max_seqs-$min_seqs)/$num_steps) || 1;
`multiple_rarefactions.py -i $global_working_dir/results/otu_table.txt  -o alpha_rare/rarefaction -m $min_seqs -x $max_seqs -n $num_reps -s $step_size`;
`alpha_diversity.py -i alpha_rare/rarefaction/ -o alpha_rare/alpha_div/ -m observed_species,shannon`;
`collate_alpha.py -i alpha_rare/alpha_div/ -o alpha_rare/alpha_div_collated/`;

# Rarefaction plot: All metrics command
`make_rarefaction_plots.py -i alpha_rare/alpha_div_collated/ -m $global_working_dir/QA/qiime_mapping.txt -o $global_working_dir/results/alpha_rare/alpha_rarefaction_plots/ --background_color white --resolution 75 --imagetype svg`;
#`beta_diversity.py -i otu_table.txt -t normalised.fa_rep_set_aligned_pfiltered.tre -m weighted_unifrac,unweighted_unifrac -o $global_working_dir/results/beta_diversity`;


print "Normalizing OTU table...\n";
if(0 <= $global_norm){
    `multiple_rarefactions_even_depth.py -i $global_working_dir/results/otu_table.txt -o $global_working_dir/rare_tables/ -d $global_norm -n 100 --lineages_included --k`;
    `average_tables.pl $global_working_dir/rare_tables/ $global_working_dir/results/collated_otu_table.txt`;
    
}

print "Summarizing by taxa.....\n";
`summarize_taxa.py -i $global_working_dir/results/collated_otu_table.txt -o $global_working_dir/results`;

print "Generating Genus-level heat map.....\n";
`getColors.pl $global_working_dir/results/collated_otu_table_L6.txt $global_working_dir/results/color_file.txt`;
`R --vanilla --slave --args $global_working_dir/results/collated_otu_table_L6.txt $global_working_dir/results/HeatMap.pdf $global_working_dir/results/color_file.txt < $Bin/HeatMap.R > R.stdout`;


print "Results are located in: $global_working_dir/results/\n";

######################################################################
# CUSTOM SUBS
######################################################################

sub get_max_lib_size {
   my ($stat_file) = @_; 

   # Given a QIIME library stat file (generated by per_library_stats.py)
   # get the number of sequences in the largest library

   # Num samples: 2
   #
   # Seqs/sample summary:
   #  Min: 336
   #  Max: 351
   #  Median: 343.5
   #  Mean: 343.5
   #  Std. dev.: 7.5
   #  Median Absolute Deviation: 7.5
   # Default even sampling depth in
   #  core_qiime_analyses.py (just a suggestion): 336
   #
   # Seqs/sample detail:
   #  211: 336
   #  212: 351
   #
   # Total observations (sequences): 687

   my $max_seqs = 0;
   open my $in, '<', $stat_file or die "Error: Could not read file $stat_file\n$!\n";
   while (my $line = <$in>) {
      chomp $line;
      if ( $line =~ /^\s*Max:\s*([0-9]+)\s*$/i ) {
	$max_seqs = $1;
        last;
      }
   }
   close $in;
   return $max_seqs;
}


sub parse_config_results
{
    #-----
    # parse the app config file and produce a qiime mappings file
    #
    open my $conf_fh, "<", $options->{'config'} or die "Error: Could not read config file ".$options->{'config'}."\n$!\n";
    while(<$conf_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;
        
        # save the MID for later
        $global_samp_ID_list{$fields[$FNA{'SampleID'}]} = $fields[$FNA{'USE'}];
    }
    
    # user options section
    while(<$conf_fh>)
    {
        chomp $_;
        my @fields = split /=/, $_;
        if($#fields > 0)
        {
            if($fields[0] eq "NORMALISE")
            {
                # is this guy set?
                if($fields[1] ne "")
                {
                    $global_norm = int($fields[1]);
                }
            }
            elsif($fields[0] eq "MUL_RARE_M")
            {
                # is this guy set?
                if($fields[1] ne "")
                {
                    $global_rare_M = int($fields[1]);
                }
            }
            elsif($fields[0] eq "MUL_RARE_X")
            {
                # is this guy set?
                if($fields[1] ne "")
                { 
                    $global_rare_X = int($fields[1]);
                }
            }
            elsif($fields[0] eq "MUL_RARE_S")
            {
                # is this guy set?
                if($fields[1] ne "")
                { 
                    $global_rare_S = int($fields[1]);
                }
            }
            elsif($fields[0] eq "MUL_RARE_N")
            {
                # is this guy set?
                if($fields[1] ne "")
                { 
                    $global_rare_N = int($fields[1]);
                }
            }
        }
    }    
    close $conf_fh;
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "config|c:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    if(!exists $options{'config'} ) { print "**ERROR: you MUST give a config file\n"; exec("pod2usage $0"); }
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
   print '';
}

__DATA__

=head1 NAME

    app_make_results.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Insert detailed description here

=head1 SYNOPSIS

    app_make_results.pl -c|config CONFIG_FILE [-help|h]

      -c CONFIG_FILE               app config file to be processed
      [-help -h]                   Displays basic usage information
         
=cut

