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
#### GLOBALS
my $global_working_dir = "~";
my $global_conf_base = basename($options->{'config'});
my %global_samp_ID_list = ();
my $global_norm = 0;
my $global_rare_M = 10;
my $global_rare_X = 100;
my $global_rare_S = 10;
my $global_rare_N = 10;

print "Checking if all the config checks out...\t\t";

#### Get the conf base
my @cb_1 = split /_/, $global_conf_base;
my @cb_2 = split /\./, $cb_1[1];
$global_conf_base = $cb_2[0];

#### GET THE WORKING DIR
# working dir is the dir of the config file
# get the present dir
my $pwd = `pwd`;
chomp $pwd;
$global_working_dir = dirname($pwd."/".$options->{'config'});

#### Override values from config
parse_config_results();

#### Start the results pipeline!
print "All good!\n";

#### NORMALISE
# TODO do this!
#print "Normalising...\n";
#if(0 >= $global_norm)
`cp $global_working_dir/QA/denoised_acacia/acacia_out__all_tags.seqOut $global_working_dir/processing/normalised.fa`; 
#else
#{ `app_normalise_reads.pl -in $global_working_dir/QA/denoised_acacia/acacia_out__all_tags.seqOut -out $global_working_dir/processing/normalised.fa -num $global_norm`; }

chdir "$global_working_dir/processing";

print "Clustering OTUs...\n";
`pick_otus.py -i normalised.fa`;

print "Picking OTU representative sequences...\n";
`pick_rep_set.py -i uclust_picked_otus/normalised_otus.txt -f normalised.fa`;

print "Assigning taxonomy...\n";
`assign_taxonomy.py -i normalised.fa_rep_set.fasta -t /srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/taxonomies/otu_id_to_greengenes.txt -r /srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010.fasta -m blast`;

#print "Treeing...\n";
#`align_seqs.py -i normalised.fa_rep_set.fasta -t /srv/whitlam/bio/db/gg/qiime_default/core_set_aligned.fasta.imputed -p 0.6`;
#`filter_alignment.py -i pynast_aligned/normalised.fa_rep_set_aligned.fasta`;
#`make_phylogeny.py -i normalised.fa_rep_set_aligned_pfiltered.fasta -r midpoint`;

print "Making OTU table...\n";
`make_otu_table.py -i uclust_picked_otus/normalised_otus.txt -t blast_assigned_taxonomy/normalised.fa_rep_set_tax_assignments.txt -o otu_table.txt`;

print "Normalizing otu table...\n";
if(0 <= $global_norm){
    `multiple_rarefactions_even_depth.py -i otu_table.txt -o $global_working_dir/rare_tables/ -d $global_norm -n 100 --lineages_included --k`;
    `average_tables.pl $global_working_dir/rare_tables/ $global_working_dir/results/collated_otu_table.txt`;
    
}

print "Summarizing by taxa.....\n";
`summarize_taxa.py -i $global_working_dir/collated_otu_table.txt -o $global_working_dir/results`;

print "Generating Genus-level heat map.....\n";
`getColors.pl $global_working_dir/results/collated_otu_table_L6.txt $global_working_dir/results/color_file.txt`;
`R --vanilla --slave --args $global_working_dir/results/collated_otu_table_L6.txt $global_working_dir/results/HeatMap.pdf $global_working_dir/results/color_file.txt < $global_working_dir/HeatMap.R`;


print "Rarefaction and diversity...\n";
#`multiple_rarefactions.py -i otu_table.txt -o rarefied_otu_tables/ -m $global_rare_M -x $global_rare_X -s $global_rare_S -n $global_rare_N`;
#`alpha_diversity.py -i rarefied_otu_tables/ -o rarefied_alpha/ -m chao1,chao1_confidence,observed_species,simpson,shannon`;
#`collate_alpha.py -i rarefied_alpha/ -o $global_working_dir/results/alpha_diversity/`;
# find step size for rarefaction
`per_library_stats.py -i otu_table.txt -o $global_working_dir/results/otu_table_stats.txt`;
#need to grep for params
my $min_seqs=0;
my $max_seqs=400;
my $num_steps=10;
my $step_size=($max_seqs-$min_seqs)/$num_steps;
`multiple_rarefactions.py -i otu_table.txt  -o $global_working_dir/results/alpha_rare/rarefaction -m $min_seqs -x $max_seqs -n $num_steps -s $step_size`;
`alpha_diversity.py -i $global_working_dir/results/alpha_rare/rarefaction/ -o $global_working_dir/results/alpha_rare/alpha_div/ -m observed_species,shannon`;
`collate_alpha.py -i $global_working_dir/results/alpha_rare/alpha_div/ -o $global_working_dir/results/alpha_rare/alpha_div_collated/`;

# Rarefaction plot: All metrics command
`make_rarefaction_plots.py -i $global_working_dir/results/alpha_rare/alpha_div_collated/ -m $global_working_dir/QA/qiime_mapping.txt -o $global_working_dir/results/alpha_rare/alpha_rarefaction_plots/ --background_color white --resolution 75 --imagetype svg`;
#`beta_diversity.py -i otu_table.txt -t normalised.fa_rep_set_aligned_pfiltered.tre -m weighted_unifrac,unweighted_unifrac -o $global_working_dir/results/beta_diversity`;

print "Tidy up...\n";
`cp otu_table.txt $global_working_dir/results/`;
#`cp normalised.fa_rep_set_aligned_pfiltered.tre $global_working_dir/results/`;

print "Results are located in: $global_working_dir/results/\n";

######################################################################
# CUSTOM SUBS
######################################################################
sub parse_config_results
{
    #-----
    # parse the app config file and produce a qiime mappings file
    #
    open my $conf_fh, "<", $options->{'config'} or die "Error: Could not read config file ".$options->{'config'}."\n$!\n;
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
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
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

