#! /usr/bin/env perl


=head1 NAME

opp_analyze.pl

=head1 DESCRIPTION

Analyse 16S rRNA sequence data: produce OTU tables, normalize data by
number of sequences, generate heat maps and rarefaction curves. Most
steps use QIIME scripts.

=head1 REQUIRED ARGUMENTS

=over

=item -c <config_file> | --config_file <config_file>

OPP config file to use

=for Euclid:
   config_file.type: readable

=back

=head1 OPTIONS

=over

=item -s <sample_size>... | --sample_size <sample_size>...

If normalization by library size it to be performed, provide a list of preferred
sample sizes (i.e. number of sequences) to use for normalization, in decreasing
order of preference. Default: sample_size.default

=for Euclid:
   sample_size.type: int, sample_size > 0
   sample_size.default: [ 10000, 5000, 1000, 500, 100 ]

=item -n <num_reps> | --num_reps <num_reps>

Number of repetitions to do for the normalization by library size. Use 0 to not
do any normalization. Default: num_reps.default

=for Euclid:
   num_reps.type: int, num_reps >= 0
   num_reps.default: 1000

=back

=head1 DEPENDENCIES

You run this program, you need Perl and the following CPAN Perl modules:

=over

=item *

Getopt::Euclid

=back

In addition, you need these bioinformatic programs installed:

=over

=item *

QIIME

=back
         
=head1 COPYRIGHT

Copyright 2011 Florent Angly and Dana Willner

Originally based on the APP scripts, copyright Michael Imelfort

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

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw( min );

# From CPAN
use Getopt::Euclid qw( :minimal_keys );

# Local OPP helper module from Perl script folder location
use FindBin qw($Bin);
use lib "$Bin";
use OPPConfig;

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# GLOBALS
my $global_working_dir = "~";
my $global_conf_full = $ARGV{'config_file'};
my $global_conf_base = basename $global_conf_full;
my %global_samp_ID_list = ();

print "Checking config...\n";

# Get the conf base
my @cb_1 = split /_/, $global_conf_base;
my @cb_2 = split /\./, $cb_1[1];
$global_conf_base = $cb_2[0];

# GET THE WORKING DIR
# working dir is the dir of the config file
# get the present dir
my $pwd = `pwd`;
chomp $pwd;
$global_working_dir = dirname($pwd."/".$global_conf_full);

# Override values from config
getSamples( $global_conf_full );

# get the working directories
getWorkingDirs($global_conf_full);
my $full_proc_dir = "$global_working_dir/$proc_dir";


# make the output directories
makeOutputDirs("");

# NORMALISE
`cp $global_working_dir/QA/denoised_acacia/acacia_out__all_tags.seqOut $full_proc_dir/normalised.fa`; 

print "Clustering OTUs...\n";
`pick_otus.py -i $full_proc_dir/normalised.fa -o $full_proc_dir/uclust_picked_otus`;

print "Picking OTU representative sequences...\n";
`pick_rep_set.py -i $full_proc_dir/uclust_picked_otus/normalised_otus.txt -f $full_proc_dir/normalised.fa`;

print "Assigning taxonomy...\n";
`assign_taxonomy.py -i $full_proc_dir/normalised.fa_rep_set.fasta -o $full_proc_dir/blast_assigned_taxonomy/ -t /srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/taxonomies/otu_id_to_greengenes.txt -r /srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010.fasta -m blast -e 1e-50`;

#print "Treeing...\n";
#`align_seqs.py -i $full_proc_dir/normalised.fa_rep_set.fasta -t /srv/whitlam/bio/db/gg/qiime_default/core_set_aligned.fasta.imputed -p 0.6`;
#`filter_alignment.py -i $full_proc_dir/pynast_aligned/normalised.fa_rep_set_aligned.fasta`;
#`make_phylogeny.py -i $full_proc_dir/normalised.fa_rep_set_aligned_pfiltered.fasta -r midpoint`;

print "Making OTU table...\n";
`make_otu_table.py -i $full_proc_dir/uclust_picked_otus/normalised_otus.txt -t $full_proc_dir/blast_assigned_taxonomy/normalised.fa_rep_set_tax_assignments.txt -o $global_working_dir/results/otu_table.txt`;


print "Rarefaction and diversity...\n";
`per_library_stats.py -i $global_working_dir/results/otu_table.txt > $full_proc_dir/otu_table_stats.txt`;

#need to grep for params

my ($lib_min_seqs, $lib_max_seqs) = get_lib_stats( "$full_proc_dir/otu_table_stats.txt" );

my $rare_min_seqs  = 1; # unfortunately, will not compute values for zero
my $rare_max_seqs  = $lib_max_seqs;
my $rare_num_reps  = 20;
my $rare_num_steps = 30;
my $rare_step_size = int(($rare_max_seqs - $rare_min_seqs)/$rare_num_steps) || 1;
`multiple_rarefactions.py -i $global_working_dir/results/otu_table.txt  -o $full_proc_dir/alpha_rare/rarefaction -m $rare_min_seqs -x $rare_max_seqs -n $rare_num_reps -s $rare_step_size`;
`alpha_diversity.py -i $full_proc_dir/alpha_rare/rarefaction/ -o $full_proc_dir/alpha_rare/alpha_div/ -m observed_species,shannon`;
`collate_alpha.py -i $full_proc_dir/alpha_rare/alpha_div/ -o $full_proc_dir/alpha_rare/alpha_div_collated/`;

# Rarefaction plot: All metrics command
`make_rarefaction_plots.py -i $full_proc_dir/alpha_rare/alpha_div_collated/ -m $global_working_dir/QA/qiime_mapping.txt -o $global_working_dir/results/alpha_rare/alpha_rarefaction_plots/ --background_color white --resolution 75 --imagetype svg`;
#`beta_diversity.py -i otu_table.txt -t normalised.fa_rep_set_aligned_pfiltered.tre -m weighted_unifrac,unweighted_unifrac -o $global_working_dir/results/beta_diversity`;


print "Normalizing OTU table...\n";
my $norm_num_reps = $ARGV{num_reps};
if( $norm_num_reps > 0) {
    my $norm_sample_size = pick_best_sample_size($lib_min_seqs, $ARGV{'sample_size'});
    `multiple_rarefactions_even_depth.py -i $global_working_dir/results/otu_table.txt -o $full_proc_dir/rare_tables/ -d $norm_sample_size -n $norm_num_reps --lineages_included --k`;
    `average_tables.pl $full_proc_dir/rare_tables/ $global_working_dir/results/normalized_otu_table.txt`;
    
} else {
   die "PROBABLY NEED TO DO SOMETHING WHEN NO NORMALIZATION REQUIRED\n";
}

print "Summarizing by taxa.....\n";
`summarize_taxa.py -i $global_working_dir/results/normalized_otu_table.txt -o $global_working_dir/results`;

print "Generating Genus-level heat map.....\n";
`getColors.pl $global_working_dir/results/normalized_otu_table_L6.txt $global_working_dir/results/color_file.txt`;
`R --vanilla --slave --args $global_working_dir/results/normalized_otu_table_L6.txt $global_working_dir/results/HeatMap.pdf $global_working_dir/results/color_file.txt < $Bin/HeatMap.R > $full_proc_dir/R.stdout`;


print "Results are located in: $global_working_dir/results/\n";


#------------------------------------------------------------------------------#


sub pick_best_sample_size {
   # Given the number of sequences in the smallest libraries, pick the first
   # number of sequences smaller than that in a list of possible values
   my ($smallest_lib, $sample_sizes) = @_;
   my $best;
   for my $sample_size (@$sample_sizes) {
      if ($sample_size < $smallest_lib) {
         $best = $sample_size;
         last;
      }
   }
   if (not defined $best) {
      my $min_sample_size = min @$sample_sizes;
      die "Error: Could not find a suitable sample size. The smallest sample ".
         "requested for normalization was $min_sample_size sequences but your ".
         "smallest library had only $smallest_lib sequences.\n"; 
   }
   return $best;
}


sub get_lib_stats {
   my ($stat_file) = @_; 

   # Given a QIIME library stat file (generated by per_library_stats.py)
   # get the number of sequences in the smallest and largest library

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

   my ($min_seqs, $max_seqs);
   open my $in, '<', $stat_file or die "Error: Could not read file $stat_file\n$!\n";
   while (my $line = <$in>) {
      chomp $line;
      if ( $line =~ /^\s*Max:\s*([0-9]+)\s*$/i ) {
	$max_seqs = $1;
        last;
      } elsif ( $line =~ /^\s*Min:\s*([0-9]+)\s*$/i ) {
        $min_seqs = $1;
      }
   }
   close $in;
   return $min_seqs, $max_seqs;
}


sub getSamples
{
    #-----
    # parse the app config file and return which samples to use
    #
    my ($conf_file) = @_;
    my $global_norm;
    open my $conf_fh, "<", $conf_file or die "Error: Could not read config file $conf_file\n$!\n";
    while(<$conf_fh>) {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;
        
        # save the MID for later
        $global_samp_ID_list{$fields[$FNA{'SampleID'}]} = $fields[$FNA{'USE'}];
    }
    
    close $conf_fh;
    return $global_norm;
}

