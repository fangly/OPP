#! /usr/bin/env perl


=head1 NAME

opp_qa.pl

=head1 DESCRIPTION

Pre-process 454 16S rRNA sequence data:

=over

=item 1.

Split libraries by MID (QIIME)

=item 2.

Remove chimeras (UCHIME)

=item 3.

Correct sequencing errors (ACACIA)

=back

=head1 REQUIRED ARGUMENTS

=over

=item -c <config_file> | --config <config_file>

OPP config file to use

=for Euclid:
   config_file.type: readable

=back

=head1 OPTIONS

=over

=item -a <acacia_conf> | --acacia <acacia_conf>

Alternate Acacia config file to use

NOTE: If you specify a different acacia config file, then you must use the
following values, or this script will break!
      
    FASTA_LOCATION=good.fasta
    OUTPUT_DIR=denoised_acacia
    OUTPUT_PREFIX=acacia_out_
    SPLIT_ON_MID=FALSE

=for Euclid:
   acacia_conf.type: readable

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

=item *

UCHIME

=item *

ACACIA

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


print "Checking if all the config checks out...\t\t";
# acacia config file
if(exists $ARGV{'acacia'})
{
    # user supplied config file
    $global_acacia_config = $ARGV{'acacia'};
}
if (!(-e $global_acacia_config)) { die "Acacia config file: $global_acacia_config does not exist!\n"; }

# get the Job_ID we're working on
my $job_ID = basename($ARGV{'config'});
my @cb_1 = split /_/, $job_ID;
my @cb_2 = split /\./, $cb_1[1];
$job_ID = $cb_2[0];

# get the working directories
getWorkingDirs($ARGV{'config'});

# make the output directories
makeOutputDirs("");

# parse the config file
parseConfigQA($ARGV{'config'});

print "All good!\n";

#### start the $QA_dir pipeline!
chdir "$global_working_dir/$QA_dir" or die "Error: Could not cd to directory $global_working_dir/$QA_dir\n$!\n";

print "Global working dir is $global_working_dir\n";

splitLibraries($job_ID);
removeChimeras();
denoise();

#### Fix the config file
print "Fixing read counts...\n";
getReadCounts();
updateConfigQA($ARGV{'config'});
