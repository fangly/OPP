###############################################################################
#
#    OppConfig.pm
#    
#    Makes more useful names for the fields in the config file
#    The app_* scripts should include this file first
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
package OPPConfig;
require Exporter;
use File::Basename;

our @ISA = qw(Exporter);
our @EXPORT=qw(
   %FNB
   %FNA
   $FNB_HEADER
   $FNA_HEADER
   $FNA_LINE_FINISHER
   $FNA_FOOTER
   $OPP_ROOT
   $OPP_RAW
   $OPP_BYJOB
   $OPP_BYRUN
   %global_samp_ID_list
   %global_raw_counts
   %global_chimer_counts
   %global_acacia_counts
   $global_acacia_config
   $global_barcode_length
   $QA_dir
   $proc_dir
   $res_dir
   $global_acacia_output_dir
   $global_working_dir
   $global_mapping_file 
   $QIIME_split_out
   getWorkingDirs
   makeOutputDirs
   splitLibraries 
   removeChimeras
   denoise
   getReadCounts
   parseConfigQA
   updateConfigQA
);

#
# A file is created in PyroDB which can be used to split the sff file and 
# make all the relavant job dirs etc... XXX.pdbm
# This file is basically a qiime mapping file, BUT it has the same name as
# the sff. (or fasta and qual). The format is given below:
#
# SampleID	BarcodeSequence	LinkerPrimerSequence	Description
# <JID.SID> MID             acgggcggtgtgtRc         <PDB sample name>
# ...
#
our $FNB_HEADER = "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription";
our %FNB = ();
$FNB{'SampleID'} = 0;
$FNB{'BarcodeSequence'} = 1;
$FNB{'LinkerPrimerSequence'} = 2;
$FNB{'Description'} = 3;

#
# Once the sff has been munged, each job will be placed into a folder in the by_jobID dir
# and given an app config file. The format is given below:
#
# #SampleID	BarcodeSequence	LinkerPrimerSequence	Description	        RAW	CHIME	ACC	USE
# <SID>     MID             acgggcggtgtgtRc         <PDB sample name>   XX  XX      XX  XX
# ...
#
our $FNA_HEADER = "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tRAW\tCHIME\tACC\tUSE";
our $FNA_LINE_FINISHER = "\tXX\tXX\tXX\t1\n";
our $FNA_FOOTER = "@@\
NORMALISE=\
DB=0\
MUL_RARE_M=\
MUL_RARE_X=\
MUL_RARE_S=\
MUL_RARE_N=";

our %FNA = ();
$FNA{'SampleID'} = 0;
$FNA{'BarcodeSequence'} = 1;
$FNA{'LinkerPrimerSequence'} = 2;
$FNA{'Description'} = 3;
$FNA{'RAW'} = 4;
$FNA{'CHIME'} = 5;
$FNA{'ACC'} = 6;
$FNA{'USE'} = 7;

#
# The OPP_ROOT environment variable should be set. there are a number of dirs we need to get from there
#
our $OPP_ROOT = `echo \$OPP_ROOT`;
chomp $OPP_ROOT;
our $OPP_RAW = $OPP_ROOT."/raw";
our $OPP_BYJOB = $OPP_ROOT."/by_jobid";
our $OPP_BYRUN = $OPP_ROOT."/by_run";

#
# We make a number of directories durig the process. Store their names here
#
our $QA_dir = "QA";
our $proc_dir = "processing";
our $res_dir = "results";

our $global_acacia_output_dir = "UNSET";
our $global_working_dir = "UNSET";
our $global_mapping_file = "UNSET";

#
# We use a number of config file which are located all over the place. Just in case we move them., we can store them here...
#
our $global_acacia_config = "/srv/whitlam/bio/apps/sw/app/beta/app_acacia.config";
our $global_barcode_length = "variable_length";

#
# Some programs make special output files. Store these filenames here
#
my $QIIME_map_file = "qiime_mapping.txt";
our $QIIME_split_out = "seqs.fna";
my $CHIME_gg_file = "/srv/whitlam/bio/db/gg/qiime_default/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010.fasta";
my $CHIME_good_file = "good.fasta";
my $CHIME_bad_file = "ch.fasta";
my $ACACIA_out_file = "acacia_out__all_tags.seqOut";

#
# Some global variables we use to store read counts and sample IDs
# 
our %global_samp_ID_list = ();          # list of sample IDs
our %global_raw_counts = ();            # raw reads in $QA_dir/$QIIME_split_out
our %global_chimer_counts = ();
our %global_acacia_counts = ();

######################################################################
# SHARED SUBS
######################################################################

sub getWorkingDirs
{
    #-----
    # Set a number of output directories
    #
    my ($config_prefix) = @_;
    
    # get the acacia denoised directory
    my $acc_dir_raw = `grep OUTPUT_DIR $global_acacia_config`;
    chomp $acc_dir_raw;
    my @acc_dir_fields = split /=/, $acc_dir_raw;
    $global_acacia_output_dir = $acc_dir_fields[1];

    # Get the working dir
    # working dir is the dir of the config file
    # get the present dir
    my $pwd = `pwd`;
    chomp $pwd;
    $global_working_dir = dirname("$pwd/$config_prefix");
    
    # set the mapping file
    $global_mapping_file = "$global_working_dir/$QA_dir/$QIIME_map_file";
}

sub makeOutputDirs
{
    #-----
    # Directories must be made before we can put files there
    #
    my ($job_dir) = @_;
    `mkdir -p $global_working_dir$job_dir/$QA_dir`;
    `mkdir -p $global_working_dir$job_dir/$QA_dir/$global_acacia_output_dir`;
    `mkdir -p $global_working_dir$job_dir/$proc_dir`;
    `mkdir -p $global_working_dir$job_dir/$res_dir`;
}

sub splitLibraries
{
    #-----
    # Wrapper for Qiime split libraries
    #
    my ($job_ID) = @_;
    print "Splitting libraries...\n";
    `split_libraries.py -m $QIIME_map_file -f ../$job_ID.fna -b $global_barcode_length`;
}

sub removeChimeras
{
    #-----
    # Remove chimeras using uclust
    #
    print "Removing chimeras...\n";
    `usearch --uchime seqs.fna --db $CHIME_gg_file --nonchimeras $CHIME_good_file --chimeras $CHIME_bad_file`;
}

sub denoise
{
    #-----
    # run acacia on the data
    #
    print "Denoising using acacia...\n";
    `java -jar \$ACACIA -c $global_acacia_config`;
    `sed -i -e "s/all_tags_[^ ]* //" $global_acacia_output_dir/$ACACIA_out_file`;
}

sub getReadCounts
{
    #-----
    # get the read counts for raw, chimera removed and acacia filtered sequences
    #
    # this guy is called from within the QA dir so the files are local
    #
    # the three files to parse are:
    # $QIIME_split_out
    # $CHIME_good_file
    # $global_acacia_output_dir/$ACACIA_out_file
    #
    open my $tmp_fh, "<", $QIIME_split_out or die $!;
    while(<$tmp_fh>)
    {
        next if(!($_ =~ /^>/));
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $fl =$f2[1];  
        foreach my $uid (keys %global_samp_ID_list)
        {
            if($fl =~ /^$uid/)
            {
                # this guy begins with the exact MID
                $global_raw_counts{$uid}++;
                next;
            }
        }
    }
    
    open $tmp_fh, "<", $CHIME_good_file or die $!;
    while(<$tmp_fh>)
    {
        next if(!($_ =~ /^>/));
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $fl =$f2[1];  
        foreach my $uid (keys %global_samp_ID_list)
        {
            if($fl =~ /^$uid/)
            {
                # this guy begins with the exact MID
                $global_chimer_counts{$uid}++;
                next;
            }
        }
    }    

    open $tmp_fh, "<", "$global_acacia_output_dir/$ACACIA_out_file" or die $!;
    while(<$tmp_fh>)
    {
        next if(!($_ =~ /^>/));
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $fl =$f2[1];  
        foreach my $uid (keys %global_samp_ID_list)
        {
            if($fl =~ /^$uid/)
            {
                # this guy begins with the exact MID
                $global_acacia_counts{$uid}++;
                next;
            }
        }
    }
}

sub parseConfigQA
{
    #-----
    # parse the app config file and produce a qiime mappings file
    #
    my ($config_prefix) = @_;
    open my $conf_fh, "<", $config_prefix or die $!;
    open my $mapping, ">", $global_mapping_file or die $!;
    print $mapping "$FNB_HEADER\n";
    while(<$conf_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;

        # save the MID for later
        $global_samp_ID_list{$fields[$FNA{'SampleID'}]} = $fields[$FNA{'USE'}];
        $global_raw_counts{$fields[$FNA{'SampleID'}]} = 0;
        $global_chimer_counts{$fields[$FNA{'SampleID'}]} = 0;
        $global_acacia_counts{$fields[$FNA{'SampleID'}]} = 0;

        if("1" eq $fields[$FNA{'USE'}])
        {
            print $mapping "$fields[$FNA{'SampleID'}]\t$fields[$FNA{'BarcodeSequence'}]\t$fields[$FNA{'LinkerPrimerSequence'}]\t$fields[$FNA{'Description'}]\n";
        }
    }
    close $conf_fh;
    close $mapping;
}

sub updateConfigQA
{
    #-----
    # parse the app config  and update to include read counts
    # this guy is called from within the QA dir so we need to do a ../ on the file names
    #
    my ($config_prefix) = @_;
    open my $conf_fh, "<", "../$config_prefix" or die $!;
    open my $conf_fh_tmp, ">", "../$config_prefix.tmp" or die $!;
    while(<$conf_fh>)
    {
        if($_ =~ /^#/) { print $conf_fh_tmp $_; next; }
        if($_ =~ /^@/) { print $conf_fh_tmp $_; last; }
        chomp $_;
        my @fields = split /\t/, $_;
        print $conf_fh_tmp "$fields[$FNA{'SampleID'}]\t$fields[$FNA{'BarcodeSequence'}]\t$fields[$FNA{'LinkerPrimerSequence'}]\t$fields[$FNA{'Description'}]\t$global_raw_counts{$fields[$FNA{'SampleID'}]}\t$global_chimer_counts{$fields[$FNA{'SampleID'}]}\t$global_acacia_counts{$fields[$FNA{'SampleID'}]}\t$fields[$FNA{'USE'}]\n";
    }
    
    # just print out the rest of the file
    while(<$conf_fh>)
    {
        print $conf_fh_tmp $_;
    }

    close $conf_fh;
    close $conf_fh_tmp;
    
    my $mv_string  = "mv ../$config_prefix.tmp ../$config_prefix";
    `$mv_string`;
}


1;
