#! /usr/bin/env perl


=head1 NAME

opp.pl

=head1 DESCRIPTION

Pre-process and analyze 454 16S rRNA sequence data:

=over

=item 1.

Run opp_qa.pl

=item 2.

Run opp_analyze.pl

=back

=head1 REQUIRED ARGUMENTS

=over

=item <config_file>

OPP config file to use

=for Euclid:
   config_file.type: readable

=back

=head1 DEPENDENCIES

You run this program, you need Perl and the following CPAN Perl modules:

=over

=item *

Getopt::Euclid

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

# Run the pipeline
my $config_file = $ARGV{'config_file'};
run( "opp_qa.pl -c $config_file" );
run( "opp_analyze.pl -c $config_file" );

exit;

