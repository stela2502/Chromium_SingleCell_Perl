#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-02-21 Stefan Lang

  This program is free software; you can redistribute it 
  and/or modify it under the terms of the GNU General Public License 
  as published by the Free Software Foundation; 
  either version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License 
  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1  SYNOPSIS

    SplitBam4Quant.pl
       -infile       :the input bam file cerated from a fastq file processed by SplitToCells.pl
       -cell_ids     :the per_cell_read_count.xls table from SplitToCells.pl
       -outpath      :the outpath for the split bam files
       -min_reads    :which cutoff should I apply to the reads?


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Take a bam file and the per_cell_read_count.xls file and create usable bam files with at least min_reads number of reads in.

  To get further help use 'SplitBam4Quant.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $cell_ids, $outpath, $min_reads);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-cell_ids=s"    => \$cell_ids,
	 "-outpath=s"    => \$outpath,
	 "-min_reads=s"    => \$min_reads,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infile) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $cell_ids) {
	$error .= "the cmd line switch -cell_ids is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $min_reads) {
	$error .= "the cmd line switch -min_reads is undefined!\n";
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
	print "$errorMessage.\n";
	pod2usage(q(-verbose) => 1);
}



my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/SplitBam4Quant.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -cell_ids '$cell_ids'" if (defined $cell_ids);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= " -min_reads '$min_reads'" if (defined $min_reads);



mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_SplitBam4Quant.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

# you know what - I quantify the thing right here...

die "Not implemented - use the QuantifyBamFile.pl instead!\n";
