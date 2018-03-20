#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-03-20 Stefan Lang

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

=head1 CREATED BY
   
   binCreate.pl from git@gitlab:stefanlang/Stefans_Lib_Esentials.git commit d6027a714e3cf82cabf695899e4c58b469f305b4
   

=head1  SYNOPSIS

    toFastq.pl
       -file       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  A tool to get 10xpipeline read names fastq file from the cellranger aligned bam.

  To get further help use 'toFastq.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;
use stefans_libs::BAMfile;
use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $file);

Getopt::Long::GetOptions(
	 "-file=s"    => \$file,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $file) {
	$error .= "the cmd line switch -file is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/toFastq.pl';
$task_description .= " -file '$file'" if (defined $file);



my ( $runs,$sample_name, $UMI );
sub sample_and_umi_cellranger {
  my (@bam_line) = @_;
  my ( $sample_name, $UMI );
  foreach (@bam_line) {
	#$sample_name = $1 if ( $_ =~m/CB:Z:([ACTGN]+)-?\d*$/);
	$sample_name = $1 if ( $_ =~ m/CR:Z:([AGCTN]+)$/ );
	$UMI         = $1 if ( $_ =~ m/UR:Z:([ACTGN]+)$/ );
	}
  return ( $sample_name, $UMI );
}

my $bam_file = stefans_libs::BAMfile->new();

sub filter {
	shift;
	my @bam_line = split( "\t", shift );
	return if ( $bam_line[0] =~ m/^\@/ );

	#die if ( $runs == 100);
	$runs++;
	( $sample_name, $UMI ) = sample_and_umi_cellranger(@bam_line);
	print "@".$bam_line[0].join("_",":S","GATCTCAG","C",$sample_name,$UMI)."\n"
		.$bam_line[9]."\n+\n".$bam_line[10]."\n";

}
$bam_file->filter_file( $file, \&filter );


## Do whatever you want!

