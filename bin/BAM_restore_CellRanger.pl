#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-08-30 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
   

=head1  SYNOPSIS

    BAM_restore_CellRanger.pl
       -infile       :<please add some info!>
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Convert my cell id and UMI store to CellRanger UMI and cell id store

  To get further help use 'BAM_restore_CellRanger.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $outfile);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infile) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/BAM_restore_CellRanger.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outfile '$outfile'" if (defined $outfile);



use stefans_libs::Version;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version Chromium_SingleCell_Perl '.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG '#library version stefanl_libs-BAMfile'. $V->version( 'stefanl_libs-BAMfile' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
my $worker = stefans_libs::BAMfile -> new();
my ( @bam_line, $sample_name, $UMI, @tmp );
my $filter = sub { 
	my ( $BamFile, $line ) = @_;
	unless ( $line =~m/^@/ ){
		@bam_line = split( "\t", $line );
		( $sample_name, $UMI ) = sample_and_umi_my(@bam_line);
		@tmp = split(":", $bam_line[0]);
		pop(@tmp);
		$bam_line[0] = join(":", @tmp);
		$line = join("\t", @bam_line, 'CB:Z:'.$sample_name, 'UB:Z:'.$UMI );
	}
	#print $line."\n";
	$line;
};

$worker -> apply_function (
	$infile ,  
	$outfile, 
	$filter 
);

sub sample_and_umi_my {
	my (@bam_line) = @_;

#NB501227:57:HGJM5BGX2:1:23206:13379:6742:S_GATCTCAG_C_TATGCCCTCCATTCTA_CAGGTGAAGC
	my ( $sample_name, $UMI );
	my @matching_IDs = split( ":", $bam_line[0] );
	@matching_IDs = split( "_", pop(@matching_IDs) );

	if ( ( !defined $matching_IDs[2] ) or ( !defined $matching_IDs[3] ) ) {
		Carp::confess(
			"Sorry I lack the UMI information - have you used SplitToCells.pl?"
		);
	}
	$sample_name = $matching_IDs[3];
	$UMI = $matching_IDs[4];
	return ( $sample_name, $UMI );
}

