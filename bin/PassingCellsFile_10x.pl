#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-06-08 Stefan Lang

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
   
   binCreate.pl from git@gitlab:stefanlang/Stefans_Lib_Esentials.git commit e7545a27e262074f058adf34953b79b184a19248
   

=head1  SYNOPSIS

    PassingCellsFile_10x.pl
       -infile       :<please add some info!>
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  parse a bam 10x file and create the passing cell file required from the quantification using the perl script.

  To get further help use 'PassingCellsFile_10x.pl -help' at the comman line.

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

$task_description .= 'perl '.$plugin_path .'/PassingCellsFile_10x.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outfile '$outfile'" if (defined $outfile);



use stefans_libs::Version;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version '.$V->version( 'Chrmium_SingleCell_Perl' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
use stefans_libs::BAMfile;

my $bam_file = stefans_libs::BAMfile->new();
my ( $sample_name, $UMI, $counter, $runs );

$runs = 0;

my $function = sub  {
	shift; ## get rid of the bam object
	my @bam_line = split( "\t", shift );
	return if ( $bam_line[0] =~ m/^\@/ );
	( $sample_name, $UMI ) = &sample_and_umi_my(@bam_line);
	#print "I got sample '$sample_name'" if ($debug);
	$counter->{$sample_name} ||= 0;
	$counter->{$sample_name}++;
	$runs ++;
};

$bam_file->filter_file( $infile, $function );


open( OUT, ">$outfile" )
  or die
"I could not open the file '$outfile'\n$!\n";
my $statement =  "#total of $runs reads processed\n".
"#no sample/UMI information in 0\n";

print OUT $statement;
print OUT "Cell_ID\tcount\n";
foreach my $key ( sort keys %$counter ) {
	print OUT "$key\t$counter->{$key}\n";
}
close(OUT);

print $statement;


sub sample_and_umi_my {
	my (@bam_line) = @_;

#NB501227:57:HGJM5BGX2:1:23206:13379:6742:S_GATCTCAG_C_TATGCCCTCCATTCTA_CAGGTGAAGC
	my ( $sample_name, $UMI );
	my @matching_IDs = split( ":", $bam_line[0] );
	@matching_IDs = split( "_", pop(@matching_IDs) );
	#print "Matching IDS:" . join(" -- ", @matching_IDs)."\n";
	if ( ( !defined $matching_IDs[2] ) or ( !defined $matching_IDs[3] ) ) {
		return sample_and_umi_cellranger(@bam_line);
		Carp::confess(
			"Sorry I lack the UMI information - have you used SplitToCells.pl?"
		);
	}
	$sample_name = $matching_IDs[3];

#$sample_name = $self->{'cell_id_to_name'}->{$sample_name} if ( defined $self->{'cell_id_to_name'}->{$sample_name});
	$UMI = $matching_IDs[4];
	return ( $sample_name, $UMI );
}

sub sample_and_umi_cellranger {
	my (@bam_line) = @_;
	my ( $sample_name, $UMI );
	foreach (@bam_line) {

		#$sample_name = $1 if ( $_ =~m/CB:Z:([ACTGN]+)-?\d*$/);
		$sample_name = $1 if ( $_ =~ m/CR:Z:([AGCTN]+)$/ );
		if ( $_ =~ m/UR:Z:([ACTGN]+)$/ ){
			$UMI  = $1 ;
		}elsif ( $_ =~ m/UB:Z:([ACTGN]+)$/ ) {
			$UMI  = $1 ;
		}
		
	}
	return ( $sample_name, $UMI );
}

