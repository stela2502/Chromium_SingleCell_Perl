#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-02-21 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c59a74e96bf8f751bb43a81e8e7cd1cef7c84747
   

=head1  SYNOPSIS

    DropUselessFastqReads.pl
       -path       :The output path from SplitTiCells.pl
       -sampleID   :The sample ID and SplitTiCells option oname
       -minReads   :Minimum read count for the usage of the cellID in one fastq file (default nextSeq 100)
       -outfile    :the new cleaned fastq file


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Use the output from SplitToCells and drop all reads from the fastq files 
  that do not belong to a cell with more than 100 reads in any of the analysies.

  To get further help use 'DropUselessFastqReads.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::FastqFile;
use IO::File;
use File::Spec::Functions;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $path, $sampleID, $minReads, $outfile);

Getopt::Long::GetOptions(
	 "-path=s"    => \$path,
	 "-sampleID=s"    => \$sampleID,
	 "-minReads=s"    => \$minReads,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $path) {
	$error .= "the cmd line switch -path is undefined!\n";
}
unless ( defined $sampleID) {
	$error .= "the cmd line switch -sampleID is undefined!\n";
}
unless ( defined $minReads) {
	$minReads = 100;
}
unless ( defined $outfile) {
	$outfile = File::Spec->catfile( $path,$sampleID );
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

$task_description .= 'perl '.$plugin_path .'/DropUselessFastqReads.pl';
$task_description .= " -path '$path'" if (defined $path);
$task_description .= " -sampleID '$sampleID'" if (defined $sampleID);
$task_description .= " -minReads '$minReads'" if (defined $minReads);
$task_description .= " -outfile '$outfile'" if (defined $outfile);


die "This script is not worth it's running time as it not even drops 10% of the reads.\n"
  . "Rather use the 'usable sample' information in the processing of the mapped data.\n";


use stefans_libs::Version;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version '.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

## first I need to get all xls files and read them:
## use only cell ids with more than $minReads reads

opendir(DIR, $path ) or die "I could not oopen the path $path\n$!";

my @samplefiles = grep { !/^\./ }  grep { $sampleID } readdir(DIR) ;
closedir( DIR );

my ($OK, @tmp, $reads);
foreach my $file ( grep { /xls$/ } @samplefiles ) {
	print "I process file $file\n";
	open ( IN , "<$file" ) or die "I could not open file '$file'\n$!\n";
	while( <IN> ) {
		next if $_ =~ m/^#/;
		chomp;
		@tmp = split("\t", $_ );
		next if ( $tmp[1] eq "count");
		if ( $tmp[1] >= $minReads ) {
			$reads = $tmp[1];
			@tmp = split("_", $tmp[0] );
			$OK->{$tmp[3]} += $reads;
		}
	}
	close ( IN );
}
my $sum = 0;
map{ $sum += $_ } values %$OK;

print "I got ".scalar( keys %$OK ) . " passing samples with a total of $sum reads\n";

## so now lets filter the fastq file - shall we?
my ($cellID, $passing, $filtered_out, $OUT, $files );

$passing = $filtered_out = $files =0;

$outfile .= ".passing.fastq.gz" unless( $outfile=~m/passing.fastq.gz$/ );

open( $OUT, "| /bin/gzip -c > $outfile" )
		  or die
		  "I could not open the out pipe '| /bin/gzip -c > $outfile'\n$!\n";
print "Established output pipe '/bin/gzip -c > $outfile'\n";
my $func = sub {
	my ( $fastqfile, $entry ) = @_;

	#@NB501227:131:HJ2GYBGX5:1:11101:9587:1077:S_ATTCCGAT_C_TGGTTCCGTCGGATCC_CTGGATTGGG 2:N:0:ATTCCGAT
	@tmp = split( "_", $entry->name() );
	$cellID = $tmp[3];
    	
	if (  $OK->{$cellID} ) {
		#warn "$cellID did match and had $OK->{$cellID} reads.\n";
		$entry->filter_low_quality(20);
		
		if ( length($entry->sequence()) > 80 ){
			$entry->write($OUT);
			$passing ++;
		}else {
			#warn "$cellID faled\n"; 
			$filtered_out++;
		}
		#die if ( $passing > 50)
	}
	$entry->clear();
	
};


my $worker = stefans_libs::FastqFile->new();
my $i = 0;
foreach my $file (  grep { /fastq.gz$/ } @samplefiles) {
	print"Processing fastq file '$file'\n";
	$files ++;
	$worker->filter_file( $func, $file );
	print "Done with $file I have kept $passing reads and dropped $filtered_out ones (summary of all past files).\n";
	
}


close ( $OUT );

print "Done with $files fastq files I have kept $passing reads and dropped $filtered_out ones.\n";
print "Results are in '$outfile'\n";
