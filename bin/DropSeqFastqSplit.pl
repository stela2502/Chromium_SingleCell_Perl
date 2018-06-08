#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-06-07 Stefan Lang

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

    DropSeqFqastSplit.pl
       -fastq       :the Drop-Seq fastq file (PMID26000488)
       -outpath     :The outpath for the thousands of fastq files created


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  This script splits one drop-seq fastq file into different cells based on a 12 + 8 cell + UMI read start. Poly A sequences willl be removed and only seuquences with more that 20 bp usable sequence will be reported.

  To get further help use 'DropSeqFqastSplit.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::FastqFile;
use stefans_libs::root;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $fastq, $outpath);

Getopt::Long::GetOptions(
	 "-fastq=s"    => \$fastq,
	 "-outpath=s"    => \$outpath,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $fastq) {
	$error .= "the cmd line switch -fastq is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/DropSeqFqastSplit.pl';
$task_description .= " -fastq '$fastq'" if (defined $fastq);
$task_description .= " -outpath '$outpath'" if (defined $outpath);



mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_DropSeqFqastSplit.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );

## Do whatever you want!

## This will move the cell id and UMI into the read name so that the fastq file can be alinged in one run
## the postprocessing should be quite straight forward and can be performed f9ollowing approximately the 10x way.

my $worker = stefans_libs::FastqFile->new();

my $fm = root->filemap( $fastq );
my $ofile = "$fm->{'filename_core'}.fastq.gz";
$ofile =~ s/\.f?a?s?t?q\.fastq/.fastq/;
$ofile = File::Spec->catfile( $outpath, $ofile);
open ( my $OUT , "| gzip -9 > $ofile") or die $!;

my $i = 0;
my $OK = 0;
my $tmp;
my $runs = 0;

my $function = sub { 
	
	my ( $worker, $entry ) = @_;
	$runs++;
	if ( $runs % 1e4 == 0 ) {
		print ".";
	}
	my $cell =  substr($entry->sequence(),0,12) ;
	my $umi  =  substr($entry->sequence(),12, 8);
	$entry -> trim( 'start', 20 );
	
	if ( $entry->sequence() =~ m/([Aa]{5}[Aa]+)$/ ) {
		$entry->trim( 'end', length( $entry->sequence() ) - length($1) );
	}
	if ( $entry->sequence() =~ m/([Tt]{5}[Tt]+)$/ ) {
		$entry->trim( 'end', length( $entry->sequence() ) - length($1) );
	}
	$tmp = "$cell.$umi";
	unless (  $tmp =~ m/[nN]/  ){
		#warn "both no N's\n";
		if (  length( $entry->sequence()) > 20 ) {
			$OK ++;
			$entry->name(join ("_", $entry->name(),"N","N", $cell, $umi ) ); ## This way the 10x pipeline can handle the bam results!
			$entry->write($OUT);
			#warn "used entry $i\n";
		}
	}else {
		$tmp =~ tr/Nn/--/;
		#warn "$tmp +  ".$entry->sequence()." $i ".length( $entry->sequence())."\n";
	}
	$entry->clear();
	
};

if ( -t STDOUT ) { ## for the progress bar
	$| = 1;
}

$worker->filter_file( $function, $fastq );


close ( $OUT );

print "\n\nAll $OK fastq entries have been renamed and stored in '$ofile' (".($runs- $OK)." have been filtered due to Ns in the cell or umi part or short read length)\n";

