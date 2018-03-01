#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-02-26 Stefan Lang

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

    ReshuffleBamFiles.pl
       -bams       :the list of bam files to re-shuffle into chr specific, sorted bam files
       -outpath    :the outpath for the chr soreted bam files
       -coverage   :the genome coverage file already used for the mapping scripts
       -sampleID   :the name for this 10x sample
       -tmp_path   :a fast tnp path (default = '$SNIC_TMP')

       -n   :number of processors to use

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Take a list of sorted bam files and the genome lengh file to re-shuffle all bam files into per chromosome sorted summary bam files - important for the 10 pipeline.

  To get further help use 'ReshuffleBamFiles.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::root;
use Parallel::ForkManager;
use File::Spec::Functions;
use DateTime;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help,     $debug,    $database,   @bams, $outpath,
	$coverage, $tmp_path, $sampleID, $n
);

Getopt::Long::GetOptions(
	"-bams=s{,}"    => \@bams,
	"-outpath=s"    => \$outpath,
	"-coverage=s"   => \$coverage,
	"-sampleID=s" => \$sampleID,
	"-n=s"          => \$n,
	"-tmp_path=s"   => \$tmp_path,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( -f $bams[0] ) {
	$error .= "the cmd line switch -bams is undefined!\n";
}
unless ( defined $outpath ) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( -f $coverage ) {
	$error .= "the cmd line switch -coverage is undefined!\n";
}
unless ( defined $sampleID ) {
	$error .= "the cmd line switch -sampleID is undefined!\n";
}
unless ($n) {
	$n = 3;
}
unless ($tmp_path) {
	$tmp_path = '$SNIC_TMP';
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	print "$errorMessage.\n";
	pod2usage( q(-verbose) => 1 );
}

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/ReshuffleBamFiles.pl';
$task_description .= ' -bams "' . join( '" "', @bams ) . '"'
  if ( defined $bams[0] );
$task_description .= " -outpath '$outpath'"       if ( defined $outpath );
$task_description .= " -coverage '$coverage'"     if ( defined $coverage );
$task_description .= " -sampleID '$sampleID'" if ( defined $sampleID );
$task_description .= " -n $n";

my $start = DateTime->now();

print "$task_description\nSTART " . $start->time() . "\n";

mkdir($outpath) unless ( -d $outpath );
open( LOG, ">$outpath/" . $$ . "_ReshuffleBamFiles.pl.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

unless ( -d $tmp_path ) {
	mkdir($tmp_path) or die $!;
}

## Do whatever you want!

open( IN, "<$coverage" )
  or die "I could not open the coverage file '$coverage'\n$!\n";
my @chrs;
while (<IN>) {
	push( @chrs, ( split( "\t", $_ ) )[0] );
}
close(IN);

my $pm = Parallel::ForkManager->new($n);

my @ofiles;

FILES:
foreach my $file (@bams) {
	my $pid = $pm->start and next FILES;
	my $fm = root->filemap($file);
	my ( $cmd, $ofile );
	unless ( -f "$file.bai" ) {
		system("samtools index $file") unless ( $debug );
	}
	my $in;
	foreach my $chr (@chrs) {
		unless ( $debug ) {
		$ofile = File::Spec->catfile( $tmp_path,
			$chr . "_OP_" . $fm->{'filename'} . ".bam" );
		}else {
			$ofile = File::Spec->catfile( $tmp_path,'FAKE_DEBUG_'.
			$chr . "_OP_" . $fm->{'filename'} . ".bam" );
		}
		if ( !-f $ofile ) {
			$cmd = "samtools view -b $file $chr > $ofile";
			print "cmd run (" . DateTime->now()->time() . "): " . $cmd . "\n";
			push( @ofiles, $ofile );
			## at the moment I need to restart the tool quite often and I do not want to re-create all from scratch.
			unless ( $debug ){
				system($cmd );
			}else {
				system( "touch $ofile");
			}
		}
	}

	$pm->finish;               # Terminates the child process
}

$pm->wait_all_children;

MERGES:
foreach my $chr (@chrs) {
	my $pid = $pm->start and next MERGES;
	my ( $cmd, $ofile, $ifiles );
	$ifiles = File::Spec->catfile( $tmp_path, $chr . "_OP_*.bam" );
	$cmd =
	    "samtools merge "
	  . File::Spec->catfile( $outpath, $chr ."_". $sampleID. ".sorted.bam" )
	  . " $ifiles";
	print "cmd run (" . DateTime->now()->time() . "): " . $cmd . "\n";
	if ( !$debug ){
		system($cmd ) 
	}
	$pm->finish; 
}

$pm->wait_all_children;

## clean up
foreach my $file (@ofiles) {
	unlink($file);
}
my $end = DateTime->now();
print "(" . $end->time() . "):Finished\n";
print "Duration: " .  join(":",$end->subtract_datetime($start)->in_units('days', 'hours', 'seconds'))  . "\n";

