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
use stefans_libs::SLURM;
use DateTime;

use stefans_libs::file_readers::gtf_file;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help,     $debug,    $database, @bams, $outpath,
	$coverage, $tmp_path, $sampleID, $gtf,  $n
);

Getopt::Long::GetOptions(
	"-bams=s{,}"  => \@bams,
	"-outpath=s"  => \$outpath,
	"-coverage=s" => \$coverage,
	"-sampleID=s" => \$sampleID,
	"-gtf=s"      => \$gtf,
	"-n=s"        => \$n,
	"-tmp_path=s" => \$tmp_path,

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
$task_description .= " -outpath '$outpath'"   if ( -d $outpath );
$task_description .= " -coverage '$coverage'" if ( -f $coverage );
$task_description .= " -sampleID '$sampleID'" if ( defined $sampleID );
$task_description .= " -n $n";
$task_description .= " -gtf '$gtf'"           if ( -f $gtf );

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
my ( @chrs, $chr_splices, @tmp );

while (<IN>) {
	chomp;
	@tmp = split( "\t", $_ );
	push( @chrs, $tmp[0] );
	$chr_splices->{ $tmp[0] } = ["$tmp[0]:0-$tmp[1]"];
}
close(IN);

my $pm = Parallel::ForkManager->new($n);

my @ofiles;
if ( -f $gtf ) {
	print "I am reading the genome annotation file\n";
	unless ( -f "$outpath/chr_splits.txt" ) {
		my $cmd =
		    "IdentifyGeneFreeRegions.pl"
		  . " -gtf  $gtf"
		  . " -coverage $coverage"
		  . " -outfile $outpath/chr_splits.txt"
		  . " -step 50000000";
		print $cmd;
		system($cmd );
	}
	$chr_splices = undef;
	open( IN, "<$outpath/chr_splits.txt" )
	  or die "I could not open the chr_splits.txt file\n";
	while (<IN>) {
		chomp;
		if ( $_ =~ m/(.+):(\d+)-(\d+)/ ) {
			$chr_splices->{$1} ||= [];
			push( @{ $chr_splices->{$1} }, $_ );
		}
	}
	close(IN);
	print "genome annotation file reading - finished\n";
	@chrs = (keys %$chr_splices);
}

warn "I am going to merge/split the ". scalar(@bams)." bam files into ".scalar(keys %$chr_splices)." chromosomal areas\n";
warn  "tmp folder $tmp_path\n";
sub byFileSize {
	-s $b <=> -s $a;
}

sub FileSizeOrder {
	my @files = @_;
	my $i     = 0;
	my $order = { map { $_ => $i++ } @files };
	my @ret;
	foreach ( sort byFileSize @files ) {
		push( @ret, $order->{$_} );
	}
	return @ret;
}

my $SLURM = stefans_libs::SLURM->new()
  ;    ## just for the $SLURM->get_files_from_path( $path, @matches );

FILES:
foreach my $file ( sort byFileSize @bams ) {
	#warn "Processing file $file\n";
	my $pid;
	unless ( $debug ) {
		$pid = $pm->start and next FILES;
	} 
	my $fm = root->filemap($file);
	my ( $cmd, $ofile );
	unless ( -f "$file.bai" ) {
		print "samtools index $file\n";
		system("samtools index $file") unless ($debug);
	}
	my $in;
	#warn "split over ".scalar(@chrs). " chromosomes\n";
	foreach my $chr (@chrs) {
		unless ( $chr =~ m/^chr/ ) {
			#$chr = "chr$chr";    ## we use --add-chrname in the hisat2 call
		}
		unless ( defined $chr_splices->{$chr} ){
			Carp::confess("ReshuffleBamFiles internal error - \$chr_splices does not contain key $chr\n");
		}
		foreach my $slice ( @{ $chr_splices->{$chr} } ) {
			#warn "processing slice $slice\n";
			unless ($debug) {
				$ofile = File::Spec->catfile( $tmp_path,
					$chr . "_OP_" . $fm->{'filename'} . ".$slice.bam" );
			}
			else {
				$ofile = File::Spec->catfile( $tmp_path,
					    'FAKE_DEBUG_'
					  . $chr . "_OP_"
					  . $fm->{'filename'}
					  . ".$slice.bam" );
			}
			unless ( -f $ofile ) {    ## the final outfile exists

				$cmd = "samtools view -b $file $slice > $ofile";
				#warn "starting $cmd\n";
				push( @ofiles, $ofile );
				## at the moment I need to restart the tool quite often and I do not want to re-create all from scratch.
				unless ($debug) {
					print "creating $slice\n$cmd\n";
					system($cmd );
				}
				else {
					system("touch $ofile");
					print "cmd run ("
					  . DateTime->now()->time() . "): "
					  . $cmd . "\n";
				}
			}
			else {
				warn "ReshuffleBamFile outfile $ofile exists\n";
			}
		}
	}

	unless ( $debug ) {
		$pm->finish;    # Terminates the child process
	}
}
unless ( $debug ) {
	$pm->wait_all_children;
}

MERGES:
foreach my $chr (@chrs) {
	my $pid = $pm->start and next MERGES;

	#unless ( $chr =~ m/^chr/ ) {
	#	$chr = "chr$chr";    ## we use --add-chrname in the hisat2 call
	#}
	
	unless ( defined $chr_splices->{$chr} ){
		Carp::confess("ReshuffleBamFiles internal error - \$chr_splices does not contain key $chr\n");
	}
	foreach my $slice ( @{ $chr_splices->{$chr} } ) {
		next
		  if ( -f $outpath . $chr . "_" . $sampleID . ".$slice.sorted.*bam" )
		  ;                  ## the final outfile exists

		my ( $cmd, $ofile, $ifiles, @tmp );
		$ofile =
		  File::Spec->catfile( $outpath,
			$chr . "_" . $sampleID . ".$slice.sorted.bam" );

		if (
			scalar(
				@tmp = $SLURM->get_files_from_path(
					$outpath, $chr . "_OP_*$slice.*bam"
				)
			) == 1
		  )
		{
			$cmd = "mv $tmp[0] $ofile\n";
		}
		else {
			$ifiles =
			  File::Spec->catfile( $tmp_path, $chr . "_OP_*$slice.bam" );

			$cmd = "samtools merge " . $ofile . " $ifiles\n";
			$cmd .= "if  [ -f $ofile ] \&\&[ -s $ofile ]; ";
			$cmd .= "then" . "\nrm -f $ifiles \nfi\n";

		}

		print "cmd run (" . DateTime->now()->time() . "): " . $cmd . "\n";
		if ( !$debug ) {
			system($cmd );
		}

	}
	$pm->finish;
}

$pm->wait_all_children;

## clean up
foreach my $file (@ofiles) {
	unlink($file);
}

## now I have some that are extremely big - they slow down the whole process..

my $end = DateTime->now();
print "(" . $end->time() . "):Finished\n";

print "Duration: "
  . join( ":",
	$end->subtract_datetime($start)->in_units( 'days', 'hours', 'seconds' ) )
  . "\n";

sub get_files_from_path {
	my ( $path, @matches ) = @_;
	return $SLURM->get_files_from_path( $path, @matches );
}

