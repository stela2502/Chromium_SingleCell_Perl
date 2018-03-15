#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-03-15 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit 70d6934e575eaa13cb8390e5c050eeb0d38bba81
   

=head1  SYNOPSIS

    IdentifyGeneFreeRegions.pl
       -gtf       :the gtf file defining the gene areas that should not be touched
       -outfile   :the outfile containing a list of chr regions USCS format
       -coverage  :the chromosome length file 
       -step      :the step size in bp (defualt 1e+7)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  A script returning a list of chromosome areas that can be used to split a chromosme into pieces without splitting a gene in 2 parts.

  To get further help use 'IdentifyGeneFreeRegions.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::file_readers::gtf_file;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, $gtf, $outfile, $coverage, $step );

Getopt::Long::GetOptions(
	"-gtf=s"      => \$gtf,
	"-outfile=s"  => \$outfile,
	"-coverage=s" => \$coverage,
	"-step=s"     => \$step,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( -f $gtf ) {
	$error .= "the cmd line switch -gtf is undefined!\n";
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ( -f $coverage ) {
	$error .= "the cmd line switch -coverage is undefined!\n";
}
unless ( defined $step ) {
	$step = 1e+7;
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

$task_description .= 'perl ' . $plugin_path . '/IdentifyGeneFreeRegions.pl';
$task_description .= " -gtf '$gtf'" if ( defined $gtf );
$task_description .= " -outfile '$outfile'" if ( defined $outfile );
$task_description .= " -coverage '$coverage'" if ( defined $coverage );
$task_description .= " -step '$step'" if ( defined $step );

use stefans_libs::Version;
my $V  = stefans_libs::Version->new();
my $fm = root->filemap($outfile);
mkdir( $fm->{'path'} ) unless ( -d $fm->{'path'} );

open( LOG, ">$outfile.log" ) or die $!;
print LOG '#library version ' . $V->version('Chromium_SingleCell_Perl') . "\n";
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!

## read the gtf file
my $gtf_obj = stefans_libs::file_readers::gtf_file->new();

my $filter = sub {
	my ( $self, $array, $i ) = @_;
	return ( @$array[2] eq "gene" );
};

$gtf_obj->read_file( $gtf, $filter );

print "gtf file read (".$gtf_obj->Lines()." entries)\n";

## read the chromosome length file
my ( $chr_info, @tmp );
open( IN, "<$coverage" )
  or die "I could not open the coverage file '$coverage'\n$!\n";
while (<IN>) {
	chomp;
	@tmp = split( /\s+/, $_ );
	$chr_info->{ $tmp[0] } = $tmp[1];
}
close(IN);
print "I got info on " . scalar( keys %$chr_info ) . " chromosomes\n";

my ( $last, $chr, $max_pos, @result );

while ( ( $chr, $max_pos ) = each(%$chr_info) ) {
	my $next_start = $last = 0;
	print "Woring on chr $chr\n";
  POSITION:
	for ( my $i = $step ; $i < $max_pos + $step ; $i += $step ) {
		#print "gene in area $chr:$i - test 1.\n";
		$next_start = &goOn( $i, $max_pos );
		unless ( defined $next_start ) {
			#print "\nposition $i was empty!\n";
			next POSITION;
		}
		#print "gene in area $chr:$i - try next empty area.\n";
		while ( defined $next_start ) {
			#print "gene in area $chr:$next_start - try next empty area.\n";
			$next_start = &goOn( $next_start + 51, $max_pos );
		}
		#print "position $last was empty!\n";
	}
}

open ( OUT , ">$outfile" ) or die "I could not create the outfile '$outfile'\n$!\n";
print OUT join("\n",@result );
close ( OUT );

print "The outfile contains ".scalar(@result)." chromosome positions to split on.\n";

sub goOn {
	my ( $pos, $max_pos ) = @_;
	my @ids =
	  $gtf_obj->efficient_match_chr_position( $chr, $pos - 50, $pos + 50 );
	if ( scalar(@ids) == 0 ) {
		## take that?
		$pos = $max_pos if ( $pos > $max_pos );
		push( @result, "$chr:$last-$pos" );
		print "\npushed $chr:$last-$pos\n";
		$last = $pos;
		return undef;
	}
	my $r = &max( &get_chr_end_4_ids( $gtf_obj, @ids ) )
	  ;    ## the most likely next empty place.
	if ( $r == 0 ) {
		## what a crap"!
		die "I have the ids "
		  . join( ", ", @ids )
		  . " translating to starts at "
		  . join( ", ", &get_chr_end_4_ids( $gtf_obj, @ids ) )
		  . " an I still get a max as $r??\n";
	}
	return $r;
}

sub get_chr_end_4_ids {
	my ( $gtf, @lines ) = @_;
	return unless ( scalar(@lines) );
	return map {
		if ( defined $_ ) {
			print "mapper line identified $_: "
			  . join( ";", @{ @{ $gtf_obj->{'data'} }[$_] } ) . "\n";
			@{ @{ $gtf_obj->{'data'} }[$_] }[4];
		}
	} @lines;
}

sub max {
	my $max = 0;
	map { $max = $_ if ( $max < $_ ) } @_;
	return $max;
}

