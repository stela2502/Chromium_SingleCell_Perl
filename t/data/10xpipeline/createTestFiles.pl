#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-02-28 Stefan Lang

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

    createTestFiles.pl
       -R1       :<please add some info!>
       -R2       :<please add some info!>
       -I1       :<please add some info!>
       -outpath       :<please add some info!>
       -reads       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  create a 10x test set

  To get further help use 'createTestFiles.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $R1, $R2, $I1, $outpath, $reads);

Getopt::Long::GetOptions(
	 "-R1=s"    => \$R1,
	 "-R2=s"    => \$R2,
	 "-I1=s"    => \$I1,
	 "-outpath=s"    => \$outpath,
	 "-reads=s"    => \$reads,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';


my ( @R1, @R2, @I1);

unless ( defined $R1 ) {
	$error .= "the cmd line switch -R1 is undefined!\n";
}elsif ( $R1 =~m/.txt$/ ){
	open ( IN, "<$R1") or die $!;
	my @tmp;
	while ( <IN> ) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_);
		}
	}
	close ( IN );
	@R1 = @tmp;
}
unless ( defined $R2 ) {
	$error .= "the cmd line switch -R2 is undefined!\n";
}elsif ( $R2 =~m/.txt$/ ){
	open ( IN, "<$R2") or die $!;
	my @tmp;
	while ( <IN> ) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_);
		}
	}
	close ( IN );
	@R2 = @tmp;
}
unless ( defined $I1 ) {
	$error .= "the cmd line switch -I1 is undefined!\n";
}elsif ( $I1 =~m/.txt$/ ){
	open ( IN, "<$I1") or die $!;
	my @tmp;
	while ( <IN> ) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_);
		}
	}
	close ( IN );
	@I1 = @tmp;
}

unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $reads) {
	$error .= "the cmd line switch -reads is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/createTestFiles.pl';
$task_description .= " -R1 '$R1'" if (defined $R1);
$task_description .= " -R2 '$R2'" if (defined $R2);
$task_description .= " -I1 '$I1'" if (defined $I1);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= " -reads '$reads'" if (defined $reads);



mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_createTestFiles.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
## zcat HJ2GYBGX5/outs/fastq_path/YUAN10x/5/Ctrl-LSK_S5_L001_R1_001.fastq.gz | head -n 4000 | gzip -9  > Ctrl-LSK_S5_L001_R1_part.fastq.gz

open ( I1, ">$outpath/I1.txt" ) or die $!;
my ($ofile, $fm, @tmp );
foreach $I1 ( @I1 ) {
	$fm = root->filemap( $I1 ) ;
	@tmp = split( "/", $fm->{'path'} );
	$ofile = "$outpath/$tmp[0]/$fm->{'filename'}";
	mkdir ( "$outpath/$tmp[0]" ) unless ( -d "$outpath/$tmp[0]");
	system( "zcat $I1 | head -n $reads | gzip -9 > $ofile");
	print I1 "$ofile\n";
}
close ( I1 );


open ( R1, ">$outpath/R1.txt" ) or die $!;
my ($ofile, $fm );
foreach $R1 ( @R1 ) {
        $fm = root->filemap( $R1 ) ;
        @tmp = split( "/", $fm->{'path'} );
        $ofile = "$outpath/$tmp[0]/$fm->{'filename'}";
        system( "zcat $R1 | head -n $reads | gzip -9 > $ofile");
        print R1 "$ofile\n";
}
close ( R1 );

open ( R2, ">$outpath/R2.txt" ) or die $!;
my ($ofile, $fm );
foreach $R2 ( @R2 ) {
        $fm = root->filemap( $R2 ) ;
        @tmp = split( "/", $fm->{'path'} );
        $ofile = "$outpath/$tmp[0]/$fm->{'filename'}";
        system( "zcat $R2 | head -n $reads | gzip -9 > $ofile");
        print R2 "$ofile\n";
}
close ( R2 );

print "Done\n";
