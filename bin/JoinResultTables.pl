#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-02-27 Stefan Lang

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

    JoinResultTables.pl
       -paths       :the output paths from either a list of cellranger run or QuantifyBamFile.pl runs
       -outfile     :the sqlite database filename
       -createTable :option to create not only the databse, but also the summary tables


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  combine all result tables into one sqlite database.

  To get further help use 'JoinResultTables.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::result_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @paths, $outfile, $createTable);

Getopt::Long::GetOptions(
       "-paths=s{,}"    => \@paths,
	 "-outfile=s"    => \$outfile,
       "-createTable"    => \$createTable,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $paths[0]) {
	$error .= "the cmd line switch -paths is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
# createTable - no checks necessary


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

$task_description .= 'perl '.$plugin_path .'/JoinResultTables.pl';
$task_description .= ' -paths "'.join( '" "', @paths ).'"' if ( defined $paths[0]);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -createTable " if ( $createTable);



use stefans_libs::Version;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version '.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
$outfile .= ".sqlite" unless ( $outfile =~m/\.sqlite$/);

my $result_table = stefans_libs::result_table ->new( { 'filename' => $outfile } );

## There might be the possibility, that these paths do not contain the info we need.
## As this is not supported by the object I need to check the paths here.

for ( my $i = $#paths; $i >= 0; $i-- ) {
	unless ( -f "$paths[$i]/matrix.mtx" ) {
		warn "path $paths[$i] does not contain the required matrix.mtx file - skipped\n";
		splice(@paths,$i,1);
	}
}
$result_table->import_tables(@paths);

if ( $result_table->{'data_storage_spliced'}){
	print "Data is stored in with splice information";
}

$result_table->write_table( $outfile );

if ( $createTable ) {
	$result_table->print2table ( $outfile );
}

print "Finished\n";
