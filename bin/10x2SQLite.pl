#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-01-17 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit 2b3dac48125abf3cf3d1c692eb96a10f20faf457
   

=head1  SYNOPSIS

    10x2SQLite.pl
       -path       :<please add some info!>
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Convert the 10X test sparse matrix files into a usable SQLite database

  To get further help use '10x2SQLite.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::database::Chromium_SingleCell::datavalues;
use stefans_libs::database::SQLiteBatch;
use File::Spec;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $path, $outfile);

Getopt::Long::GetOptions(
	 "-path=s"    => \$path,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $path) {
	$error .= "the cmd line switch -path is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/10x2SQLite.pl';
$task_description .= " -path '$path'" if (defined $path);
$task_description .= " -outfile '$outfile'" if (defined $outfile);



use stefans_libs::Version;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version '.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $obj = stefans_libs::database::Chromium_SingleCell::datavalues -> new({ 'file' =>$outfile } );

my $batch = stefans_libs::database::SQLiteBatch->new();
my $file;

foreach $file ( 'matrix.mtx', 'genes.tsv', 'barcodes.tsv' ) {
	unless ( -f File::Spec->catfile( $path, $file ) ) {
		$error .= "File '$file' missing from the path $path\n";
	}
}
if ( $error =~ m/\w/ ) {
	die $error;
}




my $genes = data_table->new();
$genes->Add_2_Header(['gname', 'gene.name']);
$genes -> read_file( File::Spec->catfile( $path, 'genes.tsv' ), undef, 1);
$genes->add_column('id', [1..$genes->Rows()]);
$genes->define_subset( 'ok', ['id', 'gname'] );
$genes = $genes -> GetAsObject('ok');
$batch -> batch_import ( $obj->{'genes'}, $genes );


my $samples = data_table->new();
$samples -> Add_2_Header( ['sname']);

$samples -> read_file( File::Spec->catfile( $path, 'barcodes.tsv' ), undef, 1);
$samples -> add_column('id', [1..$samples->Rows()]);
$samples->define_subset( 'ok', ['id', 'sname'] );
$samples = $samples -> GetAsObject('ok');
$batch -> batch_import ( $obj->{'samples'}, $samples );

my $data = data_table->new();
$data -> Add_2_Header( ['id', 'sample_id','gene_id',  'value']);
my $data_id = 1;
open ( IN, "<". File::Spec->catfile( $path, 'matrix.mtx' ) );
my $ok = my $entries = 0;

while ( <IN> ) {
	if ( $_ =~ m/^%/) {
		next;
	}
	unless ( $ok ) {
		$ok = 1;
		next; # skip the dimension info
	}
	chomp;
	push( @{$data->{'data'}}, [$data_id++, (split(/\s+/,$_))[1,0,2] ]);
	$entries ++;
}

$batch -> batch_import ( $obj, $data );

print "$entries values stored in the database '$outfile'\n";


