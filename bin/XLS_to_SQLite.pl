#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-05-23 Stefan Lang

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
   
   binCreate.pl from  commit 
   

=head1  SYNOPSIS

    XLS_to_SQLite.pl
       -infile       :the input tab separated table file
       -outfile      :the sqlite database name
	   -gene_first   :the gene is in the first column not in the last (default)

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Here you can convert a data matrix into a SQLite database of type stefans_libs::database::Chromium_SingleCell::datavalues.

  To get further help use 'XLS_to_SQLite.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::database::Chromium_SingleCell::datavalues;
use stefans_libs::database::SQLiteBatch;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $outfile, $gene_first);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-outfile=s"    => \$outfile,
	 "-gene_first" => \$gene_first,

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

$task_description .= 'perl '.$plugin_path .'/XLS_to_SQLite.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -gene_first" if ( $gene_first );


use stefans_libs::Version;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version'.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $obj = stefans_libs::database::Chromium_SingleCell::datavalues -> new( $outfile );

my $batch = stefans_libs::database::SQLiteBatch->new();

if ( $infile =~ m/.gz$/ ) {
	open ( IN, "zcat $infile |" ) or die "I could not read from '$infile'\n$!\n";
}else {
	open ( IN, "<$infile" ) or die "I could not read from '$infile'\n$!\n";
}
my (@samples,@tmp, $gene_id,$gname, $samples, $genes, $data);
my $entries = 0;

## you might want to re-code that using https://stackoverflow.com/questions/364017/faster-bulk-inserts-in-sqlite3
$gene_id = 1;
$genes = data_table->new();
$genes->Add_2_Header(['id', 'gname']);

$samples = data_table->new();
$samples ->Add_2_Header( ['id', 'sname']);

$data = data_table->new();
$data -> Add_2_Header( ['id', 'sample_id', 'gene_id', 'value']);
my $data_id = 1;

while ( <IN> ) {
	chomp();
	if ( $_ =~ m/^#/ || @samples == 0) { 
		$_ =~ s/^#//;
		@tmp = 	split("\t", $_);
		print "Samples like: ".join("  ", @tmp[1..10])."\n";
		if ( $gene_first ){
			shift(@tmp)
		}else {
			pop(@tmp)
		}
		my $i = 1;
		map { 
			push(@{$samples->{'data'}}, [$i++, $_ ] ) 
		} @tmp;
		$batch -> batch_import ( $obj->{'samples'}, $samples );
		@samples = @{$samples->GetAsArray('id')};
		#$obj -> batch_mode(1);
		next;
	}
	@tmp = 	split("\t", $_);
	if ( $gene_first ){
		$gname = shift(@tmp)
	}else {
		$gname = pop(@tmp)
	}
	push(@{$genes->{'data'}}, [$gene_id, $gname]);
		
	for ( my $i = 0; $i < @tmp; $i++ ) {
		if ( $tmp[$i] =~ m/\d/ ) {
			$entries ++;
			my $a = [$data_id ++, $i+1, $gene_id, $tmp[$i]];
			push( @{$data->{'data'}}, $a ) if ( $tmp[$i] != 0);
			#$obj -> AddDataset( { 'gene_id' => $gene_id, 'sample_id' => $samples[$i], 'value' => $tmp[$i]} );
		}
		#$obj -> commit();
	}
	$gene_id ++;
	print "done with line $gname ( $gene_id genes and a total of $entries expression values )\n" if ($gene_id  % 100 == 0);
}
$batch -> batch_import ( $obj->{'genes'}, $genes );
$batch -> batch_import ( $obj, $data );


print "$entries values stored in the database '$outfile'\n";
