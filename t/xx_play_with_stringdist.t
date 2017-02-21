#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use Text::Levenshtein qw(distance);

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, );

my $outpath = "$plugin_path/data/output/play_with_stringdist";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");

}



my $infile = "$plugin_path/data/HCS.per_cell_read_count.xls";
ok ( -f $infile, 'infile present');

my $data_table = data_table->new( {'no_doubble_cross'=> 1});

$data_table->read_file( $infile );

ok ($data_table->Rows() == 730591, $data_table->Rows()." different cell tags" );

## now I would like to collaps the data using a Levenshtein distance in the beginning - I might need to implement a weight based function later on

## this part should probably go into a second run of the tool itself - which are the acceptable tags?





#print "\$exp = ".root->print_perl_var_def($value ).";\n";