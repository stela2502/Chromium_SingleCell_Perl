#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 4;
BEGIN { use_ok 'stefans_libs::database::Chromium_SingleCell::datavalues' }

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $path );

$path = "$plugin_path/data/output/Chromium_SingleCelldb/";
if ( -d $path ) {
	system( "rm $path*" );
}else {
	system( "mkdir -p $path" );
}

my $obj = stefans_libs::database::Chromium_SingleCell::datavalues -> new( { file=> $path."database.db"} );

is_deeply ( ref($obj) , 'stefans_libs::database::Chromium_SingleCell::datavalues', 'simple test of function stefans_libs::database::Chromium_SingleCell::datavalues -> new()' );

$value = $obj->insert( { 'value' => 100, 'gene' => "Totally uninteresting gene", "sample" =>  "also extremely uninteresting"  });

ok( $value == 1, "we added one entry ($value)" );

$value = $obj->get_data_table_4_search (
{
	'search_columns' => ['gname', 'sname', 'value'],
	'where' => [ ],
}
);

#print "\$exp = ".root->print_perl_var_def( [split( /[\t\n]/ ,$value->AsString())] ).";\n";
$exp = [ '#gname', 'sname', 'value', 'Totally uninteresting gene', 'also extremely uninteresting', '100' ];

is_deeply ( $exp, [split( /[\t\n]/ ,$value->AsString())], "Data in == data out" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";




