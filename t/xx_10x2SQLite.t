#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 5;
use stefans_libs::flexible_data_structures::data_table;
use stefans_libs::database::Chromium_SingleCell::datavalues;


use File::Spec;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $path, $outfile, );

my $exec = $plugin_path . "/../bin/10x2SQLite.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/10x2SQLite";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$path = File::Spec->catdir($plugin_path, 'data');
$outfile = File::Spec->catfile($outpath,"sqlite3_test.db" );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -path " . $path 
. " -outfile " . $outfile 
. " -debug";
my $start = time;
system( $cmd. " > /dev/null" );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";

ok ( -f $outfile ,' outfile exists' );

my $obj = stefans_libs::database::Chromium_SingleCell::datavalues -> new({ 'file' =>$outfile } );

$value = $obj->get_data_table_4_search( {
	'search_columns' => ['datavalues.id','sname', 'gname', 'value'],
	'where' => [ ]
});

#warn $obj->{'complex_search'}."\n";
ok( $value->Lines() == 29997, "all lines stored");
 
#warn $value->AsTestString();

## double checked using sqlite3
$exp = [ '1', 'AAACCTGAGAAACCGC-1', 'ENSMUSG00000025366', '1' ];
is_deeply( $value->{'data'}[0] , $exp, "data(1) == 1	AAACCTGAGAAACCGC-1	ENSMUSG00000025366	1");

$exp = [ '29993', 'AAAGTAGAGTGACATA-1', 'ENSMUSG00000032328', '1' ];
is_deeply( $value->{'data'}[29992] , $exp, "data(29992) == 29993	AAAGTAGAGTGACATA-1	ENSMUSG00000032328	1");


#print "\$exp = ".root->print_perl_var_def($value->{'data'}[0] ).";\n";


#print "\$exp = ".root->print_perl_var_def($value->{'data'}[29992] ).";\n";



#$obj->getArray_of_Array_for_search
