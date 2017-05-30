#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outfile, $infile, );

my $exec = $plugin_path . "/../bin/XLS_to_SQLite.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/XLS_to_SQLite";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
}else {
	system("mkdir -p $outpath" );
}

$infile  = $plugin_path."/data/HSC_data.xls.gz";
$outfile = "$outpath/test.db";


my $cmd =
    "perl -I $plugin_path/../lib -I  $plugin_path/../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib $exec "
. " -outfile " . $outfile 
. " -infile " . $infile 
. " -debug";
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";