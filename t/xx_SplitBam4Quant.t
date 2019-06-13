#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $min_reads, $infile, $outpath, $cell_ids, );

my $exec = $plugin_path . "/../bin/SplitBam4Quant.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/SplitBam4Quant";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -min_reads " . $min_reads 
. " -infile " . $infile 
. " -outpath " . $outpath 
. " -cell_ids " . $cell_ids 
. " -debug";
system( $cmd );
#print "\$exp = ".root->print_perl_var_def($value ).";\n";