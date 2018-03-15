#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $coverage, $step, $gtf, $outfile, );

my $exec = $plugin_path . "/../bin/IdentifyGeneFreeRegions.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/IdentifyGeneFreeRegions";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$gtf = $plugin_path . "/data/ReshuffleBamFiles/genes.gtf.gz";
$step = 1e+7;
$coverage =$plugin_path . "/data/ReshuffleBamFiles/GRCm38.p5_chrom_sizes.txt";

$outfile = "$outpath/chr_splits.txt";
my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -coverage " . $coverage 
. " -step " . $step 
. " -gtf " . $gtf 
. " -outfile " . $outfile 
#. " -debug"
. " -gtf "
;

print "$cmd\n";
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";


## run this once to check the general function
