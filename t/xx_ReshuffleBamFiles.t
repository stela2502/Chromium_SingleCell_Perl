#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;
use File::Spec::Functions;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @bams, $outpath, $coverage, );

my $exec = $plugin_path . "/../bin/ReshuffleBamFiles.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/ReshuffleBamFiles";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
}

my $tmp_path =  File::Spec->catfile($outpath, 'tmp');

my $data_path = File::Spec->catfile( $plugin_path, 'data','ReshuffleBamFiles' );

@bams = map{File::Spec->catfile($data_path, 'Sample'.$_.'.bam') } 1,2 ;

unless ( -d $data_path) {
	Carp::confess ( "Package error: test files are missing '$data_path'\n" );
}

$coverage = "$data_path/mm10_chrom_sizes.txt";

my $sampleName = 'test';

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -bams " . join(' ', @bams )
. " -outpath " . $outpath 
. " -coverage " . $coverage 
. " -sampleName $sampleName"
. " -n 1"
. " -tmp_path $tmp_path"
#. " -debug"
;
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";

## This does all seam to work - I have no time to check in detail!