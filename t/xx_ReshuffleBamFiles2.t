#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 289;
use stefans_libs::flexible_data_structures::data_table;
use File::Spec::Functions;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @bams, $outpath, $coverage, @outfiles);

my $exec = $plugin_path . "/../bin/ReshuffleBamFiles.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/ReshuffleBamFiles";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
}

my $tmp_path =  File::Spec->catfile($outpath, 'tmp');

my $data_path = File::Spec->catfile( $plugin_path, 'data','ReshuffleBamFiles' );

ok (-f "$data_path/exprected_outfiles_using_gtf.txt" );
open ( IN, "<$data_path/exprected_outfiles_using_gtf.txt" ) or die "package not complete - important infile missing:$data_path/exprected_outfiles_no_gtf.txt\n$!\n";
while( <IN> ) {
	chomp;
	push( @outfiles, $_)
}
close ( IN );

@bams = map{File::Spec->catfile($data_path, 'Sample'.$_.'.bam') } 1,2 ;

unless ( -d $data_path) {
	Carp::confess ( "Package error: test files are missing '$data_path'\n" );
}

$coverage = "$data_path/GRCm38.p5_chrom_sizes.txt";
my $gtf = "$data_path/genes.gtf.gz";

my $sampleName = 'test';

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -bams " . join(' ', @bams )
. " -outpath " . $outpath 
. " -coverage " . $coverage 
. " -sampleName $sampleName"
. " -gtf $gtf"
. " -n 4"
. " -tmp_path $tmp_path"
. " -sampleID TEST"
#. " -debug"
;

system("rm -Rf $outpath/*");

my $start = time;
system( $cmd );
my $duration = time - $start;
print "Run time: $duration s\n";


foreach ( @outfiles ) {
	ok ( -f "$outpath/$_", "outfile $_" );
}

## This does all seam to work - I have no time to check in detail!