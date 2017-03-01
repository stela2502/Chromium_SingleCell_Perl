#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 10;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $gtf_file, $infile, $outfile, @options, );

my $exec = $plugin_path . "/../bin/QuantifyBamFile.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/QuantifyBamFile";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$gtf_file = "$plugin_path/data/GOI.annotation.genecode.v12.gff3";
ok ( -f $gtf_file, 'test gtf file exists');

$infile = "$plugin_path/data/Chromium_hisat2.bam";
ok ( -f $infile, 'test infile exists');

$outfile = "$outpath/test.xls";

#improvement by using UMIs to merge samples:
#no usage: 61773 samples and 4 genes
#with umi: 25261 and 4 genes
#with umi and sample merge: 16224 samples and 4 genes detected

# with only 4 geens we are down to only half the samples.

@options = ('min_UMIs', 10 );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -gtf_file " . $gtf_file 
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -options " . join(' ', @options )
. " -debug"
;

my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";


#Finished with mapping: 25261 samples and 4 gene_id's detected
#Execution time: 18 s

## when dropping samples with no hit to the transcriptome
#Finished with mapping: 9000 samples and 4 gene_id's detected
#Execution time: 19 s

## when dropping samples with less than 100 hits to the transcriptome
#Finished with mapping: 15 samples and 4 gene_id's detected
#Done
#Execution time: 25 s

## After using PDL for the computational stuff and min_umi 100
##Finished with mapping: 161 samples and 4 gene_id's detected
##Done
##Execution time: 17 s

## After using PDL for the computational stuff and min_umi 10
##Finished with mapping: 646 samples and 4 gene_id's detected
##Done
##Execution time: 19 s

@values = undef;
foreach ('test.samples.xls', 'test.xls', 'test.spliced.xls','test.merge.log','test.original.xls', 'test.xls.log',
'test.original_merged.xls'
) {
	ok ( -f "$outpath/$_", "outfile $_");
	push ( @values , "$outpath/$_");
}
system( "wc ". join(" ", @values));



#print "\$exp = ".root->print_perl_var_def($value ).";\n";