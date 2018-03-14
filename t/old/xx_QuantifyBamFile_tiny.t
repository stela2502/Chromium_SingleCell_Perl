#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 15;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $gtf_file, $infile, $outfile, @options, );

my $exec = $plugin_path . "/../bin/QuantifyBamFile.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/QuantifyBamFile_tiny";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$gtf_file = "$plugin_path/data/GOI.annotation.genecode.v12.gff3";
ok( -f $gtf_file, 'test gtf file exists' );

$infile = "$plugin_path/data/Chromium_hisat2_tiny.bam";
ok( -f $infile, 'test infile exists' );

$outfile = "$outpath/test.xls";

#improvement by using UMIs to merge samples:
#no usage: 61773 samples and 4 genes
#with umi: 25261 and 4 genes
#with umi and sample merge: 16224 samples and 4 genes detected

# with only 4 geens we are down to only half the samples.
@options = ( 'min_UMIs', 10 );

my $cmd = "perl -I $plugin_path/../lib " 
  . " -I  $plugin_path/../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib "
  . " -I  $plugin_path/../../stefans_libs-BAMfile/lib "
  . " -I  $plugin_path/../../stefans_libs-GenomeDB/lib "
  . " $exec "
  . " -gtf_file "
  . $gtf_file
  . " -infile "
  . $infile
  . " -outfile "
  . $outfile
  . " -options "
  . join( ' ', @options )
  . " -debug";

my $start = time;
print $cmd."\n";
system($cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";

#Finished with mapping: 5613 samples and 2 gene_id's detected
#Execution time: 3 s

#when dropping samples with not a single match to the transcriptome
#Finished with mapping: 2123 samples and 2 gene_id's detected
#Execution time: 3 s

#Finished with mapping: 445 samples and 2 gene_id's detected
#Saving data to '/home/stefanl/git/Chromium_SingleCell_Perl/t/data/output/QuantifyBamFile_tiny/test.xls' and '/home/stefanl/git/Chromium_SingleCell_Perl/t/data/output/QuantifyBamFile_tiny/test.spliced.xls' took: 0 hours, 0 min and 0 seconds s
#Total run took: 0 hours, 0 min and 3 seconds s
#Done
#Execution time: 3 s

@values = ();
foreach (
	'test.samples.xls',         'test.xls',
	'test.spliced.xls',         'test.merge.log',
	'test.original.xls',        'test.xls.log',
	'test.original_merged.xls', 'test.merge_report.txt'
  )
{
	ok( -f "$outpath/$_", "outfile $_" );
	push( @values, "$outpath/$_" );
}
system( "wc " . join( " ", @values ) );

open( RSC, ">$outpath/load.R" )
  or die "Could not create the sample R load script.\n$!\n";
print RSC "library(StefansExpressionSet)\n"
  . "data <- read.delim('$outpath/test.xls' )\n"
  . "samples <- read.delim('$outpath/test.samples.xls' )\n"
  . " t <- SingleCellsNGS( dat = data, Samples=samples, namecol='X.sample.tag')\n";

close(RSC);

my $data = data_table->new( { 'filename' => "$outpath/test.xls" } );
ok( $data->Columns == 426,
	"the right column count (" . $data->Columns() . ")" );
ok( $data->Rows == 2, "the right row count (2) " );

my $samples = data_table->new( { 'filename' => "$outpath/test.samples.xls" } );

ok(
	$data->Columns - 1 == $samples->Rows(),
	"the sample table has one row for each data column ("
	  . ( $data->Columns - 1 ) . " == "
	  . $samples->Rows() . ")"
);

my $OK = 1;
for ( my $i = 0 ; $i < $data->Columns() - 1 ; $i++ ) {
	@{ @{ $data->{'data'} }[0] }[$i] ||= 0;
	@{ @{ $data->{'data'} }[1] }[$i] ||= 0;
	unless (
		@{ @{ $data->{'data'} }[0] }[$i] + @{ @{ $data->{'data'} }[1] }[$i] >=
		10 )
	{
		$OK = 0;
		warn "Sample @{$data->{'header'}}[$i] has only "
		  . ( @{ @{ $data->{'data'} }[0] }[$i] +
			  @{ @{ $data->{'data'} }[1] }[$i] )
		  . " reads\n";
	}
}
ok( $OK, "all cells have more than 10 reads in the 2 genes\n" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
