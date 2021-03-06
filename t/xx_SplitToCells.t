#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 8;
use stefans_libs::flexible_data_structures::data_table;
use stefans_libs::FastqFile;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, $I1, $R2, $outpath, $R1, );

my $exec = $plugin_path . "/../bin/SplitToCells.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/SplitToCells";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}
$R1 = "$plugin_path/data/R1.fastq.gz";
ok( -f $R1, $R1 );
$R2 = "$plugin_path/data/R2.fastq.gz";
ok( -f $R2, $R2 );
$I1 = "$plugin_path/data/I1.fastq.gz";
ok( -f $I1, $I1 );

@options = ( 'oname', 'test_analysis' );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . " -options "
  . join( ' ', @options ) . " -I1 "
  . $I1 . " -R2 "
  . $R2
  . " -outpath "
  . $outpath . " -R1 "
  . $R1
  . " -debug";
system($cmd );

my $tester = stefans_libs::FastqFile->new();

my $reporter = sub {
	my ( $fastqfile, $entry ) = @_;
	push( @values, $entry->copy() );
	$entry->clear();
};
my $outfile = "$outpath/test_analysis.annotated.fastq.gz";
ok( -f $outfile, "fastq outfile created" );

$tester->filter_multiple_files( $reporter, $outfile );

$value = { map { $_->Get_UMI_Tag() => $_->sequence() } @values };

#print "\$exp = " . root->print_perl_var_def($value) . ";\n";

$exp = {
	'S_ATGAATCT_C_AGCGTCGGTTCAACCA_TTCATGCACG' =>
'ATTGACTCCTTCATCCTAGGAGGAAGTGGATCTTGATTGTACACATATGACTTCTCAACTAACGATGATATAGAATAAATAATCAAATCCA',
	'S_ATGAATCT_C_CCGTGGAGTGTTGAGG_CGCCCGGAAA' =>
	  'CCACTAGTTACCACCACTGTAGGGGCTGCCCCATCTTCCCCCCTAACTGTTCGTAACTA',
	'S_ATGAATCT_C_CCTTCGAGTCTAGCGC_AGGATAAGGT' =>
'GTGCACTGCGGGTGAGCAGCCGGGCCCGGAGTCGCATCCTCAAGGCTGGGGGTAAGATCCTCACCTGTGACCAGCTGGCCCTGGAGT',
	'S_ATGAATCT_C_GAGCAGAGTAAATGTG_ATACTCCTTA' =>
'CTGTGCTGATCACACGGCAATACAGTACCTTTTCCTCTGTTCTCTGTCTGTGTCCTGTTATGCTTCTTTGGGGATAGGATGAGCAGCCACCATT',
	'S_ATGAATCT_C_GCGCGATGTCTTGTCC_AGCCTGCCTC' =>
'CTTTGTGATCTTAAGATAATCAAGAAATGCACAGTAACTGATGTAACTTGGGCTGCAGTCTGTAATCCACCTCTTCAGCTATAAGTAGGGGACATGT',
	'S_ATGAATCT_C_GGACATTGTTGGGACA_GCGCCTGCCA' =>
'GTTGCCCAGTTCCCCAGTTCTGGCAGATTCTGGAAACTGTACATTGAAGCAGAGGTTAATAGTTTATTTATTTTTTCTTATATAGCATCTGCTG',
	'S_ATGAATCT_C_GTCTCGTAGCGTCAAG_CGTACCGTTC' =>
'GTCTCGTAGCCTCAAGAGTACCGTTAAAAAATAAATTTTATCTATAAAGTAACTGAAAGGAAGTAAGACATTAGGTAACTCTGAAGAATAAATCAGA',
	'S_ATGAATCT_C_GTTCGGGGTCAATACC_GCCCGATCTG' =>
'CACAGCTAGTGGCAAGGATAGTGGGTCCTTGTATATTCTGCCTTCTTACCCCAATACCTGGGCAAGGAACCTAGGGTGAGTTGGGGGAAAATGAATCA',
	'S_ATGAATCT_C_TGAGGGAAGTCAAGGC_ACGAAGTTAG' =>
'GCTGCGTCAGCACCTGGAGAGGCTGAAGAAAATTCGAGCCCATAGAGGGCTGCGCCACATTTCGCGCCTTCGTCTCCTCGGTCAGC',
	'S_ATGAATCT_C_TGTTCCGGTAAATGTG_CTCCTGGGGC' =>
'GACAATCAGCTGAGCATAGAGGAGAGTATGATGCTGATGGCAAAGTTGATCTATGCCTGTCATGAGAAGCTGCATGAGAACAACCCACGTCCGCAT'
};
is_deeply( $value, $exp, "sequence to sample info OK" );

#print $cmd. "\n";

ok( -f "$outpath/test_analysis.per_cell_read_count.xls",
	'per cell read count created' );
use stefans_libs::flexible_data_structures::data_table;

my $data_table = data_table->new( { 'no_doubble_cross' => 1 } );
my $file = "$outpath/test_analysis.per_cell_read_count.xls";
my ( $reads, @tmp );
$value = undef;
open( IN, "<$file" ) or die "I could not open file ''\n$!\n";
while (<IN>) {
	next if $_ =~ m/^#/;
	chomp;
	@tmp = split( "\t", $_ );
	next if ( $tmp[1] eq "count" );
	if ( $tmp[1] >= 1 ) {
		#$reads = $tmp[1];
		#@tmp = split( "_", $tmp[0] );
		$value->{ $tmp[0] } += $tmp[1];
	}
}
close(IN);

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
$exp = {
	'S_ATGAATCT_C_AGCGTCGGTTCAACCA' => '1',
	'S_ATGAATCT_C_CCGTGGAGTGTTGAGG' => '1',
	'S_ATGAATCT_C_CCTTCGAGTCTAGCGC' => '1',
	'S_ATGAATCT_C_GAGCAGAGTAAATGTG' => '1',
	'S_ATGAATCT_C_GCGCGATGTCTTGTCC' => '1',
	'S_ATGAATCT_C_GGACATTGTTGGGACA' => '1',
	'S_ATGAATCT_C_GTCTCGTAGCGTCAAG' => '1',
	'S_ATGAATCT_C_GTTCGGGGTCAATACC' => '1',
	'S_ATGAATCT_C_TGAGGGAAGTCAAGGC' => '1',
	'S_ATGAATCT_C_TGTTCCGGTAAATGTG' => '1'
};

is_deeply( $value, $exp, "cell_id_to count data" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
