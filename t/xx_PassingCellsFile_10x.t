#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 4;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $infile, $outfile, );

my $exec = $plugin_path . "/../bin/PassingCellsFile_10x.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/PassingCellsFile_10x";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$infile  = "$plugin_path/data/Spliced_Reads.bam";
$outfile = "$outpath/outfile.xls";

ok( -f $infile, 'bam infile' );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . " -infile "
  . $infile
  . " -outfile "
  . $outfile
  . " -debug";
my $start = time;
system($cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";

ok( -f $outfile, "outfile created" );

open( IN, "<$outfile" ) or die $!;
@values = map { chomp; $_ } <IN>;
close(IN);
#print "\$exp = " . root->print_perl_var_def( \@values ) . ";\n";
$exp = [
	'#total of 106 reads processed',
	'#no sample/UMI information in 0',
	'Cell_ID	count',
	'AAACCTGAGAGGACGG	1',
	'AAAGCAAAGGAGCGTT	1',
	'AACACGTCATAAAGGT	5',
	'ACCTTTATCTCAACTT	5',
	'ACGAGCCGTACGCACC	5',
	'ACGGAGAAGAGTGACC	1',
	'ACTGAGTCACAGGTTT	1',
	'AGCGTATGTTCCACAA	2',
	'AGTGGGATCTTTACAC	4',
	'AGTTGGTAGGACGAAA	2',
	'ATAGACCGTGGTGTAG	1',
	'ATTACTCCATCGGAAG	1',
	'ATTTCTGTCTAGAGTC	1',
	'CAAGTTGTCGCATGGC	1',
	'CACAAACCACTTCTGC	1',
	'CACCACTCTTGACGTT	1',
	'CATCAAGGTGCAGTAG	2',
	'CATCAGAGTGTGTGCC	1',
	'CATTCGCGTAGCAAAT	1',
	'CCACTACCAGTATGCT	2',
	'CCGTGGACATAACCTG	2',
	'CCTCAGTGTTGACGTT	1',
	'CCTTTCTCAAACGTGG	1',
	'CGAACATTCATTCACT	1',
	'CGCTTCAAGGTGCACA	1',
	'CGGACTGAGCTGAAAT	2',
	'CGGACTGTCCTGCTTG	1',
	'CGGTTAAGTGCACGAA	1',
	'CGTTAGACACGGTTTA	4',
	'CGTTCTGAGGTGCACA	1',
	'CGTTCTGTCGACCAGC	1',
	'CTCGAGGCATCACGAT	1',
	'CTCTTCCGATCTCTCT	1',
	'CTGAAGTGTGGGTATG	1',
	'GACACGCCAGCTCGAC	1',
	'GACGCGTCAGCTTAAC	3',
	'GATCGATAGCCACTAT	3',
	'GATGAGGCAACACGCC	1',
	'GATTCAGGTTATTCTC	1',
	'GCAGTTATCAGCTCTC	1',
	'GCGAGAACATGCCTTC	1',
	'GCGCCAATCAAAGTAG	1',
	'GCGCCAATCTTGTACT	5',
	'GGACCTTCATGCATGC	1',
	'GGGAATGAGAACAATC	1',
	'GGGACCTCATCCCACT	1',
	'GGGCACTCACGCTTTC	3',
	'GTCATTTAGTCGTTTG	1',
	'GTCCTCATCTGCAGTA	1',
	'GTGCATAGTCTAGAGG	1',
	'GTTACAGTCTTTAGGG	1',
	'GTTCATTTCCAAGTAC	1',
	'TAAGCGTGTGGGTATG	3',
	'TACAGTGGTTATTCTC	1',
	'TCAATCTCATTTGCTT	1',
	'TCATTTGAGACTACAA	1',
	'TCGCGAGGTGAACCTT	1',
	'TCGCGAGTCCCATTAT	1',
	'TGATTTCCAAGCCATT	1',
	'TGCGTGGAGCTACCGC	3',
	'TGCTGCTCAAACCTAC	5',
	'TTCTCAATCAAGGTAA	2',
	'TTTGTCAGTAGCTGCC	1'
];

is_deeply( \@values, $exp, "right cell counts");

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

