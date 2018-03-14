#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 6;
use stefans_libs::flexible_data_structures::data_table;

use stefans_libs::database::Chromium_SingleCell::datavalues;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $gtf_file, $infile, $outfile, @options, );

my $exec = $plugin_path . "/../bin/QuantifyBamFile.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/QuantifyBamFile_CellRanger_Mpo";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$gtf_file = "$plugin_path/data/mm10.Mpo_large.gtf";

#$gtf_file = "/home/stefanl/lunarc/genomes/mouse/mm10/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf";

ok( -f $gtf_file, 'test gtf file exists' );

$infile =
  "$plugin_path/data/mm10.singleCells.CellRanger.chr11_87760834_87856494.bam";
ok( -f $infile, 'test infile exists' );

$outfile = "$outpath/singleCells.CellRanger.chr11_87760834_87856494_quantifed";

#improvement by using UMIs to merge samples:
#no usage: 61773 samples and 4 genes
#with umi: 25261 and 4 genes
#with umi and sample merge: 16224 samples and 4 genes detected

# with only 4 genes we are down to only half the samples.

@options = (
	'min_UMIs', 1

	  #	, 'quantify_on', 'transcript'
);

my $cmd =
    "perl "
  . "-I $plugin_path/../lib "
  . "-I $plugin_path/../../stefans_libs-BAMfile/lib/ "
  . "-I ../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib "
  . "-I ../../SLURM_bioinformatics_scripts/lib/ "
  . "-I ../../SLURM_bioinformatics_scripts/lib/ "
  . "-I ../../stefans_libs-GenomeDB/lib "
  . "$exec "
  . " -gtf_file "
  . $gtf_file
  . " -infile "
  . $infile
  . " -outfile "
  . $outfile
  . " -options "
  . join( ' ', @options )

  #  . " -debug"
  ;

my $start = time;
system($cmd );
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

my $db = stefans_libs::database::Chromium_SingleCell::datavalues->new(
	{
		'file' =>
"$outpath/singleCells.CellRanger.chr11_87760834_87856494_quantifed.original_merged.db",
		'data_storage_spliced' => 1
	}
);

$value = $db->get_data_table_4_search(
	{
		'search_columns' => [ 'sname', 'gname', 'value' ],
		'where'          => [],
	}
);

ok( $value->Rows() == 205, "205 samples(" . $value->Rows() . ")" );

print "\$exp = "
  . root->print_perl_var_def( [ sort keys %{ $value->createIndex('gname') } ] )
  . ";\n";

$exp = [
	'1',                     'ENSMUSG00000009350.13',
	'ENSMUSG00000034121.12', 'ENSMUSG00000034156.16'
];

is_deeply( [ sort keys %{ $value->createIndex('gname') } ],
	$exp, "only one gene ENSMUSG00000022103.9" );

print "\$exp = " . root->print_perl_var_def( $value->{'data'} ) . ";\n";
$exp = [
	[ 'CCAATCCAGCTACCGC', 'ENSMUSG00000034156.16', '4' ],
	[ 'ACCCACTGTTACCAGT', 'ENSMUSG00000034156.16', '2' ],
	[ 'TGGCTGGCACCTCGGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'AGTGTCATCACCGGGT', 'ENSMUSG00000034156.16', '5' ],
	[ 'CGCGGTACAGCCTTGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'TTGTAGGGTTTGCATG', 'ENSMUSG00000034156.16', '2' ],
	[ 'ATTATCCAGCCCAACC', 'ENSMUSG00000034156.16', '3' ],
	[ 'GGCAATTGTACCGCTG', 'ENSMUSG00000034156.16', '3' ],
	[ 'GACAGAGGTACGCTGC', 'ENSMUSG00000034156.16', '5' ],
	[ 'CACTCCAGTTGGTTTG', 'ENSMUSG00000034156.16', '3' ],
	[ 'AAGACCTGTGCACCAC', 'ENSMUSG00000034156.16', '6' ],
	[ 'GCAGTTAAGACGACGT', 'ENSMUSG00000034156.16', '3' ],
	[ 'GCGGCTGTTTGGCGCG', 'ENSMUSG00000034156.16', '2' ],
	[ 'AGGTCCGGTTCACCTC', 'ENSMUSG00000034156.16', '3' ],
	[ 'GCGACCAAGTCCTCCT', 'ENSMUSG00000034156.16', '3' ],
	[ 'AACTGGTTCCACTGGG', 'ENSMUSG00000034156.16', '3' ],
	[ 'ATGGGAGGTCTAGTGT', 'ENSMUSG00000034156.16', '10' ],
	[ 'GGAATAATCACGGTTA', 'ENSMUSG00000034156.16', '3' ],
	[ 'TTCGAAGAGTTGTCGT', 'ENSMUSG00000034156.16', '3' ],
	[ 'TTTGGTTCACACGCTG', 'ENSMUSG00000034156.16', '4' ],
	[ 'GATTCAGAGCTATGCT', 'ENSMUSG00000034156.16', '2' ],
	[ 'ACGGCCAGTTCCTCCA', 'ENSMUSG00000034156.16', '4' ],
	[ 'CATCCACGTAAGTGGC', 'ENSMUSG00000034156.16', '4' ],
	[ 'TGCTACCAGTGTCCAT', 'ENSMUSG00000034156.16', '3' ],
	[ 'TCAGGTATCAGCACAT', 'ENSMUSG00000034156.16', '3' ],
	[ 'AAAGTAGGTACCAGTT', 'ENSMUSG00000034156.16', '16' ],
	[ 'CGAACATCATCCTAGA', 'ENSMUSG00000034156.16', '4' ],
	[ 'TCACGAATCACGCATA', 'ENSMUSG00000034156.16', '6' ],
	[ 'AAGGTTCTCCTCAATT', 'ENSMUSG00000034156.16', '1' ],
	[ 'ACACCGGGTTGTTTGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'CTCCTAGTCTACGAGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'GCTCCTAGTCCGTGAC', 'ENSMUSG00000034156.16', '3' ],
	[ 'CTTGGCTGTTTCCACC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CAGCATAAGTTTCCTT', 'ENSMUSG00000034156.16', '2' ],
	[ 'GATCGTACATCGATTG', 'ENSMUSG00000034156.16', '2' ],
	[ 'ACATGGTAGATCTGCT', 'ENSMUSG00000034156.16', '5' ],
	[ 'CATCGAACATGCTGGC', 'ENSMUSG00000034156.16', '2' ],
	[ 'TACCTTATCCCATTAT', 'ENSMUSG00000034156.16', '5' ],
	[ 'CAACTAGAGAGGTTGC', 'ENSMUSG00000034156.16', '3' ],
	[ 'CAGCGACGTTCACCTC', 'ENSMUSG00000034156.16', '4' ],
	[ 'CCAGCGAGTAGCCTCG', 'ENSMUSG00000034156.16', '8' ],
	[ 'GTCACGGAGATCCGAG', 'ENSMUSG00000034156.16', '3' ],
	[ 'CAAGCCATAATCCAAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'AGTGGGACAAACTGTC', 'ENSMUSG00000034156.16', '4' ],
	[ 'CTCGGGATCTTGTATC', 'ENSMUSG00000034156.16', '6' ],
	[ 'CGTCTACTCACCGTAA', 'ENSMUSG00000034156.16', '8' ],
	[ 'CTCTGGTCAGCTCGAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CGGCTAGTCGGACAAG', 'ENSMUSG00000034156.16', '5' ],
	[ 'CCTACACAGCGTGAAC', 'ENSMUSG00000034156.16', '2' ],
	[ 'CGAGCCATCTCGCATC', 'ENSMUSG00000034156.16', '2' ],
	[ 'ACTTTCATCAATAAGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'AACTGGTTCATCGATG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CCGGTAGTCAATCTCT', 'ENSMUSG00000034156.16', '3' ],
	[ 'TTAGGCAGTAGGGTAC', 'ENSMUSG00000034156.16', '5' ],
	[ 'ACGAGGAGTGCCTGGT', 'ENSMUSG00000034156.16', '4' ],
	[ 'GCAGCCAAGCTCAACT', 'ENSMUSG00000034156.16', '1' ],
	[ 'AGGCCGTCAGATCCAT', 'ENSMUSG00000034156.16', '3' ],
	[ 'CTCGTCAAGGGAGTAA', 'ENSMUSG00000034156.16', '10' ],
	[ 'GCAATCAGTAAGTGGC', 'ENSMUSG00000034156.16', '4' ],
	[ 'TCAATCTTCTGATACG', 'ENSMUSG00000034156.16', '5' ],
	[ 'CTCACACGTACTTAGC', 'ENSMUSG00000034156.16', '4' ],
	[ 'TACAGTGGTACGCACC', 'ENSMUSG00000034156.16', '2' ],
	[ 'TCCCGATTCGTCTGAA', 'ENSMUSG00000034156.16', '3' ],
	[ 'GTAACTGCAGTATCTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GCTGCAGGTCCATCCT', 'ENSMUSG00000034156.16', '3' ],
	[ 'GTAGGCCTCTTCGAGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'CTTATTACATTCCTCG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CGATGTACACACAGAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GAACCTAGTAGCCTAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'CTGATCCCACGGTAGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGGGAAGCAGCCTGTG', 'ENSMUSG00000034156.16', '3' ],
	[ 'CGGAGTCGTTCGGCAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CTGCTGTAGAGTAAGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CATTATCCACGTAAGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GGATGTTGTAGCGTGA', 'ENSMUSG00000034156.16', '2' ],
	[ 'TGCCCATGTCAGAAGC', 'ENSMUSG00000034156.16', '3' ],
	[ 'GCAATCACAAATTGCC', 'ENSMUSG00000034156.16', '3' ],
	[ 'GCGACCATCATAGCAC', 'ENSMUSG00000034156.16', '3' ],
	[ 'CTCGTACCAGTTCATG', 'ENSMUSG00000034156.16', '4' ],
	[ 'ATCACGAAGCCACTAT', 'ENSMUSG00000034156.16', '2' ],
	[ 'TGCTACCGTCTTGATG', 'ENSMUSG00000034156.16', '4' ],
	[ 'CAGTCCTAGACCTAGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CCTACCAAGTCCTCCA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TTTACTGGTTGACGTT', 'ENSMUSG00000034156.16', '1' ],
	[ 'GAATAAGTCTTCCTTC', 'ENSMUSG00000034156.16', '2' ],
	[ 'GGCCCATGTCAGAAGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CATCGAAAGTGTCCAT', 'ENSMUSG00000034156.16', '4' ],
	[ 'GTCACAACACGTGAGA', 'ENSMUSG00000034156.16', '2' ],
	[ 'TGCTGCTTCTGCAGTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'CCTACCAAGTCCTCCT', 'ENSMUSG00000034156.16', '3' ],
	[ 'GTCACGGTCGAATGCT', 'ENSMUSG00000034156.16', '1' ],
	[ 'TCACGAAGTTCCACGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CAAGAAAAGGCCATAG', 'ENSMUSG00000034156.16', '3' ],
	[ 'ACTTACTAGTGTCCCG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GACCTGGTCACCATAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'CACCAGGCAGCGTAAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TTATGCTGTTAGATGA', 'ENSMUSG00000034156.16', '2' ],
	[ 'TGCGCAGTCGTTGCCT', 'ENSMUSG00000034156.16', '2' ],
	[ 'CCAATCCCATCCTAGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'CAGCGACGTATGAAAC', 'ENSMUSG00000034156.16', '2' ],
	[ 'CACCAGGGTACTCGCG', 'ENSMUSG00000034156.16', '2' ],
	[ 'GACTACAGTCTCACCT', 'ENSMUSG00000034156.16', '2' ],
	[ 'CGTCAGGGTGCCTGTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TACTCGCAGTGAATTG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CATGACAGTCATGCAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'ACTTGTTTCAACGCTA', 'ENSMUSG00000034156.16', '2' ],
	[ 'CAGGTGCGTGGGTATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGGGAAGTCACCCGAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'ATCACGAGTGTGCCTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'CCGTTCATCTCCAGGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'CTGCCTAAGTTATCGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'GCATGTAGTGGTGTAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GCTGGGTCATATACGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGTGTTTTCAGTCAGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'TTGTAGGGTGCGCTTG', 'ENSMUSG00000034156.16', '3' ],
	[ 'GAGGTGCGTGGGTATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GGCTCTCGTTCAACCA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGCCCATGTCATAAGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'AAGGTTCGTCGAACAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'CTGATCCAACGGTAGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'CAAGATCTCAGCATGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'CACATAGAGCTAGCCC', 'ENSMUSG00000034156.16', '1' ],
	[ 'TCACGAATGTCTGATA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGTGGTATCCCTAATT', 'ENSMUSG00000034156.16', '1' ],
	[ 'CGCTATCTCTTGAGGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'CGACCTTTCCGCGGTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TTATGCTGTGGTCTCG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TCCCGATCAGTGAGTG', 'ENSMUSG00000034156.16', '2' ],
	[ 'GTACTCCTCTGGTATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGCGTGGGTTTGTGTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'ATCATCTCATGGTCAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'ACACTGAGTTGAACTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'AAGGTTCCACAAGACG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TAAGCGTGTTCGCGAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CGGACTGGTGCCTGGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'TCATTTGGTGTGCCTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'ATCACGAGTTGGGACA', 'ENSMUSG00000034156.16', '1' ],
	[ 'ATTACTCGTATCAGTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CATGGCGTCACGGTTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TAAGAGAAGTTAGCGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGACTAGTCAACGGCC', 'ENSMUSG00000034156.16', '1' ],
	[ 'TCTGAGATCGCCATAA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TGACGGCGTCAATGTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'GACGCGTGTCCAACTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'GCGCGATTCTTTAGGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'GGGAGATGTCAAACTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'CGGGTCACATTAACCG', 'ENSMUSG00000034156.16', '1' ],
	[ 'CACCACTAGTTATCGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'GTGCAGCAGCGAGAAA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TTCTCCTAGTGAACAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'CCAGCGAGTAAGAGGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'TACTCATAGCAACGGT', 'ENSMUSG00000034156.16', '2' ],
	[ 'AAGACCTGTGCACCAC', '1',                     '1' ],
	[ 'GCTCCTAGTCCGTGAC', '1',                     '1' ],
	[ 'CGTCTACTCACCGTAA', '1',                     '1' ],
	[ 'ACTTGTTTCAACGCTA', '1',                     '4' ],
	[ 'TACTCATAGATGTAAC', '1',                     '1' ],
	[ 'TCTGAGATCGCCATAA', '1',                     '1' ],
	[ 'GACCATTAGCTGCTTG', '1',                     '1' ],
	[ 'CCTTGTAGGCTCTGTC', '1',                     '1' ],
	[ 'TTGTAGGCTCTGTCCA', '1',                     '3' ],
	[ 'GATCAGTAGCTGCTTG', '1',                     '1' ],
	[ 'GTTAAGCAGTTAACGA', '1',                     '1' ],
	[ 'TAGTGGTAGCCTCGTG', '1',                     '1' ],
	[ 'GCAGCTTCCTCTTCAG', '1',                     '1' ],
	[ 'GCATGTAGTTCCAACA', '1',                     '1' ],
	[ 'GCTGCAGCACGGCCAT', '1',                     '1' ],
	[ 'TTGTAGGAGGTGCACA', '1',                     '4' ],
	[ 'AGAGTGGGTCCGTGAC', '1',                     '1' ],
	[ 'ACAGCCGGTAACGTTC', '1',                     '1' ],
	[ 'CCCAATCGTGACTCAT', '1',                     '1' ],
	[ 'GACCTGGCACCAACCG', '1',                     '1' ],
	[ 'TTAACTCGTTCCCATG', '1',                     '1' ],
	[ 'GCGGGCTTCACCAGGC', '1',                     '1' ],
	[ 'CGGAGTCCAGGAACGT', '1',                     '1' ],
	[ 'CCTGTTTTACTTCAAG', '1',                     '4' ],
	[ 'CACCTTGAGCTAAACA', '1',                     '1' ],
	[ 'TTAACTCTCTTCATGT', '1',                     '2' ],
	[ 'GGGACCTCACCAGGCT', '1',                     '1' ],
	[ 'GACAGAGGTACGCTGC', 'ENSMUSG00000009350.13', '1' ],
	[ 'GTCACGGAGATCCGAG', 'ENSMUSG00000009350.13', '1' ],
	[ 'CGTCTACTCACCGTAA', 'ENSMUSG00000009350.13', '1' ],
	[ 'CCTATTAGGCTCTGTC', 'ENSMUSG00000009350.13', '2' ],
	[ 'TCTGAGACATATACCG', 'ENSMUSG00000009350.13', '1' ],
	[ 'CAGCATAAGCACCGTC', 'ENSMUSG00000009350.13', '1' ],
	[ 'GACCATTAGCTGCTTG', 'ENSMUSG00000009350.13', '1' ],
	[ 'GTCACAAGCTGCTTGG', 'ENSMUSG00000009350.13', '2' ],
	[ 'CCTTGTAGGCTCTGTC', 'ENSMUSG00000009350.13', '1' ],
	[ 'GTGCAGCTCTGTCCAC', 'ENSMUSG00000009350.13', '2' ],
	[ 'TTGTAGGCTCTGTCCA', 'ENSMUSG00000009350.13', '2' ],
	[ 'GACGCTCTCTTCTGGA', 'ENSMUSG00000009350.13', '2' ],
	[ 'CCATGTCTCTGTCCGT', 'ENSMUSG00000009350.13', '3' ],
	[ 'GCTGCAGCACGGCCAT', 'ENSMUSG00000009350.13', '2' ],
	[ 'TTCCCAGAGAGTGACC', 'ENSMUSG00000009350.13', '1' ],
	[ 'TTGGAACGTGACTCAT', 'ENSMUSG00000009350.13', '1' ],
	[ 'CGAGCACAGAGGTTGC', 'ENSMUSG00000009350.13', '1' ],
	[ 'CTCTACGTCTGGTATG', 'ENSMUSG00000009350.13', '1' ],
	[ 'TACCTATAGTCAAGCG', 'ENSMUSG00000009350.13', '1' ],
	[ 'AAGCAGTGGTATCAAC', 'ENSMUSG00000009350.13', '2' ],
	[ 'TTGACTTTCCGAATGT', 'ENSMUSG00000009350.13', '1' ],
	[ 'GAGCAGGGTAGAGGAA', 'ENSMUSG00000009350.13', '1' ],
	[ 'CGTAGCGAGACCACGA', 'ENSMUSG00000009350.13', '1' ],
	[ 'CATCAGAGTTGGAGGT', 'ENSMUSG00000009350.13', '1' ],
	[ 'GATCGTACAGTCCTTC', 'ENSMUSG00000009350.13', '1' ],
	[ 'GGGACCTCACCAGGCT', 'ENSMUSG00000034121.12', '5' ]
];
is_deeply( $value->{'data'}, $exp, "all data OK" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
