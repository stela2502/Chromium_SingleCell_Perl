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
my $outpath =
  "$plugin_path/data/output/QuantifyBamFile_chr11_87760834_87856494";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$gtf_file = "$plugin_path/data/mm10.Mpo_large.gtf";

#$gtf_file = "/home/stefanl/lunarc/genomes/mouse/mm10/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf";

ok( -f $gtf_file, 'test gtf file exists' );

$infile =
  "$plugin_path/data/mm10.singleCells.merged.chr11_87760834_87856494.bam";
ok( -f $infile, 'test infile exists' );

$outfile = "$outpath/singleCells.merged.chr11_87760834_87856494_quantifed";

#improvement by using UMIs to merge samples:
#no usage: 61773 samples and 4 genes
#with umi: 25261 and 4 genes
#with umi and sample merge: 16224 samples and 4 genes detected

# with only 4 geens we are down to only half the samples.

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
"$outpath/singleCells.merged.chr11_87760834_87856494_quantifed.original_merged.db",
		'data_storage_spliced' => 1
	}
);

$value = $db->get_data_table_4_search(
	{
		'search_columns' => [ 'sname', 'gname', 'value' ],
		'where'          => [],
	}
);

ok( $value->Rows() == 158, "158 samples(" . $value->Rows() . ")" );

#print "\$exp = ".root->print_perl_var_def([keys %{$value->createIndex('gname')}]).";\n";

$exp = [ '1', 'ENSMUSG00000009350.13', 'ENSMUSG00000034156.16' ];

is_deeply( [ sort keys %{ $value->createIndex('gname') } ],
	$exp, "only one gene ENSMUSG00000022103.9" );

print "\$exp = " . root->print_perl_var_def( $value->{'data'} ) . ";\n";

$exp = [
	[ 'C_CCAATCCAGCTACCGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACCCACTGTTACCAGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGGCTGGCACCTCGGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGCGGTACAGCCTTGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_AGTGTCATCACCGGGT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TTGTAGGGTTTGCATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GGCAATTGTACCGCTG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_ATTATCCAGCCCAACC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GACAGAGGTACGCTGC', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CACTCCAGTTGGTTTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_AAGACCTGTGCACCAC', 'ENSMUSG00000034156.16', '3' ],
	[ 'C_GCAGTTAAGACGACGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCGACCAAGTCCTCCT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_AACTGGTTCCACTGGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GGAATAATCACGGTTA', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TTTGGTTCACACGCTG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_ATGGGAGGTCTAGTGT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TTCGAAGAGTTGTCGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACGGCCAGTTCCTCCA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGCTACCAGTGTCCAT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CATCCACGTAAGTGGC', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_AGGTCCGGTTCACCTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TCAGGTATCAGCACAT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CGAACATCATCCTAGA', 'ENSMUSG00000034156.16', '3' ],
	[ 'C_AAAGTAGGTACCAGTT', 'ENSMUSG00000034156.16', '7' ],
	[ 'C_TCACGAATCACGCATA', 'ENSMUSG00000034156.16', '4' ],
	[ 'C_AAGGTTCTCCTCAATT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CTCCTAGTCTACGAGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACACCGGGTTGTTTGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCTCCTAGTCCGTGAC', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CTTGGCTGTTTCCACC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAGCATAAGTTTCCTT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GATCGTACATCGATTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACATGGTAGATCTGCT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CATCGAACATGCTGGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCGGCTGTTTGGCGCG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TACCTTATCCCATTAT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CCAGCGAGTAGCCTCG', 'ENSMUSG00000034156.16', '3' ],
	[ 'C_CAGCGACGTTCACCTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAACTAGAGAGGTTGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GTCACGGAGATCCGAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAAGCCATAATCCAAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_AGTGGGACAAACTGTC', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CTCGGGATCTTGTATC', 'ENSMUSG00000034156.16', '3' ],
	[ 'C_CGTCTACTCACCGTAA', 'ENSMUSG00000034156.16', '5' ],
	[ 'C_CTCTGGTCAGCTCGAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGGCTAGTCGGACAAG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CCTACACAGCGTGAAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGAGCCATCTCGCATC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACTTTCATCAATAAGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TTAGGCAGTAGGGTAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_AACTGGTTCATCGATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CCGGTAGTCAATCTCT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACGAGGAGTGCCTGGT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_GCAATCAGTAAGTGGC', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TCAATCTTCTGATACG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CTCGTCAAGGGAGTAA', 'ENSMUSG00000034156.16', '6' ],
	[ 'C_AGGCCGTCAGATCCAT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TCCCGATTCGTCTGAA', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TACAGTGGTACGCACC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GTAACTGCAGTATCTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCTGCAGGTCCATCCT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CTCACACGTACTTAGC', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CGATGTACACACAGAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CTGCTGTAGAGTAAGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GGATGTTGTAGCGTGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCAATCACAAATTGCC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGGAGTCGTTCGGCAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ATCACGAAGCCACTAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAGTCCTAGACCTAGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TTTACTGGTTGACGTT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CCTACCAAGTCCTCCA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GGCCCATGTCAGAAGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CATCGAAAGTGTCCAT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_CTCGTACCAGTTCATG', 'ENSMUSG00000034156.16', '3' ],
	[ 'C_GTCACGGTCGAATGCT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GTCACAACACGTGAGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TCACGAAGTTCCACGG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_ACTTACTAGTGTCCCG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CACCAGGCAGCGTAAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GACCTGGTCACCATAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TTATGCTGTTAGATGA', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_TGCGCAGTCGTTGCCT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CCAATCCCATCCTAGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGCCCATGTCAGAAGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CACCAGGGTACTCGCG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_GACTACAGTCTCACCT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGTCAGGGTGCCTGTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TACTCGCAGTGAATTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CATGACAGTCATGCAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAGGTGCGTGGGTATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGGGAAGTCACCCGAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CCGTTCATCTCCAGGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGGGAAGCAGCCTGTG', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_GCATGTAGTGGTGTAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAAGAAAAGGCCATAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCTGGGTCATATACGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGTGTTTTCAGTCAGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TTGTAGGGTGCGCTTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCGACCATCATAGCAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAGCGACGTATGAAAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGCTACCGTCTTGATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGCCCATGTCATAAGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_AAGGTTCGTCGAACAG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GAATAAGTCTTCCTTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CAAGATCTCAGCATGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CACATAGAGCTAGCCC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACTTGTTTCAACGCTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGTGGTATCCCTAATT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGCTATCTCTTGAGGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGACCTTTCCGCGGTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TCCCGATCAGTGAGTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GTACTCCTCTGGTATG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ATCATCTCATGGTCAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ACACTGAGTTGAACTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TAAGCGTGTTCGCGAC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGGACTGGTGCCTGGT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ATCACGAGTTGGGACA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_ATTACTCGTATCAGTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TGACTAGTCAACGGCC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TAAGAGAAGTTAGCGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CATGGCGTCACGGTTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TCATTTGGTGTGCCTG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TCTGAGATCGCCATAA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GACGCGTGTCCAACTA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GCGCGATTCTTTAGGG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_GGGAGATGTCAAACTC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CGGGTCACATTAACCG', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CACCACTAGTTATCGC', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TTCTCCTAGTGAACAT', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_CCAGCGAGTAAGAGGA', 'ENSMUSG00000034156.16', '1' ],
	[ 'C_TACTCATAGCAACGGT', 'ENSMUSG00000034156.16', '2' ],
	[ 'C_AAGACCTGTGCACCAC', '1',                     '1' ],
	[ 'C_CGTCTACTCACCGTAA', '1',                     '1' ],
	[ 'C_GCTGCTTGTATCACCA', '1',                     '1' ],
	[ 'C_GACCATTAGCTGCTTG', '1',                     '1' ],
	[ 'C_AGAGTGGGTCCGTGAC', '1',                     '1' ],
	[ 'C_ACAGCCGGTAACGTTC', '1',                     '1' ],
	[ 'C_CCCAATCGTGACTCAT', '1',                     '1' ],
	[ 'C_GACCTGGCACCAACCG', '1',                     '1' ],
	[ 'C_GATCGTACAGTCCTTC', '1',                     '1' ],
	[ 'C_TTAACTCGTTCCCATG', '1',                     '1' ],
	[ 'C_GCGGGCTTCACCAGGC', '1',                     '1' ],
	[ 'C_GGGACCTCACCAGGCT', '1',                     '1' ],
	[ 'C_CTCGGGATCTTGTATC', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_CCTATTAGGCTCTGTC', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_GTCACAAGCTGCTTGG', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_TTGTAGGCTCTGTCCA', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_GTGCAGCTCTGTCCAC', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_ATTGGTGTCTGGTATG', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_CTCTACGTCTGGTATG', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_TTGGAACGTGACTCAT', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_CGAGCACAGAGGTTGC', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_TACCTATAGTCAAGCG', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_AAGCAGTGGTATCAAC', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_TTGACTTTCCGAATGT', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_GAGCAGGGTAGAGGAA', 'ENSMUSG00000009350.13', '1' ],
	[ 'C_CATCAGAGTTGGAGGT', 'ENSMUSG00000009350.13', '1' ]
];
is_deeply( $value->{'data'}, $exp, "all data OK" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
