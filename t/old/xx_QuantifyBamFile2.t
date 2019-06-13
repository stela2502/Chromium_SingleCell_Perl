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
my $outpath = "$plugin_path/data/output/QuantifyBamFile_small";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$gtf_file = "$plugin_path/data/mm10.Gfra2_large.gtf";

#$gtf_file = "$plugin_path/data/mm10.chr14.gtf";

#$gtf_file = "/home/stefanl/lunarc/genomes/mouse/mm10/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf";

ok( -f $gtf_file, 'test gtf file exists' );

$infile = "$plugin_path/data/mm10.singleCells.merged.Gfra2.bam";
ok( -f $infile, 'test infile exists' );

$outfile = "$outpath/singleCells.merged.Gfra2_quantifed";

#improvement by using UMIs to merge samples:
#no usage: 61773 samples and 4 genes
#with umi: 25261 and 4 genes
#with umi and sample merge: 16224 samples and 4 genes detected

# with only 4 geens we are down to only half the samples.

@options = (
	'min_UMIs', 1

	  #,'quantify_on', 'transcript'
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
		  "$outpath/singleCells.merged.Gfra2_quantifed.original_merged.db",
		'data_storage_spliced' => 1
	}
);

$value = $db->get_data_table_4_search(
	{
		'search_columns' => [ 'sname', 'gname', 'value' ],
		'where'          => [],
	}
);

ok( $value->Rows() == 48, "48 samples (".$value->Rows().")" );

#print "\$exp = ".root->print_perl_var_def([keys %{$value->createIndex('gname')}]).";\n";

is_deeply( [ keys %{ $value->createIndex('gname') } ],
	["ENSMUSG00000022103.9"], "only one gene ENSMUSG00000022103.9" );

#print "\$exp = ".root->print_perl_var_def($value->{'data'}).";\n";

$exp = [
	[ 'C_TTAGGACAGTCTCGGC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CCTAGCTAGGAATCGC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CGGACGTCAGTATAAG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_ACCCACTGTTACCAGT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TTTGTCAGTCTCCATC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_AAAGTAGTCAGCTCTC', 'ENSMUSG00000022103.9', '2' ],
	[ 'C_GTCAAGTAGACCCACC', 'ENSMUSG00000022103.9', '4' ],
	[ 'C_TCCCACTGCCACGACA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TAAGCGTGTTCGCGAC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GACTACAGTACTTCTT', 'ENSMUSG00000022103.9', '2' ],
	[ 'C_TAAGCGTCATTCTCAT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TACTCATTCCTAGAAC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GTACTCCTCCTCCTAG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_AGCAGCCAGATGCGAC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CGTAGCGGTCTGCAAT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_AGACGTTTCGGAAATA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CATGGCGGTACAGCAG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CCTATTATCGGTTCGG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TGGCTGGGTTGGGACA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CATCGAAGTTGGGACA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CCGTGGAGTAGTACCT', 'ENSMUSG00000022103.9', '3' ],
	[ 'C_CGTTCTGGTAGAAGGA', 'ENSMUSG00000022103.9', '4' ],
	[ 'C_AGAGCGAAGGCAGTCA', 'ENSMUSG00000022103.9', '2' ],
	[ 'C_TACTCATCAAGCTGGA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GTGTGCGTCACCTCGT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GTTCTCGTCTCACATT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GGGAGATTCTGTCTCG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TGCCCTACACCACCAG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GTGCTTCGTTGTCGCG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_ATTATCCCATCACAAC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_AACTCAGTCGGAAACG', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CCTTCCCTCCCATTTA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TACGGGCGTCTTTCAT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GTCAAGTAACCCCACC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CAGCTAAAGGCCCTCA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_ATCCACCGTAAACACA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TGGCTGGCACCTCGGA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GTTAAGCGTAAATGAC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GGTATTGTCTTGTACT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TGATTTCTCTACTATC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_ATGAGGGCATCCCACT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CGGACACCATCAGTCA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TTAGGCACACGCCAGT', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_CTCAGAAGTGCACTTA', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GCTCTGTGTTCATGGT', 'ENSMUSG00000022103.9', '2' ],
	[ 'C_GCTCTGTGTCAATACC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_TTGGAACCATCGGGTC', 'ENSMUSG00000022103.9', '1' ],
	[ 'C_GGAATAAGTGATGTGG', 'ENSMUSG00000022103.9', '1' ]
];


is_deeply( $value->{'data'}, $exp, "all data OK" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
