#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 3;
BEGIN { use_ok 'stefans_libs::GeneModelMatcher' }

use File::Spec::Functions;

use stefans_libs::root;
use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );

my $filename = File::Spec->catfile($plugin_path, 'data','GOI.annotation.genecode.v12.gff3');

my $OBJ = stefans_libs::GeneModelMatcher -> new({'debug' => 1, 'filename' => $filename});
is_deeply ( ref($OBJ) , 'stefans_libs::GeneModelMatcher', 'simple test of function stefans_libs::gtf_file::GeneModelMatcher -> new() ');

#print "\$exp = ".root->print_perl_var_def( $OBJ->{'data'}  ).";\n";
is_deeply( 6, scalar( @{$OBJ->{'data'}}), "6 gene models? (".scalar( @{$OBJ->{'data'}}).")" );

#print "\$exp = ".root->print_perl_var_def( $OBJ->{'gtf_file'}->{'data'}  ).";\n";
$exp = [ 
	[ 'chr19', 'HAVANA', 'gene', '5795690', '5802672', '.', '-', '.', 'ENSMUSG00000092341.2' ], 
	[ 'chr3', 'HAVANA', 'gene', '90668978', '90670035', '.', '+', '.', 'ENSMUSG00000056054.9' ], 
	[ 'chr3', 'HAVANA', 'gene', '90692632', '90695721', '.', '-', '.', 'ENSMUSG00000056071.12' ], 
	[ 'chr3', 'HAVANA', 'gene', '90726906', '90741517', '.', '+', '.', 'ENSMUSG00000042250.13' ], 
	[ 'chrM', 'ENSEMBL', 'gene', '7927', '8607', '.', '+', '.', 'ENSMUSG00000064357.1' ], 
	[ 'chrX', 'HAVANA', 'gene', '167207093', '167209315', '.', '-', '.', 'ENSMUSG00000049775.16' ] 
];

is_deeply( $OBJ->{'gtf_file'}->{'data'}, $exp, "correct gene information" );


$OBJ -> {'gtf_file'}->{'slice_length'} = 1e+5;
$OBJ -> {'gtf_file'}->GeneFreeSplits(); ## define the gene free splits in order to keep problems with missing genes as low as possible.


@values = $OBJ->genes_at_position_plus_one('chr3', '90669000',  '90669100' );

is_deeply( scalar(@values), 2 , "got two results models back" );

$value = $OBJ->match_cigar( 'chr3', 90668978, "4S120M209N10M" );

$exp = {'ENSMUSG00000056054.9' => 'spliced'};

is_deeply( $value, $exp , "match_cigar with spliced cigar info" );


#90726952        90727000
$value = $OBJ->match_cigar( 'chr3', 90726952, "48M" );

$exp = {'ENSMUSG00000042250.13' => 'exon'};

#is_deeply( $value, $exp , "match_cigar with spliced cigar info" );

print "\$exp = ".root->print_perl_var_def( $value ).";\n";
