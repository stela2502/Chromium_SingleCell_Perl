#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 18;
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

print "\$exp = ".root->print_perl_var_def( $OBJ->{'gtf_file'}->{'data'}  ).";\n";
$exp = [ 
	[ 'chr19', 'HAVANA', 'gene', '5795690', '5802672', '.', '-', '.', 'ENSMUSG00000092341.2' ], 
	[ 'chr3', 'HAVANA', 'gene', '90668978', '90670035', '.', '+', '.', 'ENSMUSG00000056054.9' ], 
	[ 'chr3', 'HAVANA', 'gene', '90692632', '90695721', '.', '-', '.', 'ENSMUSG00000056071.12' ], 
	[ 'chr3', 'HAVANA', 'gene', '90726906', '90741517', '.', '+', '.', 'ENSMUSG00000042250.13' ], 
	[ 'chrM', 'ENSEMBL', 'gene', '7927', '8607', '.', '+', '.', 'ENSMUSG00000064357.1' ], 
	[ 'chrX', 'HAVANA', 'gene', '167207093', '167209315', '.', '-', '.', 'ENSMUSG00000049775.16' ] 
];

is_deeply( $OBJ->{'gtf_file'}->{'data'}, $exp, "correct gene information" );


@values = $OBJ->genes_at_position_plus_one('chr3', '90669000',  '90669100' );
is_deeply( scalar(@values), 2 , "got two results models back" );

$exp = {
  'ID' => 'ENSMUSG00000056054.9',
  'chr' => 'chr3',
  'end' => '90670035',
  'gene_id' => 'ENSMUSG00000056054.9',
  'gene_name' => 'S100a8',
  'gene_type' => 'protein_coding',
  'havana_gene' => 'OTTMUSG00000052758.1',
  'level' => '2',
  'orientation' => '+',
  'start' => '90668978'
};
#print "\$exp = ".root->print_perl_var_def(  $values[0]->AsHash() ).";\n";
is_deeply( $values[0]->AsHash() , $exp , "right internals 1" );

$values[0]->match('90668978', '90670035');
 
#print "\$exp = ".root->print_perl_var_def(  $values[0]->{'positions'} ).";\n";
$exp = [ '90668978', '90669102', '90669311', '90669682', '90669839', '90670035' ];
is_deeply( $values[0]->{'positions'}, $exp , "right state changes" );

#print "\$exp = ".root->print_perl_var_def(  $values[1]->AsHash() ).";\n";
$exp = {
  'ID' => 'ENSMUSG00000056071.12',
  'chr' => 'chr3',
  'end' => '90695721',
  'gene_id' => 'ENSMUSG00000056071.12',
  'gene_name' => 'S100a9',
  'gene_type' => 'protein_coding',
  'havana_gene' => 'OTTMUSG00000026669.1',
  'level' => '2',
  'orientation' => '-',
  'start' => '90692632'
};
is_deeply( $values[1]->AsHash() , $exp , "right internals 2" );

## and now lets try to use the match_cigar

$value = $OBJ->match_cigar( 'chr3', 90668978, "4S120M209N10M" );

$exp = {'ENSMUSG00000056054.9' => 'spliced'};

is_deeply( $value, $exp , "match_cigar with spliced cigar info" );

$value = $OBJ->match_cigar( 'chr3', 90668978, "4S120M" );

$exp = {'ENSMUSG00000056054.9' => 'exon'};

is_deeply( $value, $exp , "match_cigar with exon cigar info" );

$value = $OBJ->match_cigar( 'chr3', 90669102, "4S120M" );

$exp = {'ENSMUSG00000056054.9' => 'primary'};

is_deeply( $value, $exp , "match_cigar with intron cigar info" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

$OBJ = stefans_libs::GeneModelMatcher -> new({'debug' => 1});

$OBJ->read_file( $filename, 'chr3' );
#print "\$exp = ".root->print_perl_var_def( $OBJ->{'data'}  ).";\n";
is_deeply( 3, scalar( @{$OBJ->{'data'}}), "3 gene models? (".scalar( @{$OBJ->{'data'}}).")" );

#print "\$exp = ".root->print_perl_var_def( $OBJ->{'gtf_file'}->{'data'}  ).";\n";
$exp = [ 
	[ 'chr3', 'HAVANA', 'gene', '90668978', '90670035', '.', '+', '.', 'ENSMUSG00000056054.9' ], 
	[ 'chr3', 'HAVANA', 'gene', '90692632', '90695721', '.', '-', '.', 'ENSMUSG00000056071.12' ], 
	[ 'chr3', 'HAVANA', 'gene', '90726906', '90741517', '.', '+', '.', 'ENSMUSG00000042250.13' ], 
#	[ 'chr19', 'HAVANA', 'gene', '5795690', '5802672', '.', '-', '.', 'ENSMUSG00000092341.2' ], 
#	[ 'chrX', 'HAVANA', 'gene', '167207093', '167209315', '.', '-', '.', 'ENSMUSG00000049775.16' ], 
#	[ 'chrM', 'ENSEMBL', 'gene', '7927', '8607', '.', '+', '.', 'ENSMUSG00000064357.1' ] 
];
is_deeply( $OBJ->{'gtf_file'}->{'data'}, $exp, "correct gene information chr3 only" );

$value = $OBJ->match_cigar( 'chr3', 90668978, "4S120M209N10M" );

$exp = {'ENSMUSG00000056054.9' => 'spliced'};

is_deeply( $value, $exp , "match_cigar with spliced cigar info" );


## Lets try another gtf file:
$filename = File::Spec->catfile($plugin_path, 'data','Spliced_Reads.gtf');

$OBJ = stefans_libs::GeneModelMatcher -> new({'debug' => 1, 'filename' => $filename});
is_deeply ( ref($OBJ) , 'stefans_libs::GeneModelMatcher', 'simple test of function stefans_libs::gtf_file::GeneModelMatcher -> new() ');

is_deeply( 34, scalar( @{$OBJ->{'data'}}), "34 gene models? (".scalar( @{$OBJ->{'data'}}).")" );

print "\$exp = ".root->print_perl_var_def( $OBJ->{'gtf_file'}->{'data'} ).";\n";
$exp = [ 
[ 'JH792828.1', 'HAVANA', 'gene', '46372', '46523', '.', '+', '.', 'ENSMUSG00000107717.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '86001', '87834', '.', '+', '.', 'ENSMUSG00000097606.2' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '97437', '98987', '.', '+', '.', 'ENSMUSG00000107958.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '102327', '103844', '.', '-', '.', 'ENSMUSG00000107963.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '109020', '111963', '.', '-', '.', 'ENSMUSG00000108159.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '115996', '116478', '.', '-', '.', 'ENSMUSG00000107974.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '120113', '120272', '.', '-', '.', 'ENSMUSG00000108182.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '128146', '144832', '.', '+', '.', 'ENSMUSG00000107547.1' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '136802', '137105', '.', '-', '.', 'ENSMUSG00000111516.1' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '137736', '137895', '.', '-', '.', 'ENSMUSG00000111076.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '154498', '169029', '.', '-', '.', 'ENSMUSG00000107871.1' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '217300', '217382', '.', '+', '.', 'ENSMUSG00000097826.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '217312', '217369', '.', '-', '.', 'ENSMUSG00000106761.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '217458', '248370', '.', '-', '.', 'ENSMUSG00000097297.3' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '217593', '217674', '.', '+', '.', 'ENSMUSG00000096925.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '217863', '217944', '.', '+', '.', 'ENSMUSG00000097782.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '217872', '217933', '.', '-', '.', 'ENSMUSG00000106921.1' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '218156', '218234', '.', '+', '.', 'ENSMUSG00000097568.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '219017', '219096', '.', '+', '.', 'ENSMUSG00000097688.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '219315', '219398', '.', '+', '.', 'ENSMUSG00000097034.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '219447', '219515', '.', '+', '.', 'ENSMUSG00000097400.2' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '263777', '264620', '.', '-', '.', 'ENSMUSG00000107910.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '283159', '284038', '.', '-', '.', 'ENSMUSG00000107437.1' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '289100', '298942', '.', '+', '.', 'ENSMUSG00000097829.2' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '302048', '329611', '.', '+', '.', 'ENSMUSG00000097449.8' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '331471', '366737', '.', '+', '.', 'ENSMUSG00000096998.3' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '336061', '336395', '.', '+', '.', 'ENSMUSG00000108051.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '376637', '388596', '.', '-', '.', 'ENSMUSG00000107636.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '389199', '414164', '.', '+', '.', 'ENSMUSG00000106728.3' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '390112', '399442', '.', '-', '.', 'ENSMUSG00000107912.1' ], 
[ 'JH792828.1', 'ENSEMBL', 'gene', '413582', '413672', '.', '+', '.', 'ENSMUSG00000097499.2' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '420181', '421961', '.', '-', '.', 'ENSMUSG00000108281.1' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '422711', '434183', '.', '+', '.', 'ENSMUSG00000097894.8' ], 
[ 'JH792828.1', 'HAVANA', 'gene', '449873', '473536', '.', '+', '.', 'ENSMUSG00000107537.1' ] 
];


is_deeply( 34, scalar( @{$OBJ->{'data'}}), "34 gene models? (".scalar( @{$OBJ->{'data'}}).")" );

$OBJ -> {'gtf_file'}->{'slice_length'} = 1e+5;
$OBJ -> {'gtf_file'}->GeneFreeSplits(); ## define the gene free splits in order to keep problems with missing genes as low as possible.


@values = $OBJ->genes_at_position_plus_one('JH792828.1', '217593',  '217674' );

is_deeply( scalar(@values), 3 , "got two results models back" );




