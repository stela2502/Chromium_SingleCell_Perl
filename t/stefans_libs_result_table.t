#! /usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Test::More tests => 31;
BEGIN { use_ok 'stefans_libs::result_table' }

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );

my $OBJ = stefans_libs::result_table->new(
	{ 'debug' => 1, 'data_storage_spliced' => 1 } );
is_deeply( ref($OBJ), 'stefans_libs::result_table',
	'simple test of function stefans_libs::result_table -> new() ' );

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 0, 0, 0 ],
	'start'
);

ok( $OBJ->{'data_storage_spliced'}, "start setting is stored" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

$OBJ->Add_2_Header(
	[ 'Gene_ID', 'Sample1', 'Sample2', 'Sample3', 'Sample1 spliced' ] );
push( @{ $OBJ->{'data'} }, [ 'Gene1', undef, undef, 10,    undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene2', 1,     undef, 10,    1 ] );
push( @{ $OBJ->{'data'} }, [ 'Gene3', undef, 1,     undef, undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene4', undef, undef, undef, undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene6', 1,     undef, undef, undef ] );

my $outpath =
  File::Spec->catfile( $plugin_path, 'data', 'output', 'result_table' );

if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
}
else {
	system("mkdir -p $outpath");
}

my $outfile = File::Spec->catfile( $outpath, "result_table_test.db" );

$OBJ->write_table($outfile);

my $obj = stefans_libs::database::Chromium_SingleCell::datavalues->new(
	{ 'data_storage_spliced' => 1, file => $outfile } );

$value = $obj->get_data_table_4_search(
	{
		'search_columns' =>
		  [ 'datavalues.id', 'sname', 'gname', 'value', 'spliced' ],
		'where' => [],
	}
);

#warn $OBJ->info();

#warn $obj->{'complex_search'};
#print "\$exp = " . root->print_perl_var_def( [ $value->{'header'}, @{ $value->{'data'} } ] ). ";\n";

$exp = [
	[ 'datavalues.id', 'sname',   'gname', 'value', 'spliced' ],
	[ '1',             'Sample3', 'Gene1', '10',    '0' ],
	[ '2',             'Sample1', 'Gene2', '1',     '1' ],
	[ '3',             'Sample3', 'Gene2', '10',    '0' ],
	[ '4',             'Sample2', 'Gene3', '1',     '0' ],
	[ '5',             'Sample1', 'Gene6', '1',     '0' ]
];
is_deeply( [ $value->{'header'}, @{ $value->{'data'} } ],
	$exp, "data stored as expected" );

$OBJ = stefans_libs::result_table->new(
	{ 'filename' => $outfile, 'data_storage_spliced' => 1 } );

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 3, 5, 5 ],
	'I expect 5 samples 5 genes and 5 data points'
);

system("rm -Rf $outpath/*");    # restart

$OBJ = stefans_libs::result_table->new(
	{ 'filename' => $outfile, 'data_storage_spliced' => 1 } );

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 0, 0, 0 ],
	'restart'
);

$OBJ->Add_2_Header(
	[ 'Gene_ID', 'Sample1', 'Sample2', 'Sample3', 'Sample1 spliced' ] );
push( @{ $OBJ->{'data'} }, [ 'Gene1', undef, undef, 10,    undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene2', 1,     undef, 10,    1 ] );
push( @{ $OBJ->{'data'} }, [ 'Gene3', undef, 1,     undef, undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene4', undef, undef, undef, undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene6', 1,     undef, undef, undef ] );

$OBJ->write_table( undef, 1 );

ok( $OBJ->Lines() == 4, "save only one gene in the db" );

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 3, 5, 1 ],
	'line 1 added 5 genes 3 samples and one value?'
);

$OBJ->write_table( $outfile, 2 );    ## now I have added 3 lines and 4 values

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 3, 5, 4 ],
	'line 2 and 3 added 5 genes 3 samples and 4 values?'
);

$OBJ->write_table($outfile);

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 3, 5, 5 ],
	'all lines added 5 genes 3 samples and 5 values?'
);

push( @{ $OBJ->{'data'} }, [ 'Gene7', 1, 5, 6, 1 ] );

#warn "Now the Gene7 should be added?\n". join(" ", @{$OBJ->GetAsArray('Gene_ID') } )."\n";

$OBJ->write_table($outfile);

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 3, 6, 8 ],
	'one new line added 6 genes 3 samples and 8 values?'
);

#print "\$exp = " . root->print_perl_var_def( $OBJ->{__gene2id__} ). ";\n";

$OBJ->Add_2_Header('Sample4');
push( @{ $OBJ->{'data'} }, [ 'Gene8', undef, undef, undef, undef, 19 ] );
$OBJ->write_table($outfile);
is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 4, 7, 9 ],
	'one new line added 7 genes 4 samples and 9 values?'
);

## now restart and check the AddDataset usability by asking the tool to reload the database:

#system("rm -Rf $outpath/*"); # restart

$OBJ = stefans_libs::result_table->new( { 'data_storage_spliced' => 1 } );

## lest add one value...

$OBJ->AddDataset(
	{ 'Gene_ID' => 'Gene13', 'Sample_ID' => 'Sample45', 'value' => 12 }, 1 );

is_deeply( $OBJ->{'data'}, [ [ 'Gene13', '12' ] ], "Data added correctly" );

$OBJ->AddDataset(
	{ 'Gene_ID' => 'Gene13', 'Sample_ID' => 'Sample45', 'value' => 12 }, 1 );

is_deeply( $OBJ->{'data'}, [ [ 'Gene13', '24' ] ], "Data added (+) correctly" );

$OBJ->AddDataset(
	{ 'Gene_ID' => 'Gene13', 'Sample_ID' => 'Sample45', 'value' => 12 }, 0 );

is_deeply( $OBJ->{'data'}, [ [ 'Gene13', '12' ] ], "Data replaced correctly" );

$OBJ->AddDataset(
	{ 'Gene_ID' => 'Gene13', 'Sample_ID' => 'Sample5', 'value' => 12 }, 0 );

is_deeply(
	$OBJ->{'data'},
	[ [ 'Gene13', '12', '12' ] ],
	"Sample added correctly"
);

$OBJ->AddDataset(
	{ 'Gene_ID' => 'Gene3', 'Sample_ID' => 'Sample5', 'value' => 12 }, 0 );

is_deeply(
	$OBJ->{'data'},
	[ [ 'Gene13', '12', '12' ], [ 'Gene3', undef, 12 ] ],
	"Gene added correctly"
);

$OBJ = stefans_libs::result_table->new( { 'data_storage_spliced' => 1 } );

$OBJ = $OBJ->import_database($outfile);

is_deeply( ref($OBJ), 'stefans_libs::result_table', 'the right object' );

is_deeply(
	[ $OBJ->Columns(), $OBJ->Rows() ],
	[ 6,               6 ],
	"got 6 rows and 6 columns in the data table"
);

## wher IS THE FUCKING PROBLEM ???!! Why is nothing stored in the object?
is_deeply(
	$OBJ->{'header'},
	[
		'Gene_ID', 'Sample1', 'Sample2', 'Sample3', 'Sample4',
		'Sample1 spliced'
	],
	'Correct header in the samples table'
);

#print $OBJ -> AsString();
$exp = [
	[ 'Gene1', undef, undef, '10' ],
	[ 'Gene2', '1',   undef, '10', undef, '1' ],
	[ 'Gene3', undef, '1' ],
	[ 'Gene6', '1' ],
	[ 'Gene7', '1',   '5',   '6',  undef, '1' ],
	[ 'Gene8', undef, undef, undef, '19' ]
];

#print "\$exp = " . root->print_perl_var_def(  $OBJ->{'data'}   ). ";\n";
is_deeply( $OBJ->{'data'}, $exp, 'Correct data in the samples table' );

is_deeply(
	[
		$OBJ->{'__lastSampleID__'}, $OBJ->{'__lastGeneID__'},
		$OBJ->{'__lastDataID__'}
	],
	[ 4, 7, 9 ],
	'database imported correctly'
);

## now test the result tables way - should be way faster for saving initailly
$OBJ = stefans_libs::result_table->new( { 'data_storage_spliced' => 1 } );
## quick and simple create object
$OBJ->Add_2_Header(
	[ 'Gene_ID', 'Sample1', 'Sample2', 'Sample3', 'Sample1 spliced' ] );
push( @{ $OBJ->{'data'} }, [ 'Gene1', undef, undef, 10,    undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene2', 1,     undef, 10,    1 ] );
push( @{ $OBJ->{'data'} }, [ 'Gene3', undef, 1,     undef, undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene4', undef, undef, undef, undef ] );
push( @{ $OBJ->{'data'} }, [ 'Gene6', 1,     undef, undef, undef ] );

my $ofile = File::Spec->catfile( $outpath, "table_test.db" );

$OBJ->print2table($ofile);

$ofile = File::Spec->catfile( $outpath, "table_test" );
ok( -d $ofile, "outpath for table export created" );

foreach my $name ( 'barcodes.tsv', 'genes.tsv', 'matrix.mtx' ) {
	$exp = File::Spec->catfile( $ofile, $name );
	ok( -f $exp, "outfile $name created ($exp)" );
}

$exp = File::Spec->catfile( $ofile, 'genes.tsv' );
if ( -f $exp ) {
	open( IN, "<$exp" );
	$value = [ map { chomp; [ split( " ", $_ ) ] } <IN> ];

	#print "\$value = " . root->print_perl_var_def( $value ). ";\n";

}
else {
	$value = [];
}
$exp = [ ['Gene1'], ['Gene2'], ['Gene3'], ['Gene4'], ['Gene6'] ];

is_deeply( $value, $exp, "file 'genes.tsv'" );

$exp = File::Spec->catfile( $ofile, 'barcodes.tsv' );
if ( -f $exp ) {
	open( IN, "<$exp" );
	$value = [ map { chomp; [ split( " ", $_ ) ] } <IN> ];

	#print "\$value = " . root->print_perl_var_def( $value ). ";\n";
}
else {
	$value = [];
}
$exp = [ ['Sample1'], ['Sample2'], ['Sample3'] ];

is_deeply( $value, $exp, "file 'barcodes.tsv'" );

$exp = File::Spec->catfile( $ofile, 'matrix.mtx' );
if ( -f $exp ) {
	open( IN, "<$exp" );
	$value = [ map { chomp; [ split( " ", $_ ) ] } <IN> ];

	#print "\$value = " . root->print_perl_var_def( $value ). ";\n";
}
else {
	$value = [];
}
$exp = [
	[ '1', '3', '10', '0' ],    #Gene1, Sample3, 10, 0
	[ '2', '1', '1',  '1' ],    #Gene2, Sample1, 1,  1
	[ '2', '3', '10', '0' ],    #Gene2, Sample3, 10, 0
	[ '3', '2', '1',  '0' ],    #Gene3, Sample2, 1,  0
	[ '5', '1', '1',  '0' ]     #Gene6, Smaple1, 1,  0
];

is_deeply( $value, $exp, "file 'matrix.mtx'" );

$exp = stefans_libs::result_table->new( { 'data_storage_spliced' => 1 } );
## quick and simple create object
$exp->Add_2_Header(
	[ 'Gene_ID', 'Sample1', 'Sample2', 'Sample3', 'Sample1 spliced' ] );
push( @{ $exp->{'data'} }, [ 'Gene1', undef, undef, 10,    undef ] );
push( @{ $exp->{'data'} }, [ 'Gene2', 1,     undef, 10,    1 ] );
push( @{ $exp->{'data'} }, [ 'Gene3', undef, 1,     undef, undef ] );
push( @{ $exp->{'data'} }, [ 'Gene4', undef, undef, undef, undef ] );
push( @{ $exp->{'data'} }, [ 'Gene6', 1,     undef, undef, undef ] );

$OBJ = stefans_libs::result_table->new();

$OBJ->import_tables($ofile);

for ( my $i = 0 ; $i < 4 ; $i++ ) {
	my $array = @{ $OBJ->{'data'} }[$i];
	for ( my $a = scalar(@$array) - 1 ; $a < 4 ; $a++ ) {
		push( @$array, undef );
	}
}
splice( @{ $exp->{'data'} }, 3, 1 );    ## get rid of the empty line of Gene4

is_deeply(
	[ $OBJ->{'header'}, $OBJ->{'data'} ],
	[ $exp->{'header'}, $exp->{'data'} ],
	"import_tables one path"
);

## now for the finish try to read the same data from 2 different folders
my @tmp = (
	File::Spec->catfile( $outpath, 'chr1.db' ),
	File::Spec->catfile( $outpath, 'chr2.db' )
);

$exp = stefans_libs::result_table->new( { 'data_storage_spliced' => 1 } );
## quick and simple create object
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene1',
		'Sample_ID' => 'Sample3',
		'value'     => 10
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene2',
		'Sample_ID' => 'Sample1',
		'value'     => 1
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene2',
		'Sample_ID' => 'Sample1 spliced',
		'value'     => 1
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene2',
		'Sample_ID' => 'Sample3',
		'value'     => 10
	},
	0
);
$exp->print2table( $tmp[0], );

$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene3',
		'Sample_ID' => 'Sample2',
		'value'     => 1
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene6',
		'Sample_ID' => 'Sample1',
		'value'     => 1
	},
	0
);

$exp->print2table( $tmp[1] );

## now I need to regenerate my exp table....
$exp = stefans_libs::result_table->new( { 'data_storage_spliced' => 0 } );

$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene1',
		'Sample_ID' => 'Sample3',
		'value'     => 10
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene2',
		'Sample_ID' => 'Sample1',
		'value'     => 1
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene2',
		'Sample_ID' => 'Sample1 spliced',
		'value'     => 1
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene2',
		'Sample_ID' => 'Sample3',
		'value'     => 10
	},
	0
);

$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene3',
		'Sample_ID' => 'Sample2',
		'value'     => 1
	},
	0
);
$exp->AddDataset(
	{
		'Gene_ID'   => 'Gene6',
		'Sample_ID' => 'Sample1',
		'value'     => 1
	},
	0
);

$OBJ = stefans_libs::result_table->new();
@tmp = map { $_ =~ s/.db$//; $_ } @tmp;
$OBJ->import_tables(@tmp);

#$exp = [ [ 'Gene_ID', 'Sample3', 'Sample1', 'Sample1 spliced', 'Sample2' ], [ 'Gene1', '10' ], [ 'Gene2', '10', '1', '1' ] ];

#is_deeply(
#	[ $OBJ->{'header'}, $OBJ->{'data'} ],
#	[ $exp->{'header'}, $exp->{'data'} ],
#	"import_tables from two different folders"
#);

#print "\$exp = " . root->print_perl_var_def( [ $OBJ->{'header'}, @{ $OBJ->{'data'} } ] ). ";\n";
