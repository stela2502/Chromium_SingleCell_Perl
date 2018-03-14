#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 15;
use stefans_libs::flexible_data_structures::data_table;
use stefans_libs::database::Chromium_SingleCell::datavalues;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outfile, @paths, );

my $exec = $plugin_path . "/../bin/JoinResultTables.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/JoinResultTables";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

@paths = map { "$plugin_path/data/JoinResultTables/$_" } 'A', 'B';

$outfile = $outpath . "/test.sqlite";

foreach (@paths) {
	ok( -d $_,                "path $_ exists" );
	ok( -f "$_/matrix.mtx",   "matrix file exists" );
	ok( -f "$_/genes.tsv",    "genes file exists" );
	ok( -f "$_/barcodes.tsv", "barcodes file exists" );
}

ok( !-f $outfile, "outfile not existing" );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . " -outfile "
  . $outfile
  . " -paths "
  . join( ' ', @paths )
  . " -createTable "    # or not?
  . " -force"

  #. " -debug"
  ;
my $start = time;
system($cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

ok( -f $outfile, "outfile created" );

my $datavalues = stefans_libs::database::Chromium_SingleCell::datavalues->new(
	{ 'data_storage_spliced' => 1, 'file' => $outfile } );

$datavalues->{'use_this_sql'} = "select gname from genes";

$value = $datavalues->get_data_table_4_search(
	{
		'search_columns' => ['gname'],
		'where'          => [],
	}
);

is_deeply(
	$datavalues->{'complex_search'},
	"select gname from genes",
	"used the correct search '$datavalues->{'complex_search'}'"
);

#print "\$exp = ".root->print_perl_var_def($value->GetAsArray('gname') ).";\n";

$exp = [
	qw( fake AC132586.6 Myadm AC217111.2 Cacng8 AC164883.8 Fzd6 AC164883.6 AC164883.4 Slc25a32)
];

is_deeply( $value->GetAsArray('gname'), $exp, "right genes in the database" );

## last gene in folder A is 'Cacng8'
$exp = {
'TAAGCGTGTGGGTATG' => 4,
'CCGTGGACATAACCTG' => 16,
'GTGTTAGCAGAGCCAA' => 12,
};

$value = $datavalues->get_data_table_4_search(
	{
		'search_columns' => ['sname', 'value'],
		'where'          => [['gname', '=', 'my_value'] ],
	}, 'Cacng8'
);

is_deeply( $value->GetAsHash('sname', 'value' ), $exp, "right data for a gene in folder A");


## this is the entry in the test matrix B for gene 5 'AC164883.4'
#5 36 2 0
#5 37 4 0
#5 38 2 0
#5 39 2 0
#5 40 2 0
#5 41 2 0
#5 42 2 0
#5 43 2 0

## so now I want to get these values from the database...
$value = $datavalues->get_data_table_4_search(
	{
		'search_columns' => ['sname', 'value'],
		'where'          => [['gname', '=', 'my_value'] ],
	}, 'AC164883.4'
);

#print $value->AsString();
$exp = {
	'CCTAGCTTCAACCATG' => 2, 
	'GGGTTGCGTCTAGCCG' => 4,
	'GATGAGGAGCGCCTCA' => 2,
	'CGACCTTAGACTTTCG' => 2,
	'CAAGATCAGGCAAAGA' => 2,
	'CACCAGGGTGGTAACG' => 2,
	'ATCCACCTCGGGAGTA' => 2,
	'ATGAGGGCAGATGAGC' => 2,
};

is_deeply( $value->GetAsHash('sname', 'value' ), $exp, "right data for a gene in folder B");

#print "\$exp = ".root->print_perl_var_def( \@snames ).";\n";
