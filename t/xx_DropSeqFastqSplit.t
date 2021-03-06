#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 3;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $fastq, $outpath, );

my $exec = $plugin_path . "/../bin/DropSeqFastqSplit.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/DropSeqFastqSplit";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$fastq = "$plugin_path/data/simpleFastq.fastq.gz";

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . " -fastq "
  . $fastq
  . " -outpath "
  . $outpath
  . " -debug";
my $start = time;
system($cmd );
my $duration = time - $start;
print "Run time: $duration s\n";

use stefans_libs::FastqFile;

my $worker = stefans_libs::FastqFile->new();
my $ofile = File::Spec->catfile( $outpath, 'simpleFastq.fastq.gz' );
ok( -f $ofile, "outfile" );

my $function = sub {
	my ( $worker, $entry ) = @_;
	push( @values, join( "\n", @{ $entry->{'data'} } ) );
	$entry->clear();
};

$worker->filter_file( $function, $ofile );
#print "\$exp = " . root->print_perl_var_def( \@values ) . ";\n";

$exp = [ '@SRR6365735.76_N_N_CGGTTAAGATAC_AGGAGAAG
CATTTAAAANATAAANCANNGTANAAGTCACCAGATACTACTATCATTAT
+
///A//A//#EE//E#//##///#///A/<////<6A/////////////', '@SRR6365735.79_N_N_CCCCGTACTCTG_CGTAGATA
CCCCCGTACTCTGCGTNGATACCACTGCTTCCGCGGACAGGCGTGTAGA
+
///<//////A///<#///6//E/E/6<////<</A//6///E/EA///', '@SRR6365735.80_N_N_CTCTGCGTTGAT_ACCACTGC
CTCCTAGGCCACAGTNGTACTCTGCGTTGATACCACTGCTTCCGCGGAC
+
/6///A</E/E///<#/<E///A/</EE///</<A/A/A///A/<A///' ];

is_deeply( \@values, $exp, "right fastq entries" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";
