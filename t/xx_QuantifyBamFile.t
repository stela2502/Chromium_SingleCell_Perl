#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 20;
use stefans_libs::flexible_data_structures::data_table;
use stefans_libs::result_table;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $gtf_file, $infile, $outfile, @options, );

my $exec = $plugin_path . "/../bin/QuantifyBamFile.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/QuantifyBamFile";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}
my $fastqPath = "$plugin_path/data/QuantifBamFile_Fake_fastq_folder";
ok ( -d $fastqPath, "the fastqPath");
ok ( -f $fastqPath."/.passing_samples.txt", "the hidden file .passing_samples.txt");

$gtf_file = "$plugin_path/data/Spliced_Reads.gtf";
ok ( -f $gtf_file, 'test gtf file exists');

$infile = "$plugin_path/data/Spliced_Reads.bam";
ok ( -f $infile, 'test infile exists');

$outfile = "$outpath/test";

#improvement by using UMIs to merge samples:
#no usage: 61773 samples and 4 genes
#with umi: 25261 and 4 genes
#with umi and sample merge: 16224 samples and 4 genes detected

# with only 4 geens we are down to only half the samples.

@options = ('min_UMIs', 10 );

my $cmd =
    "perl "
    . "-I $plugin_path/../lib "
    . "-I $plugin_path/../../stefans_libs-BAMfile/lib/ "
    . "-I ../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../stefans_libs-GenomeDB/lib "
    ."$exec "
. " -gtf_file " . $gtf_file 
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -options " . join(' ', @options )
. " -fastqPath $fastqPath"
. " -sampleID TestSample"
#. " -debug"
;

my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";

## now that was simple...
## the outfiles are
my @outfiles = map{ $outfile."/$_" } qw(barcodes.tsv  genes.tsv  matrix.mtx);
map { ok ( -f $_, "outfile $_") } @outfiles;

my $result_table = stefans_libs::result_table->new();
$result_table -> import_tables ( $outfile );
$result_table ->line_separator(';');

#for ( my $i = 0; $i < 2; $i++ ){
#	my %t;
#	@t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[$i]};
#	print "$i : \$exp = ".root->print_perl_var_def( \%t ).";\n";
#	
#}

$exp = {
  'ATAGACCGTGGTGTAG' => '2',
  'CACAAACCACTTCTGC' => '2',
  'CCGTGGACATAACCTG' => undef,
  'CCGTGGACATAACCTG spliced' => undef,
  'GTCATTTAGTCGTTTG' => '2',
  'GTCCTCATCTGCAGTA' => '2',
  'Gene_ID' => 'ENSMUSG00000108159.1',
  'TAAGCGTGTGGGTATG' => undef,
  'TAAGCGTGTGGGTATG spliced' => undef,
  'TACAGTGGTTATTCTC' => '2',
  'TGCGTGGAGCTACCGC' => '6'
};

my $t;
@$t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[0]};

is_deeply( $t, $exp, "data line for ENSMUSG00000108159.1" ); #109020	111963

$exp = {
  'ATAGACCGTGGTGTAG' => undef,
  'CACAAACCACTTCTGC' => undef,
  'CCGTGGACATAACCTG' => '2',
  'CCGTGGACATAACCTG spliced' => '2',
  'GTCATTTAGTCGTTTG' => undef,
  'GTCCTCATCTGCAGTA' => undef,
  'Gene_ID' => 'ENSMUSG00000106728.3',
  'TAAGCGTGTGGGTATG' => '3',
  'TAAGCGTGTGGGTATG spliced' => '3',
  'TACAGTGGTTATTCTC' => undef,
  'TGCGTGGAGCTACCGC' => undef
};

$t = undef;

@$t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[1]};

is_deeply( $t, $exp, "data line for ENSMUSG00000106728.3" ); # 389199	414164

## restart #1
system("rm -Rf $outpath/*");
$outfile = "$outpath/chrJH792828_1_and_something_else.sqlite";
$cmd =
    "perl "
    . "-I $plugin_path/../lib "
    . "-I $plugin_path/../../stefans_libs-BAMfile/lib/ "
    . "-I ../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../stefans_libs-GenomeDB/lib "
    ."$exec "
. " -gtf_file " . $gtf_file 
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -options " . join(' ', @options )
. " -fastqPath $fastqPath"
. " -sampleID TestSample"
. " -drop_chr chrJH792828.1"
#. " -debug"
;

$start = time;
system( $cmd );
 $duration = time - $start;
print "Execution time: $duration s\n";

## now that was simple...
## the outfiles are
$outfile = "$outpath/chrJH792828_1_and_something_else";
@outfiles = map{ $outfile."/$_" } qw(barcodes.tsv  genes.tsv  matrix.mtx);
map { ok ( -f $_, "outfile $_") } @outfiles;

$result_table = stefans_libs::result_table->new();
$result_table -> import_tables ( $outfile );
$result_table ->line_separator(';');

#for ( my $i = 0; $i < 2; $i++ ){
#	my %t;
#	@t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[$i]};
#	print "$i : \$exp = ".root->print_perl_var_def( \%t ).";\n";
#	
#}

$exp = {
  'ATAGACCGTGGTGTAG' => '2',
  'CACAAACCACTTCTGC' => '2',
  'CCGTGGACATAACCTG' => undef,
  'CCGTGGACATAACCTG spliced' => undef,
  'GTCATTTAGTCGTTTG' => '2',
  'GTCCTCATCTGCAGTA' => '2',
  'Gene_ID' => 'ENSMUSG00000108159.1',
  'TAAGCGTGTGGGTATG' => undef,
  'TAAGCGTGTGGGTATG spliced' => undef,
  'TACAGTGGTTATTCTC' => '2',
  'TGCGTGGAGCTACCGC' => '6'
};

$t = undef;
@$t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[0]};

is_deeply( $t, $exp, "drop_chr data line for ENSMUSG00000108159.1" ); #109020	111963

$exp = {
  'ATAGACCGTGGTGTAG' => undef,
  'CACAAACCACTTCTGC' => undef,
  'CCGTGGACATAACCTG' => '2',
  'CCGTGGACATAACCTG spliced' => '2',
  'GTCATTTAGTCGTTTG' => undef,
  'GTCCTCATCTGCAGTA' => undef,
  'Gene_ID' => 'ENSMUSG00000106728.3',
  'TAAGCGTGTGGGTATG' => '3',
  'TAAGCGTGTGGGTATG spliced' => '3',
  'TACAGTGGTTATTCTC' => undef,
  'TGCGTGGAGCTACCGC' => undef
};

$t = undef;

@$t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[1]};

is_deeply( $t, $exp, "drop_chr data line for ENSMUSG00000106728.3" ); # 389199	414164


## restart
system("rm -Rf $outpath/*");
$outfile = "$outpath/chrJH792828_1:108000-114000";
$cmd =
    "perl "
    . "-I $plugin_path/../lib "
    . "-I $plugin_path/../../stefans_libs-BAMfile/lib/ "
    . "-I ../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../stefans_libs-GenomeDB/lib "
    ."$exec "
. " -gtf_file " . $gtf_file 
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -options " . join(' ', @options )
. " -fastqPath $fastqPath"
. " -sampleID TestSample"
. " -drop_chr chrJH792828.1:108000-114000" ## only ENSMUSG00000108159.1
#. " -debug"
;
system( $cmd );
$result_table = stefans_libs::result_table->new();
ok ( -d $outfile, "the chr slice#1 outpath");
$result_table -> import_tables ( $outfile );
#print "\$exp = ".root->print_perl_var_def($value ).";\n";

$exp = {
  'ATAGACCGTGGTGTAG' => '2',
  'CACAAACCACTTCTGC' => '2',
  'GTCATTTAGTCGTTTG' => '2',
  'GTCCTCATCTGCAGTA' => '2',
  'Gene_ID' => 'ENSMUSG00000108159.1',
  'TACAGTGGTTATTCTC' => '2',
  'TGCGTGGAGCTACCGC' => '6'
};

$t = undef;
@$t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[0]};

#print "\$exp = ".root->print_perl_var_def($t ).";\n";
is_deeply( $t, $exp, "data line for ENSMUSG00000108159.1 only" ); #109020	111963

ok ($result_table->Lines() == 1, "Only one gene quantified." );

## restart
system("rm -Rf $outpath/*");
$outfile = "$outpath/chrJH792828_1:389100-415000";

$cmd =
    "perl "
    . "-I $plugin_path/../lib "
    . "-I $plugin_path/../../stefans_libs-BAMfile/lib/ "
    . "-I ../../Stefans_Libs_Essentials/Stefans_Libs_Essentials/lib "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../SLURM_bioinformatics_scripts/lib/ "
    . "-I ../../stefans_libs-GenomeDB/lib "
    ."$exec "
. " -gtf_file " . $gtf_file 
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -options " . join(' ', @options )
. " -fastqPath $fastqPath"
. " -sampleID TestSample"
. " -drop_chr chrJH792828.1:389100-415000" ## only ENSMUSG00000106728.3 389199	414164
#. " -debug"
;
system( $cmd );
$result_table = stefans_libs::result_table->new();
$result_table -> import_tables ( $outfile );

$t = undef;
@$t{@{$result_table->{header}}} = @{@{$result_table->{'data'}}[0]};

$exp = {
  'CCGTGGACATAACCTG' => '2',
  'CCGTGGACATAACCTG spliced' => '2',
  'Gene_ID' => 'ENSMUSG00000106728.3',
  'TAAGCGTGTGGGTATG' => '3',
  'TAAGCGTGTGGGTATG spliced' => '3',
};

is_deeply( $t, $exp, "data line for ENSMUSG00000108159.1 only" ); #109020	111963

ok ($result_table->Lines() == 1, "Only one gene quantified." );



