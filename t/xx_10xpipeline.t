#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 88;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, $R1, $R2, $I1, $genome, $gtf, $coverage, $outpath, );

my $exec = $plugin_path . "/../bin/10xpipeline.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/10xpipeline";
if ( -d $outpath ) {
	#system("rm -Rf $outpath/*");
}

$R1 = "$plugin_path/data/10xpipeline/R1.txt";
$R2 = "$plugin_path/data/10xpipeline/R2.txt";
$I1 = "$plugin_path/data/10xpipeline/I1.txt";

ok( -f $R1, "R1 list file '$R1'");
ok( -f $R2, "R2 list file '$R2'");
ok( -f $I1, "I1 list file '$I1'");

## based on the results this could be a usable gtf....
$gtf = "$plugin_path/data/10xpipeline/10x_gencode.vM16.chr_patch_hapl_scaff.gtf";

ok( -f $gtf, "gtf file '$gtf'");

$coverage = "$plugin_path/data/10xpipeline/mm10_chrom_sizes.txt";

ok( -f $coverage, "coverage file '$coverage'");

$genome = '~/lunarc/indicies/hisat2/mouse/mm10/genome'; ## OK that is very specific - should I cerate a minimal here?

@options = qw( A lsens2017-3-2 t 02:00:00 p dell min_UMIs 1 report_on gene_name);

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -options " . join(' ', @options )
. " -R1 $R1" 
. " -R2 $R2" 
. " -I1 $I1" 
. " -genome " . $genome 
. " -gtf " . $gtf 
. " -coverage " . $coverage 
. " -fast_tmp ". "$plugin_path/data/output/10xpipeline/fast_tmp"
. " -outpath " . $outpath 
. " -sname sampleX"
#. " -debug"
;



my $start = time;
system( $cmd );
my $duration = time - $start;
print "Run time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";

ok ( -d $outpath , "outpath created" );
my $run_folder = $outpath."/10xpipeline_run_files";
ok ( -d $run_folder, "run filkes folder" );

## In the main run folder I get all the merged fastq files:

my @mergedFastq = qw( HJ2GYBGX5_Ctrl-LSK_S5_L001.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S5_L001.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S5_L002.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S5_L002.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S5_L003.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S5_L003.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S5_L004.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S5_L004.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S6_L001.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S6_L001.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S6_L002.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S6_L002.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S6_L003.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S6_L003.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S6_L004.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S6_L004.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S7_L001.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S7_L001.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S7_L002.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S7_L002.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S7_L003.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S7_L003.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S7_L004.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S7_L004.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S8_L001.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S8_L001.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S8_L002.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S8_L002.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S8_L003.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S8_L003.annotated.fastq.gz
HJ2GYBGX5_Ctrl-LSK_S8_L004.annotated.fastq.gz  HTL2KBGX3_Ctrl-LSK_S8_L004.annotated.fastq.gz);

foreach  ( @mergedFastq ) {
	ok ( -f "$run_folder/$_", "merged fastq '$_'")
}

ok ( -d "$run_folder/HISAT2_mapped", "hisat2 mapped folder" );

my @HISAT_results = qw( HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S5_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S5_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S5_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S5_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S6_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S6_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S6_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S6_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S7_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S7_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S7_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S7_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S8_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S8_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S8_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HJ2GYBGX5_Ctrl-LSK_S8_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S5_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S5_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S5_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S5_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S6_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S6_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S6_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S6_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S7_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S7_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S7_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S7_L004.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S8_L001.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S8_L002.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S8_L003.annotated.fastq_hisat.sorted.bam
HISAT2_mapped/HTL2KBGX3_Ctrl-LSK_S8_L004.annotated.fastq_hisat.sorted.bam
 ); 

foreach ( @HISAT_results ) {
 	ok ( -f "$run_folder/$_", "hisat2 mapped file $_" );
}

ok ( -d "$run_folder/reshuffled", "reshuffled folder");

my @reshuffled_bams = qw( chr1_sampleX.sorted.bam  chr2_sampleX.sorted.bam  chr3_sampleX.sorted.bam  chrX_sampleX.sorted.bam
); # using the smaller set here ;-)

foreach ( @reshuffled_bams ) {
 	ok ( -f "$run_folder/reshuffled/$_", "reshuffled bam file $_" );
}

my @quant_output = qw( sampleX.chr1_sampleX sampleX.chr2_sampleX sampleX.chr3_sampleX sampleX.chrX_sampleX  );

foreach ( @quant_output[0,1] ) {
	ok ( -d "$run_folder/$_", "quantification for chr $_ (folder)");
	foreach my $f ( qw(barcodes.tsv  genes.tsv  matrix.mtx) ) {
		ok ( -f "$run_folder/$_/$f", "quantification for chr $_ outfile $f");
	} 
}

foreach ( @quant_output[2,3] ) {
	ok ( -d "$run_folder/$_", "quantification for chr $_ (folder)");
	# files are missing as no reads mapped to any features on these.
#	foreach my $f ( qw(barcodes.tsv  genes.tsv  matrix.mtx) ) {
#		ok ( -f "$run_folder/$_/$f", "quantification for chr $_ outfile $f");
#	} 
}

## now I really want to look at the result data:


