#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 164;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, $R1, $R2, $I1, $genome, $gtf, $coverage, $outpath, );

my $exec = $plugin_path . "/../bin/10xpipeline.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/10xpipeline_debug";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
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

@options = qw( A lsens2017-3-2 t 02:00:00 p dell min_UMIs 1 report_on gene_name w nodelist=ls2-n4 );

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
. " -debug"
;

sub fcontent { 
	my $f = shift;
	return () unless ( -f $f );
	open ( OO, "<$f" ) or die "I could not open the file '$f'\n$!\n";
	my @ret = map { chomp; $_} <OO>;
	close ( OO );
	@ret;
}

my $start = time;
system( $cmd );
my $duration = time - $start;
print "Run time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";

ok ( -d $outpath , "outpath created" );
my $run_folder = $outpath."/10xpipeline_run_files";
ok ( -d $run_folder, "run filkes folder" );

## In the main run folder I get all the merged fastq files:


my @scripts1 = qw(
Ctrl-LSK_S5_L001_R1_001.fastq.sh  Ctrl-LSK_S6_L003_R1_001.fastq.sh  Ctrl-LSK_S8_L001_R1_001.fastq.sh
Ctrl-LSK_S5_L002_R1_001.fastq.sh  Ctrl-LSK_S6_L004_R1_001.fastq.sh  Ctrl-LSK_S8_L002_R1_001.fastq.sh
Ctrl-LSK_S5_L003_R1_001.fastq.sh  Ctrl-LSK_S7_L001_R1_001.fastq.sh  Ctrl-LSK_S8_L003_R1_001.fastq.sh
Ctrl-LSK_S5_L004_R1_001.fastq.sh  Ctrl-LSK_S7_L002_R1_001.fastq.sh  Ctrl-LSK_S8_L004_R1_001.fastq.sh
Ctrl-LSK_S6_L001_R1_001.fastq.sh  Ctrl-LSK_S7_L003_R1_001.fastq.sh
Ctrl-LSK_S6_L002_R1_001.fastq.sh  Ctrl-LSK_S7_L004_R1_001.fastq.sh
);

foreach  ( @scripts1 ) {
	ok ( -f "$run_folder/$_", "script file '$_'");
	#print join("\n",  &fcontent("$run_folder/$_"));
	ok (scalar( grep( /SBATCH \-w/, &fcontent("$run_folder/$_")) ) == 1, "script $_ contains -w SLURM option" );
}

## HISAT scripts:
my @hisat_scripts = qw( HJ2GYBGX5_Ctrl-LSK_S5_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L001.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S5_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L002.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S5_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L003.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S5_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L004.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L001.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L002.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L003.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L004.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S7_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L001.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S7_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L002.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S7_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L003.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S7_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L004.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S8_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L001.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S8_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L002.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S8_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L003.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S8_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L004.annotated.fastq.sh );

foreach  ( @hisat_scripts ) {
	ok ( -f "$run_folder/$_", "script file '$_'");
	ok (scalar( grep( /SBATCH \-w/, &fcontent("$run_folder/$_")) ) == 1, "script $_ contains -w SLURM option" );
}

my @hisat_fake_output = qw(FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S5_L001.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S5_L001.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S5_L002.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S5_L002.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S5_L003.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S5_L003.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S5_L004.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S5_L004.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S6_L001.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S6_L001.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S6_L002.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S6_L002.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S6_L003.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S6_L003.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S6_L004.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S6_L004.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S7_L001.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S7_L001.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S7_L002.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S7_L002.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S7_L003.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S7_L003.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S7_L004.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S7_L004.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S8_L001.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S8_L001.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S8_L002.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S8_L002.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S8_L003.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S8_L003.sorted.bam
FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S8_L004.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S8_L004.sorted.bam
);

foreach  ( @hisat_fake_output ) {
	ok ( -f "$run_folder/HISAT2_mapped/$_", "fake bam file '$_'");
}

## reshuffle output

$_ = 'reshuffled/sampleX_reshuffle.sorted.sh';

ok ( -f "$run_folder/$_", "script file '$_'");
ok (scalar( grep( /SBATCH \-w/, &fcontent("$run_folder/$_")) ) == 1, "script $_ contains -w SLURM option" );

my @reshuffle = qw(chr2_sampleX.sorted.bam  chrX_sampleX.sorted.bam
chr1_sampleX.sorted.bam         chr3_sampleX.sorted.bam  sampleX_reshuffle.sorted.sh
 );

foreach  ( @reshuffle ) {
	ok ( -f "$run_folder/reshuffled/$_", "reshuffled output bam '$_'")
}

my @QuantD = qw(
sampleX.chrX_sampleX
sampleX.chr2_sampleX
sampleX.chr1_sampleX
sampleX.chr3_sampleX
);


foreach  ( @QuantD ) {
	ok ( -d "$run_folder/$_", "Quantify output folders '$_'")
}

my @QuantF = qw(chrX_sampleX_sampleX.sh
chr2_sampleX_sampleX.sh
chr1_sampleX_sampleX.sh
chr3_sampleX_sampleX.sh
sampleX.chrX_sampleX/genes.tsv
sampleX.chrX_sampleX/barcodes.tsv
sampleX.chrX_sampleX/matrix.mtx
sampleX.chr2_sampleX/genes.tsv
sampleX.chr2_sampleX/barcodes.tsv
sampleX.chr2_sampleX/matrix.mtx
sampleX.chr1_sampleX/genes.tsv
sampleX.chr1_sampleX/barcodes.tsv
sampleX.chr1_sampleX/matrix.mtx
sampleX.chr3_sampleX/genes.tsv
sampleX.chr3_sampleX/barcodes.tsv
sampleX.chr3_sampleX/matrix.mtx
);

foreach  ( @QuantF ) {
	ok ( -f "$run_folder/$_", "Quantify output files '$_'")
}

ok ( -f "$outpath/sampleX_FAKE_DEBUG.sqlite", "main outfile" );
