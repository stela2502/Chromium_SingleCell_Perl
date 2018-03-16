#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 403;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my (
	$value, @values, $exp, @options,  $R1, $R2,
	$I1,    $genome, $gtf, $coverage, $outpath,
);

my $exec = $plugin_path . "/../bin/10xpipeline.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/10xpipeline_debug";
if ( -d $outpath ) {

	#	system("rm -Rf $outpath/*");
}

chdir($outpath);

$R1 = "$plugin_path/data/10xpipeline/R1.txt";
$R2 = "$plugin_path/data/10xpipeline/R2.txt";
$I1 = "$plugin_path/data/10xpipeline/I1.txt";

ok( -f $R1, "R1 list file '$R1'" );
ok( -f $R2, "R2 list file '$R2'" );
ok( -f $I1, "I1 list file '$I1'" );

## based on the results this could be a usable gtf....
$gtf =
  "$plugin_path/data/10xpipeline/10x_gencode.vM16.chr_patch_hapl_scaff.gtf";

ok( -f $gtf, "gtf file '$gtf'" );

$coverage = "$plugin_path/data/10xpipeline/mm10_chrom_sizes.txt";

ok( -f $coverage, "coverage file '$coverage'" );

$genome = '~/lunarc/indicies/hisat2/mouse/mm10/genome'
  ;    ## OK that is very specific - should I cerate a minimal here?

@options =
  qw( A lsens2017-3-2 t 02:00:00 p dell min_UMIs 1 report_on gene_name w nodelist=ls2-n4 );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . " -options "
  . join( ' ', @options )
  . " -R1 $R1"
  . " -R2 $R2"
  . " -I1 $I1"
  . " -genome "
  . $genome
  . " -gtf "
  . $gtf . " -n 1"
  . " -coverage "
  . $coverage
  . " -fast_tmp "
  . "$plugin_path/data/output/10xpipeline/fast_tmp"
  . " -outpath "
  . $outpath
  . " -sname sampleX"
  . " -local"
  . " -debug";

sub fcontent {
	my $f = shift;
	return () unless ( -f $f );
	open( OO, "<$f" ) or die "I could not open the file '$f'\n$!\n";
	my @ret = map { chomp; $_ } <OO>;
	close(OO);
	@ret;
}

my $start = time;
print "I am starting cmd $cmd\n";
system($cmd );
my $duration = time - $start;
print "Run time: $duration s\n";

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

ok( -d $outpath, "outpath created" );
my $run_folder = $outpath . "/10xpipeline_run_files";
ok( -d $run_folder, "run filkes folder" );

## In the main run folder I get all the merged fastq files:

my @scripts1 = qw(
HJ2GYBGX5_Ctrl-LSK_S5_L001.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S7_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L001.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S5_L002.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S7_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L002.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S5_L003.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S7_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L003.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S5_L004.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S7_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S7_L004.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L001.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S8_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L001.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L002.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S8_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L002.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L002.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L003.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S8_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L003.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L003.annotated.fastq.sh
HJ2GYBGX5_Ctrl-LSK_S6_L004.annotated.fastq.sh  HJ2GYBGX5_Ctrl-LSK_S8_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S6_L004.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S8_L004.annotated.fastq.sh
);

foreach (@scripts1) {
	ok( -f "$run_folder/$_", "script file '$_'" );

	#print join("\n",  &fcontent("$run_folder/$_"));
	ok( scalar( grep( /SBATCH \-w/, &fcontent("$run_folder/$_") ) ) == 1,
		"script $_ contains -w SLURM option" );
}

## HISAT scripts:
my @hisat_scripts = qw(
  hisat2_run.sh
  HJ2GYBGX5_Ctrl-LSK_S5_L001.annotated.fastq.sh  HTL2KBGX3_Ctrl-LSK_S5_L001.annotated.fastq.sh
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

foreach (@hisat_scripts) {
	ok( -f "$run_folder/$_", "script file '$_'" );
	ok( scalar( grep( /SBATCH \-w/, &fcontent("$run_folder/$_") ) ) == 1,
		"script $_ contains -w SLURM option" );
}

my @hisat_fake_output =
  qw(FAKE_DEBUG_HJ2GYBGX5_Ctrl-LSK_S5_L001.sorted.bam  FAKE_DEBUG_HTL2KBGX3_Ctrl-LSK_S5_L001.sorted.bam
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

foreach (@hisat_fake_output) {
	ok( -f "$run_folder/HISAT2_mapped/$_", "fake bam file '$_'" );
}

## reshuffle output

$_ = 'reshuffled/sampleX_reshuffle_script_run.sh';

ok( -f "$run_folder/$_", "script file '$_'" );
ok( scalar( grep( /SBATCH \-w/, &fcontent("$run_folder/$_") ) ) == 1,
	"script $_ contains -w SLURM option" );

my @reshuffle = qw(

  chr1_sampleX.chr1:0-10000000.sorted.bam           chr1_sampleX.chr1:90000000-100000000.sorted.bam   chr2_sampleX.chr2:90000000-100000000.sorted.bam   chrX_sampleX.chrX:100000000-110000000.sorted.bam
  chr1_sampleX.chr1:100000000-110000000.sorted.bam  chr2_sampleX.chr2:0-10000000.sorted.bam           chr3_sampleX.chr3:0-10000000.sorted.bam           chrX_sampleX.chrX:10000000-20000000.sorted.bam
  chr1_sampleX.chr1:10000000-20000000.sorted.bam    chr2_sampleX.chr2:100000000-110000000.sorted.bam  chr3_sampleX.chr3:100000000-110000000.sorted.bam  chrX_sampleX.chrX:110000000-120000000.sorted.bam
  chr1_sampleX.chr1:110000000-120000000.sorted.bam  chr2_sampleX.chr2:10000000-20000000.sorted.bam    chr3_sampleX.chr3:10000000-20000000.sorted.bam    chrX_sampleX.chrX:120000000-130000000.sorted.bam
  chr1_sampleX.chr1:120000000-130000000.sorted.bam  chr2_sampleX.chr2:110000000-120000000.sorted.bam  chr3_sampleX.chr3:110000000-120000000.sorted.bam  chrX_sampleX.chrX:130000000-140000000.sorted.bam
  chr1_sampleX.chr1:130000000-140000000.sorted.bam  chr2_sampleX.chr2:120000000-130000000.sorted.bam  chr3_sampleX.chr3:120000000-130000000.sorted.bam  chrX_sampleX.chrX:140000000-150000000.sorted.bam
  chr1_sampleX.chr1:140000000-150000000.sorted.bam  chr2_sampleX.chr2:130000000-140000000.sorted.bam  chr3_sampleX.chr3:130000000-140000000.sorted.bam  chrX_sampleX.chrX:150000000-160000000.sorted.bam
  chr1_sampleX.chr1:150000000-160000000.sorted.bam  chr2_sampleX.chr2:140000000-150000000.sorted.bam  chr3_sampleX.chr3:140000000-150000000.sorted.bam  chrX_sampleX.chrX:160000000-170000000.sorted.bam
  chr1_sampleX.chr1:160000000-170000000.sorted.bam  chr2_sampleX.chr2:150000000-160000000.sorted.bam  chr3_sampleX.chr3:150000000-160000000.sorted.bam  chrX_sampleX.chrX:170000000-171031299.sorted.bam
  chr1_sampleX.chr1:170000000-180000000.sorted.bam  chr2_sampleX.chr2:160000000-170000000.sorted.bam  chr3_sampleX.chr3:160000000-160039680.sorted.bam  chrX_sampleX.chrX:20000000-30000000.sorted.bam
  chr1_sampleX.chr1:180000000-190000000.sorted.bam  chr2_sampleX.chr2:170000000-180000000.sorted.bam  chr3_sampleX.chr3:20000000-30000000.sorted.bam    chrX_sampleX.chrX:30000000-40000000.sorted.bam
  chr1_sampleX.chr1:190000000-195471971.sorted.bam  chr2_sampleX.chr2:180000000-182113224.sorted.bam  chr3_sampleX.chr3:30000000-40000000.sorted.bam    chrX_sampleX.chrX:40000000-50000000.sorted.bam
  chr1_sampleX.chr1:20000000-30000000.sorted.bam    chr2_sampleX.chr2:20000000-30000000.sorted.bam    chr3_sampleX.chr3:40000000-50000000.sorted.bam    chrX_sampleX.chrX:50000000-60000000.sorted.bam
  chr1_sampleX.chr1:30000000-40000000.sorted.bam    chr2_sampleX.chr2:30000000-40000000.sorted.bam    chr3_sampleX.chr3:50000000-60000000.sorted.bam    chrX_sampleX.chrX:60000000-70000000.sorted.bam
  chr1_sampleX.chr1:40000000-50000000.sorted.bam    chr2_sampleX.chr2:40000000-50000000.sorted.bam    chr3_sampleX.chr3:60000000-70000000.sorted.bam    chrX_sampleX.chrX:70000000-80000000.sorted.bam
  chr1_sampleX.chr1:50000000-60000000.sorted.bam    chr2_sampleX.chr2:50000000-60000000.sorted.bam    chr3_sampleX.chr3:70000000-80000000.sorted.bam    chrX_sampleX.chrX:80000000-90000000.sorted.bam
  chr1_sampleX.chr1:60000000-70000000.sorted.bam    chr2_sampleX.chr2:60000000-70000000.sorted.bam    chr3_sampleX.chr3:80000000-90000000.sorted.bam    chrX_sampleX.chrX:90000000-100000000.sorted.bam
  chr1_sampleX.chr1:70000000-80000000.sorted.bam    chr2_sampleX.chr2:70000000-80000000.sorted.bam    chr3_sampleX.chr3:90000000-100000000.sorted.bam
  chr1_sampleX.chr1:80000000-90000000.sorted.bam    chr2_sampleX.chr2:80000000-90000000.sorted.bam    chrX_sampleX.chrX:0-10000000.sorted.bam
);

foreach (@reshuffle) {
	ok( -f "$run_folder/reshuffled/$_", "reshuffled output bam '$_'" );
}

my @QuantD = qw(
  chr1:0-10000000_sampleX           chr1:150000000-160000000_sampleX  chr1:40000000-50000000_sampleX   chr2:100000000-110000000_sampleX  chr2:160000000-170000000_sampleX  chr2:60000000-70000000_sampleX
  chr1:100000000-110000000_sampleX  chr1:160000000-170000000_sampleX  chr1:50000000-60000000_sampleX   chr2:10000000-20000000_sampleX    chr2:170000000-180000000_sampleX  chr2:70000000-80000000_sampleX
  chr1:10000000-20000000_sampleX    chr1:170000000-180000000_sampleX  chr1:60000000-70000000_sampleX   chr2:110000000-120000000_sampleX  chr2:180000000-182113224_sampleX  chr2:80000000-90000000_sampleX
  chr1:110000000-120000000_sampleX  chr1:180000000-190000000_sampleX  chr1:70000000-80000000_sampleX   chr2:120000000-130000000_sampleX  chr2:20000000-30000000_sampleX    chr2:90000000-100000000_sampleX
  chr1:120000000-130000000_sampleX  chr1:190000000-195471971_sampleX  chr1:80000000-90000000_sampleX   chr2:130000000-140000000_sampleX  chr2:30000000-40000000_sampleX
  chr1:130000000-140000000_sampleX  chr1:20000000-30000000_sampleX    chr1:90000000-100000000_sampleX  chr2:140000000-150000000_sampleX  chr2:40000000-50000000_sampleX
  chr1:140000000-150000000_sampleX  chr1:30000000-40000000_sampleX    chr2:0-10000000_sampleX          chr2:150000000-160000000_sampleX  chr2:50000000-60000000_sampleX
);

foreach (@QuantD) {
	ok( -d "$run_folder/$_", "Quantify output folders '$_'" );
}


foreach my $path (@QuantD) {
	map { ok( -f "$run_folder/$path/$_", "opath $path ofile $_" ) }
	  qw( barcodes.tsv genes.tsv matrix.mtx);
}

ok( -f "$outpath/sampleX_FAKE_DEBUG.sqlite", "main outfile" );
