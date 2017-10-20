#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-10-16 Stefan Lang

  This program is free software; you can redistribute it 
  and/or modify it under the terms of the GNU General Public License 
  as published by the Free Software Foundation; 
  either version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License 
  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1 CREATED BY
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit f411d5d0199ebcbda88fa8129a8369e8a6c53b05
   

=head1  SYNOPSIS

    10xpipeline.pl
       -fastq     :<please add some info!> you can specify more entries to that
       -gtf       :<please add some info!>
       -coverage       :<please add some info!>
       -genome       :<please add some info!>
       -outpath       :<please add some info!>
       -options     :the full set of SLURM options:
                       A lsens2017-3-2 N 1 t 02:00:00

       -fast_tmp :the nodes should have a fast local tmp folder that could be used for 
                  the intermediate files (default '$SNIC_TMP')
                  
       -sname    :the sample name analyzed here
                  
       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Run the whole perl pipeline based on HISAT2 sidestepping most of Illuminas own programm.
  
  The fastq files need to be obtained from the 10x data using the mkfastq command.
  All three I1 R1 and R2 are used here.

  To get further help use '10xpipeline.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::SLURM;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help,     $debug,  $database, @fastq,   $gtf,     $fast_tmp,
	$coverage, $genome, $outpath,  $options, @options, $sname
);

Getopt::Long::GetOptions(
	"-fastq=s{,}"   => \@fastq,
	"-gtf=s"        => \$gtf,
	"-coverage=s"   => \$coverage,
	"-genome=s"     => \$genome,
	"-fast_tmp=s"   => \$fast_tmp,
	"-outpath=s"    => \$outpath,
	"-options=s{,}" => \@options,
	"-sname=s"      => \$sname,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $fastq[0] ) {
	$error .= "the cmd line switch -fastq is undefined!\n";
}
unless ( defined $gtf ) {
	$error .= "the cmd line switch -gtf is undefined!\n";
}
unless ( defined $coverage ) {
	$error .= "the cmd line switch -coverage is undefined!\n";
}
unless ( defined $genome ) {
	$error .= "the cmd line switch -genome is undefined!\n";
}
unless ( defined $outpath ) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $options[0] ) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ($fast_tmp) {
	$fast_tmp = '$SNIC_TMP';
}
unless ($sname) {
	$warn .= "The sname was set to 'sampleX'\n";
	$sname = 'sampleX';
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	print "$errorMessage.\n";
	pod2usage( q(-verbose) => 1 );
}

### initialize default options:

#$options->{'n'} ||= 10;

###

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/10xpipeline.pl';
$task_description .= ' -fastq "' . join( '" "', @fastq ) . '"'
  if ( defined $fastq[0] );
$task_description .= " -gtf '$gtf'"           if ( defined $gtf );
$task_description .= " -coverage '$coverage'" if ( defined $coverage );
$task_description .= " -genome '$genome'"     if ( defined $genome );
$task_description .= " -outpath '$outpath'"   if ( defined $outpath );
$task_description .= ' -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################
mkdir($outpath) unless ( -d $outpath );
open( LOG, ">$outpath/" . $$ . "_10xpipeline.pl.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

my $fm = root->filemap( "$outpath/" . $$ . "_10xpipeline.pl.log");
$outpath = $fm->{'path'};

## default slurm options
$options->{'A'} ||= 'lsens2017-3-2';
$options->{'n'} ||= 1;
$options->{'N'} ||= 1;
$options->{'t'} ||= '02:00:00';

my $SLURM = stefans_libs::SLURM->new($options, 0);
$SLURM->{'debug'} = 1 if ($debug);

## Do whatever you want!

## first we need to run the SplitToCells.pl script. this should definitely be run on a blade.
my ( $cmd, $sumFastq ) = &SplitToCell();
my $SLURM_id = 0;

print "$cmd"."\n$sumFastq\n";

if ( ! -f $sumFastq  and ! $debug ) {

	$SLURM_id = $SLURM->run( $cmd, $sumFastq );
}else {
	warn "OUtfile present - SplitToCell not re-run\n";
}

&wait_for_PID($SLURM_id);
## now map the sample

#hisat2_run_aurora.pl 
#print root::get_hashEntries_as_string($options,3, "the oSLURM ptions:");

my $hisat2 =
"hisat2_run_aurora.pl -files $sumFastq -outpath $outpath/HISAT2_mapped/ -options n 5 partitition $options->{'p'} A ".$options->{'A'};
$hisat2 .=
" -genome $genome -coverage $coverage -bigwigTracks $outpath/hisat2_$sname.html";
$hisat2 .= " -fast_tmp '$fast_tmp'";

if ($debug) {
	print $hisat2."\n";
}
else {
	print $hisat2."\n";
	open( RUN, $hisat2 . " | " );
	$SLURM_id = 0;
	foreach (<RUN>) {
		if ( $_ = /Submitted batch job (\d+)/ ) {
			$SLURM_id = $1;
			last;
		}
	}

	close(RUN);
	system($cmd );
}

&wait_for_PID($SLURM_id);

## now I need to Quantify this sample

#QuantifyBamFile.pl -infile ./HVJ5GBGX2/outs/fastq_path/Tobias1/s1/CLP_S1_L001_R2_001.fastq.gz -outfile perl_requant/CLP_S1_L001_R2 -gtf_file ~/lunarc/genomes/mouse/mm10/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf
# This takes so much memory I need to run it on the frontend ;-(

$cmd = QuantifyBamFile();
if ($debug) {
	print $cmd."\n";
}
else {
	print $cmd."\n";
	system($cmd ) unless ( -f "$outpath/$sname.original_merged.db" );
}

sub QuantifyBamFile {
	my $bam_file =
	  "$outpath/HISAT2_mapped/$sname.annotated.fastq_hisat.sorted.bam";
	my $cmd = 'QuantifyBamFile.pl ';
	$cmd .= " -infile $bam_file";
	$cmd .= " -gtf_file $gtf";
	$cmd .= " -outfile $outpath/quant_$sname";
	return $cmd;
}

sub wait_for_PID {
	my $ID = shift;
	if ( $ID > 0 ) {    ## else the file has already been prepared
		while ( !$SLURM->pids_finished($ID) ) {
			sleep(10);
		}
	}
}

sub SplitToCell {
	my $cmd = 'SplitToCells.pl';

#SplitToCells.pl -I1 ./HVJ5GBGX2/outs/fastq_path/Tobias1/s1/CLP_S1_L001_I1_001.fastq.gz -R1 ./HVJ5GBGX2/outs/fastq_path/Tobias1/s1/CLP_S1_L001_R1_001.fastq.gz -R2
#./HVJ5GBGX2/outs/fastq_path/Tobias1/s1/CLP_S1_L001_R2_001.fastq.gz -outpath Joined_Tina_fastq -options oname Tina_CLP_singleCells
#gives file 'Joined_Tina_fastq/Tina_CLP_singleCells.annotated.fastq.gz'
	my $i;
	map {
		if ( $_ =~ m/\d\d\d_([IR][12])_\d\d\d/ ) {
			$cmd .= " -$1 $_";
			$i->{$1}||=0;
			$i->{$1} ++;
		}else {
			warn "$_ does not match \\d\\d\\d_[IR][12]_\\d\\d\\d?\n";
		}
	} @fastq;
	unless ( join( ' ', sort keys %$i ) eq "I1 R1 R2"
		and $i->{'I1'} + $i->{'R1'} + $i->{'R2'} == 3 )
	{
		Carp::confess(
"Sorry, but SplitToCell.pl requires exactly one I1, one R1 and one R2 fastq file\n'".
 join( ' ', sort keys %$i ). "' should be 'I1 R1 R2' (". ($i->{'I1'} + $i->{'R1'} + $i->{'R2'}). ") should be 3\n"
		);
	}
	$cmd .= " -outpath $fast_tmp";
	$cmd .= " -options oname $sname";
	$cmd .= "\ncp $fast_tmp/$sname.annotated.fastq.gz $outpath\n";
	$cmd .= "cp $fast_tmp/*.log $outpath\n";
	

	return $cmd, join( "/", $outpath, $sname . ".annotated.fastq.gz" );
}

