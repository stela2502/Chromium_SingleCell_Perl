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
       -R1       :All R1 files for the sample
       -R2       :All R2 files for the sample
       -I1       :All I1 files for the sample
       -gtf      :the gtf file containing the genes to quantify on
       -coverage :the chromosmome length file
       -genome   :the hisat2 index to use
       -outpath  :the path all outfiles should be written to
       -n        :the number of cores to use (default 3)

       -fast_tmp :the nodes should have a fast local tmp folder that can be used for 
                  the intermediate files (default '$SNIC_TMP')
                  
       -sname    :the sample name analyzed here
       
       -options  :all options have to be given as key<space>value combinations
                  like -options A lsens2017-3-2 t 02:00:00 p dell
           A     :the slurm sbatch A option
           t     :the slurm time option
           p     :the slurm partition option
        min_UMIs :the QuantifyBamFiles min_UMIs option to identify samples (default 100)
              
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
	$help,     $debug,  $database, $gtf,     $fast_tmp,
	$coverage, $genome, $outpath,  $options, @options, @R1, @R2, @I1,$sname, $n
);

Getopt::Long::GetOptions(
    "-R1=s{,}"           => \@R1,
    "-R2=s{,}"           => \@R2,
    "-I1=s{,}"           => \@I1,
	"-gtf=s"        => \$gtf,
	"-coverage=s"   => \$coverage,
	"-genome=s"     => \$genome,
	"-fast_tmp=s"   => \$fast_tmp,
	"-outpath=s"    => \$outpath,
	"-options=s{,}" => \@options,
	"-sname=s"      => \$sname,
	"-n=s"          => \$n,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $R1[0] ) {
	$error .= "the cmd line switch -R1 is undefined!\n";
}elsif ( $R1[0] =~m/.txt$/ ){
	open ( IN, "<$R1[0]") or die $!;
	my @tmp;
	while ( <IN> ) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_);
		}
	}
	close ( IN );
	@R1 = @tmp;
}
unless ( defined $R2[0] ) {
	$error .= "the cmd line switch -R2 is undefined!\n";
}elsif ( $R2[0] =~m/.txt$/ ){
	open ( IN, "<$R2[0]") or die $!;
	my @tmp;
	while ( <IN> ) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_);
		}
	}
	close ( IN );
	@R2 = @tmp;
}
unless ( defined $I1[0] ) {
	$error .= "the cmd line switch -I1 is undefined!\n";
}elsif ( $I1[0] =~m/.txt$/ ){
	open ( IN, "<$I1[0]") or die $!;
	my @tmp;
	while ( <IN> ) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_);
		}
	}
	close ( IN );
	@I1 = @tmp;
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
	$warn .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $n ) {
	$n =3;
	warn "n was not defined - set 'n'umber of pcores to use to 3\n";
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
$task_description .= ' -R1 "' . join( '" "', @R1 ) . '"'  if ( defined $R1[0] );
$task_description .= ' -R2 "' . join( '" "', @R2 ) . '"'  if ( defined $R2[0] );
$task_description .= ' -I1 "' . join( '" "', @I1 ) . '"'  if ( defined $I1[0] );
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
$options->{'min_UMIs'} ||= 100;
##############################
mkdir($outpath) unless ( -d $outpath );
open( LOG, ">$outpath/" . $$ . "_10xpipeline.pl.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);
my $outpath_orig = $outpath;
$outpath .= "/10xpipeline_run_files";
mkdir($outpath) unless ( -d $outpath );

my $fm = root->filemap( "$outpath/" . $$ . "_10xpipeline.pl.log");
$outpath = $fm->{'path'};

## turn on autoflush for the process bar:
$| = 1;

## default slurm options
$options->{'A'} ||= 'lsens2017-3-2';
$options->{'n'} ||= 1;
$options->{'N'} = 1;
$options->{'t'} ||= '02:00:00';

my $SLURM = stefans_libs::SLURM->new($options, 0);
$SLURM->{'debug'} = 1 if ($debug);

## Do whatever you want!

## first we need to run the SplitToCells.pl script. this should definitely be run on a blade.
my ( $SLURM_ids, $sumFastqs ) = &SplitToCell( $SLURM );
print "waiting for SplitToCell to finish\n";
while ( ! $SLURM->pids_finished( @$SLURM_ids ) ){
	print ".";
	sleep(30); # 30 sec waiting time (to show live)
}
print  "\n";
#hisat2_run_aurora.pl # reimplementation

## sometimes this does not finish - lets check that agin:
while (scalar(@$SLURM_ids) > 0 ) {
	print "Better check once more:\n";
	( $SLURM_ids, $sumFastqs ) = &SplitToCell( $SLURM );
	while ( ! $SLURM->pids_finished( @$SLURM_ids ) ){
		print ".";
		sleep(30); # 30 sec waiting time (to show live)
	}
	print  "\n";
}



#hisat2_run_aurora.pl 
#print root::get_hashEntries_as_string($options,3, "the oSLURM ptions:");

my $hisat2 =
"hisat2_run_aurora.pl -files '".join("' '", @$sumFastqs)."' -outpath $outpath/HISAT2_mapped/ -options n 5 partitition $options->{'p'} A ".$options->{'A'};
$hisat2 .=
" -genome $genome -coverage $coverage -bigwigTracks $outpath/hisat2_$sname.html";
$hisat2 .= " -fast_tmp '$fast_tmp'";

print "we start the hist2 mapping process:\n".$hisat2."\n";

#die "First make sure we got all files!\n";

unless ( $debug) {
	open( RUN, $hisat2 . " | " );
	$SLURM_ids = [];
	foreach (<RUN>) {
		if ( $_ = /Submitted batch job (\d+)/ ) {
			push(@$SLURM_ids,  $1 );
		}
	}
	close(RUN);
}

#die "Please check, that all ".scalar(@$sumFastqs)." mapping scripts are processed\n".join("\n",@$sumFastqs)."\n";
print "waiting for the batch jobs to finish\n";

while ( ! $SLURM->pids_finished( @$SLURM_ids ) ){
	print ".";
	sleep(30); # 30 sec waiting time (to show live)
}
print  "\n";

open ( IN, "ls $outpath/HISAT2_mapped/*.bam | wc -l |" ) or die "could not run 'ls $outpath/HISAT2_mapped/*.bam | wc -l'\n";
my @OK= map { chomp; $_} <IN>;
close ( IN );

unless ( $OK[0] == scalar(@$sumFastqs)){
	Carp::confess ( "The HISAT2 scripts did not produce the expected output\n"
	."I only got $OK[0] outfiles and expected ". scalar(@$sumFastqs)." mapped bam files\n"
	);
}

## now reshuffle the bam files

## this is only one programm call - hence use all cores we have:

$SLURM->{'options'}->add('n', $n);

my $cmd = 'ReshuffleBamFiles.pl';

$cmd .= " -bams $outpath/HISAT2_mapped/*.sorted.bam";
$cmd .= " -outpath $outpath/reshuffled/";
$cmd .= " -coverage $coverage";
$cmd .= " -sampleID $sname";
$cmd .= " -n $n";
$cmd .= " -tmp_path '$fast_tmp'\n";

print "Starting reshuffling the bam files to chromosome order\n";

unless ( -d "$outpath/reshuffled/") {
	mkdir ( "$outpath/reshuffled/"  )
}

print "$cmd\n";
unless ( -f "$outpath/reshuffled/chr1_". $sname. ".sorted.bam" ){
	@$SLURM_ids = ($SLURM->run( $cmd, "$outpath/reshuffled/". $sname. "_reshuffle.sorted.bam") );
}

print "waiting for the batch job to finish\n";

while ( ! $SLURM->pids_finished( @$SLURM_ids ) ){
	print ".";
	sleep(30); # 30 sec waiting time (to show live)
}
print  "\n";

## now I need to Quantify this sample

#QuantifyBamFile.pl -infile ./HVJ5GBGX2/outs/fastq_path/Tobias1/s1/CLP_S1_L001_R2_001.fastq.gz -outfile perl_requant/CLP_S1_L001_R2 -gtf_file ~/lunarc/genomes/mouse/mm10/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf
# This takes so much memory I need to run it on the frontend ;-(

opendir ( DIR, "$outpath/reshuffled/" ) or die "I could not read from the directory '$outpath/reshuffled/'\n";
my @chrBams = map{  "$outpath/reshuffled/$_" } grep { !/^\./ }  grep { /$sname.sorted.bam/ } readdir(DIR) ;
closedir( DIR );
my $paths;
print "Starting to quantify the bam files\n";
($SLURM_ids, $paths ) = QuantifyBamFile($SLURM, @chrBams);

print "waiting for the batch jobs to finish\n";

while ( ! $SLURM->pids_finished( @$SLURM_ids ) ){
	print ".";
	sleep(30); # 30 sec waiting time (to show live)
}
print  "\n";


$cmd = "JoinResultTables.pl";
$cmd .= " -outfile $outpath_orig/$sname";
$cmd .= " -paths '".join("' '", @$paths )."'";
$cmd .= " -createTable";

if ($debug) {
	print $cmd."\n";
}
else {
	print $cmd."\n";
	unless ( -f "$outpath/$sname.sqlite" ){
		system($cmd )
	}else {
		warn ("The sqlite databse exists - please remove if I should re-create it\n"
		."rm -f $outpath/$sname.sqlite\n");
	}
}

sub QuantifyBamFile {
	my ( $SLURM, @chrBams) = @_;
	my ( @ids, @paths);
	$SLURM->{'options'}->add('n', 1);
	
	foreach my $bam_file ( @chrBams ) {
		my $fm = root->filemap( $bam_file );
		my $opath = "$fast_tmp/$sname"."_$fm->{'filename_base'}";
		if ( -d  $opath) {
			system( 'rm -RF '.$opath );
		}
		my $cmd = "QuantifyBamFile.pl";
		$cmd .= " -infile $bam_file";
		$cmd .= " -outfile $opath.sqlite";
		$cmd .= ' -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );
		$cmd .= " -fastqPath $outpath";
		$cmd .= " -sampleID $sname.$fm->{'filename_base'}";
		$cmd .= " -gtf_file $gtf";
		
		$cmd .= "\ncp -R $opath $outpath/$sname.$fm->{'filename_base'}";
		
		push ( @paths, "$outpath/$sname.$fm->{'filename_base'}" );
		unless ( -f "$outpath/$sname.$fm->{'filename_base'}/barcodes.tsv" ) {
			push( @ids, $SLURM->run($cmd, "$outpath/$fm->{'filename_base'}"."_$sname" ) )
		}
		
	}
	return ( \@ids, \@paths );
}

sub wait_for_PID {
	my $ID = shift;
	if ( $ID > 0 ) {    ## else the file has already been prepared
		while ( !$SLURM->pids_finished($ID) ) {
			sleep(10);
		}
	}
}

sub fix_path_problems {
	my @use_as_separator = @_;
	my $fm = root->filemap( $R1[0] );
	my @path = split( "/", $fm->{'path'} );
	my @path2;
	my $problems;
	my $seen;
	for (my $i = 0; $i < @R1; $i ++ ){
		my $fm = root->filemap( $R1[$i] );
		@path2 = split( "/", $fm->{'path'} );
		my $key;
		if ( scalar(@use_as_separator) > 0 ) {
			$key = join("_", @path2[@use_as_separator], $fm->{'filename'} );
		}else {
			$key = $fm->{'filename'};
		}
		if ( $seen->{$key} ) {
		for ( my $a = 0; $a <  @path2; $a ++) {
			unless ( $path2[$a] eq $path[$a] ) {
				$problems->{$a}->{'orig'} = $path[$a];
				$problems->{$a}->{$path2[$a]} ++;
			}
			
		}
		}
		$seen->{$key} = $R1[$i];
	}
	$problems->{33} = { 'fix' => $seen};
	
	
	my @tmp = sort {$a <=> $b} keys %$problems ;
	my @fix = shift (@tmp);
	my $i = 0;
	#warn "this is not working: ".scalar(keys %$problems )."\n";
	if ( scalar(keys %$problems ) > 1 ){
		#print "problematic path parts:\n \$problems = ".root->print_perl_var_def( $problems ).";\n";
		## So now I would suppose I should try to start in the beginning to separate the files.
		$problems = &fix_path_problems( @fix );
		while ( scalar(keys %$problems ) > 1 ) {
			#print "problematic path parts with adding path [".join(', ',@fix)."] to the filenames:\n \$problems = ".root->print_perl_var_def( $problems ).";\n";
			@tmp = sort {$a <=> $b} keys %$problems ;
			push (@fix ,shift (@tmp));
			$problems = &fix_path_problems( @fix );
			if ($i ++  == 20 ) {
				#warn "Path problems cound not be fixed!\n";
				#last();
			}
		}
	}
	return $problems;
}

sub SplitToCell {
	my ( $SLURM ) = @_;
	
	## now I need to start one script for every R1 R2 I1 combination
	my (@slurmIDs, @fastqs);

	# Here I need to make sure, that the path the files are in are all the same
	
	my $rev = &fix_path_problems ()->{33}->{'fix'};
	my $fix = { map{ $rev->{$_} => $_ } keys %$rev};
	#die "this should fix the path problems::\n \$problems = ".root->print_perl_var_def( $fix ).";\n";
		
	for (my $i = 0; $i < @R1; $i ++ ){
		my $ofile = $fix ->{ $R1[$i] };
		my $f = root->filemap( $R1[$i] );
		
		$ofile =~s/.R1.*$//;
		
		my $cmd = 'SplitToCells.pl';
		$cmd .= " -I1 $I1[$i]";
		$cmd .= " -R1 $R1[$i]";
		$cmd .= " -R2 $R2[$i]";
		$cmd .= " -outpath $fast_tmp";
		$cmd .= " -options oname $ofile";
		$cmd .= "\ncp $fast_tmp/$ofile.annotated.fastq.gz $outpath\n";
		$cmd .= "\ncp $fast_tmp/$ofile.per_cell_read_count.xls $outpath\n";
		$cmd .= "cp $fast_tmp/$ofile*_SplitToCells.pl.log $outpath\n";
		unless ( -f  "$outpath/$ofile.annotated.fastq.gz" ){
			push( @slurmIDs, $SLURM->run( $cmd, "$outpath/$ofile.annotated.fastq.gz" ) );
		}
		push( @fastqs , "$outpath/$ofile.annotated.fastq.gz" );
		
	}
	return ( \@slurmIDs, \@fastqs);
}

