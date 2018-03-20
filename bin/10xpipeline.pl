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
       
 -local (option) :Do not use the SLURM blades, but calculate all on the frontend
                  or if you are working on a workstation without SLURM environment.
       
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
use Parallel::ForkManager;
use stefans_libs::file_readers::gtf_file;

use DateTime;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help,     $debug,  $database, $gtf,     $fast_tmp, $local,
	$coverage, $genome, $outpath,  $options, @options,  @R1,
	@R2,       @I1,     $sname,    $n,       $avail_chr
);

Getopt::Long::GetOptions(
	"-R1=s{,}"      => \@R1,
	"-R2=s{,}"      => \@R2,
	"-I1=s{,}"      => \@I1,
	"-gtf=s"        => \$gtf,
	"-coverage=s"   => \$coverage,
	"-genome=s"     => \$genome,
	"-fast_tmp=s"   => \$fast_tmp,
	"-outpath=s"    => \$outpath,
	"-options=s{,}" => \@options,
	"-sname=s"      => \$sname,
	"-n=s"          => \$n,
	"-local"        => \$local,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $start = DateTime->now();
my ( $end, $start_this );
print "Started at " . $start->time() . "\n";

my $warn  = '';
my $error = '';

## one other option: You only give me R1 data and I figure out where to put the files myself:

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/10xpipeline.pl';
$task_description .= ' -R1 "' . join( '" "', @R1 ) . '"' if ( defined $R1[0] );
$task_description .= ' -R2 "' . join( '" "', @R2 ) . '"' if ( defined $R2[0] );
$task_description .= ' -I1 "' . join( '" "', @I1 ) . '"' if ( defined $I1[0] );
$task_description .= " -gtf '$gtf'"           if ( defined $gtf );
$task_description .= " -coverage '$coverage'" if ( defined $coverage );
$task_description .= " -genome '$genome'"     if ( defined $genome );
$task_description .= " -outpath '$outpath'"   if ( defined $outpath );
$task_description .= ' -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );
$task_description .= " -sname $sname";
$task_description .= ' -debug' if ($debug);
$task_description .= ' -local' if ($local);

if ( -d $R1[0] ) {
	## OK that is the best way to do it :-D
	my $fm = root->filemap( $R1[0] . "/test.txtx" );
	@R2 = undef;
	@I1 = undef;
	my $path = $fm->{'path'};    ## abs path please
	open( IN, "find $path -name '*$sname*.gz' |" )
	  or die "I could not start the find subprocess\n$!\n";
	@R1 = map { chomp; split( /\s+/, $_ ) } <IN>;

	close(IN);
}
if ( -f $R1[0] and !defined $R2[0] ) {    ## I just assume I1 is also empty..
	my @tmp = @R1;
	@R2 = grep { /_R2_/ } @tmp;
	@I1 = grep { /_I1_/ } @tmp;
	@R1 = grep { /_R1_/ } @tmp;

#die "\$exp = ".root->print_perl_var_def( {'R1' => [@R1[1,2]], 'R2' => [@R2[1,2]], 'I1' => [@I1[1,2]] } ).";\n";

}

unless ( defined $R1[0] ) {
	$error .= "the cmd line switch -R1 is undefined!\n";
}
elsif ( $R1[0] =~ m/.txt$/ ) {
	open( IN, "<$R1[0]" ) or die $!;
	my @tmp;
	while (<IN>) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_ );
		}
	}
	close(IN);
	@R1 = @tmp;
}
unless ( defined $R2[0] ) {
	$error .= "the cmd line switch -R2 is undefined!\n";
}
elsif ( $R2[0] =~ m/.txt$/ ) {
	open( IN, "<$R2[0]" ) or die $!;
	my @tmp;
	while (<IN>) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_ );
		}
	}
	close(IN);
	@R2 = @tmp;
}
unless ( defined $I1[0] ) {
	$error .= "the cmd line switch -I1 is undefined!\n";
}
elsif ( $I1[0] =~ m/.txt$/ ) {
	open( IN, "<$I1[0]" ) or die $!;
	my @tmp;
	while (<IN>) {
		chomp;
		if ( -f $_ ) {
			push( @tmp, $_ );
		}
	}
	close(IN);
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
	$n = 3;
	warn "n was not defined - set 'n'umber of pcores to use to 3\n";
}
unless ($fast_tmp) {
	$fast_tmp = '$SNIC_TMP';
}
unless ($sname) {
	$error .= "The sname was set to 'sampleX'\n";
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

my $fm = root->filemap( "$outpath/" . $$ . "_10xpipeline.pl.log" );

#$outpath = $fm->{'path'};

## turn on autoflush for the process bar:
#$| = 1;

## default slurm options
$options->{'A'} ||= 'lsens2017-3-2';
$options->{'n'} ||= 1;
$options->{'N'} = 1;
$options->{'t'} ||= '02:00:00';

my $SLURM       = stefans_libs::SLURM->new( $options, 0 );
my $SLURM_local = stefans_libs::SLURM->new( $options, 0 );
$SLURM_local->{'run_local'} = 1;

my $slurmOptions;
while ( my ( $key, $value ) = each %{ $SLURM->{options}->options() } ) {
	$slurmOptions .= " '$key' '$value'";
}

$SLURM->{'debug'} = 1 if ($debug);

my $pm;
if ($local) {
	$pm = Parallel::ForkManager->new($n);
	$SLURM->{'run_local'} = 1;
}

## Do whatever you want!

$start_this = &check_time_since( $start, 'Setup' );

#die "after setup death to test input options ".root->print_perl_var_def( {'R1' => [@R1[1,2]], 'R2' => [@R2[1,2]], 'I1' => [@I1[1,2]] } ).";\n";

## first we need to run the SplitToCells.pl script. this should definitely be run on a blade.
my ( $SLURM_ids, $sumFastqs );
unless ($local) {
	( $SLURM_ids, $sumFastqs ) = &SplitToCell($SLURM);
	print "waiting for SplitToCell to finish\n";
	while ( !$SLURM->pids_finished(@$SLURM_ids) ) {
		print ".";
		sleep(30);    # 30 sec waiting time (to show live)
	}
	print "\n";

}
else {
	( $SLURM_ids, $sumFastqs ) = &SplitToCell_local($SLURM);
}

$start_this = &check_time_since( $start_this,
	'merge R1, R2 and I1 reads + polyA and low quality filter', $SLURM_ids );

#hisat2_run_aurora.pl
#print root::get_hashEntries_as_string($options,3, "the oSLURM ptions:");

my $hisat2 =
    "hisat2_run_aurora.pl -files '"
  . join( "' '", @$sumFastqs )
  . "' -outpath $outpath/HISAT2_mapped/"
  . " -mapper_options ' --score-min L,-0.0,-0.4'" ## relax the mapper efficiency from L,0.0,-0.2 # experimentall checked against cellranger results.
  . " -options n 5 partitition $options->{'p'} A "
  . $options->{'A'};
my $used = { map { $_ => 1 } qw( n A p N ) };
while ( my ( $key, $value ) = each %{ $SLURM->{options}->options() } ) {
	$hisat2 .= " '$key' '$value'" unless ( $used->{$key} );
}
$hisat2 .=
" -genome $genome -coverage $coverage -bigwigTracks $outpath/hisat2_$sname.html";
$hisat2 .= " -fast_tmp '$fast_tmp' -justMapping";
$hisat2 .= " -debug" if ($debug);
$hisat2 .= " -local" if ($local);

print "we start the hist2 mapping process:\n" . $hisat2 . "\n";

#die "First make sure we got all files!\n";

system("rm $outpath/hisat2_run.local*");

$SLURM_local->run( $hisat2, "$outpath/hisat2_run" );
## the SLURM_local should run this script on the frontend; remove all old log files and create a $outpath/hisat2_run<PID>.err and $outpath/hisat2_run<PID>.out file

## Identify the run ids we need to wait for:

my @tmp = get_files_from_path( $outpath, "hisat2_run.local.*out" );
unless ( -f $tmp[0] ) {
	Carp::confess("I did not get a usable hisat out file: '$tmp[0]'\n");
}
open( RUN, $tmp[0] )
  or die "I could not oopen the hisat2_runaurora stdout file '$tmp[0]'\n$!\n";

$SLURM_ids = [];
foreach (<RUN>) {
	if ( $_ = /Submitted batch job (\d+)/ ) {
		push( @$SLURM_ids, $1 );
	}
}
close(RUN);

if ($debug) {
	mkdir("$outpath/HISAT2_mapped/") unless ( -d "$outpath/HISAT2_mapped/" );
	## and we would expect one outfile for all fastq files - so lets add some fake ones:
	my $fm;
	foreach my $f (@$sumFastqs) {
		$fm = root->filemap($f);
		my $nfile =
		    "$outpath/HISAT2_mapped/FAKE_DEBUG_"
		  . $fm->{'filename_base'}
		  . ".sorted.bam";
		unless ( -f $nfile ) {
			system("touch $nfile");
		}
	}
}

print "We got the HISAT2 slurm IDs: " . join( " ", @$SLURM_ids ) . "\n";

$start_this = &check_time_since( $start_this, 'HISAT2 mapping', $SLURM_ids );

## check the hisat2 results
my @OK = get_files_from_path( "$outpath/HISAT2_mapped/", ".bam\$" );

unless ( scalar(@OK) == scalar(@$sumFastqs) ) {
	Carp::confess( "The HISAT2 scripts did not produce the expected output\n"
		  . "I only got "
		  . scalar(@OK)
		  . " outfiles and expected "
		  . scalar(@$sumFastqs)
		  . " mapped bam files\n" );
}
else {
	print "HISAT2 did produce the required " . scalar(@OK) . " outfiles\n";
}

## now reshuffle the bam files

## this is only one programm call - hence use all cores we have:

$SLURM->{'options'}->add( 'n', $n );

my $cmd = 'ReshuffleBamFiles.pl';

$cmd .= " -bams $outpath/HISAT2_mapped/*.sorted.bam";
$cmd .= " -outpath $outpath/reshuffled/";
$cmd .= " -coverage $coverage";
$cmd .= " -gtf $gtf";
$cmd .= " -sampleID $sname";
$cmd .= " -n $n";
$cmd .= " -tmp_path '$fast_tmp'\n";
$cmd .= " -debug" if ($debug);

print "Starting reshuffling the bam files to chromosome order\n";

unless ( -d "$outpath/reshuffled/" ) {
	mkdir("$outpath/reshuffled/");
}

print "$cmd\n";

my @chrBams = &get_files_from_path( "$outpath/reshuffled/", ".*bam\$" );
if ( scalar(@chrBams) == 0 ) {
	@$SLURM_ids = (
		$SLURM->run(
			$cmd,
			"$outpath/reshuffled/script_run",
			"$outpath/reshuffled/chr1_" . $sname . ".sorted.bam"
		)
	  )
	  ; ## the other one will not run due to the 'debug' setting, but I created fake data and want to test this here
	if ($debug) {

		$SLURM_local->run(
			$cmd,
			"$outpath/reshuffled/$sname" . "_reshuffle_script_run",
			"$outpath/reshuffled/chr1_" . $sname . ".sorted.bam"
		);

	}
}
else {
	print
"If you lack some bam file in the output please remove all bam files in the folder  $outpath/reshuffled/\n";
}
## wait and report time
warn "I wait for the slurm ids: " . join( ", ", @$SLURM_ids ) . "\n";
$start_this =
  &check_time_since( $start_this, 'Reshuffle Bam Files', $SLURM_ids );

## now I need to Quantify this sample

#QuantifyBamFile.pl -infile ./HVJ5GBGX2/outs/fastq_path/Tobias1/s1/CLP_S1_L001_R2_001.fastq.gz -outfile perl_requant/CLP_S1_L001_R2 -gtf_file ~/lunarc/genomes/mouse/mm10/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf
# This takes so much memory I need to run it on the frontend ;-(

my @chrBams = &get_files_from_path( "$outpath/reshuffled/", ".*bam\$" );

my $paths;
print "Starting to quantify the bam files\n";
if ( $local or $debug ) {
	( $SLURM_ids, $paths ) = QuantifyBamFile_local( $SLURM, @chrBams );
}
else {
	( $SLURM_ids, $paths ) = QuantifyBamFile( $SLURM, @chrBams );
}

## wait and report time
warn "I wait for the slurm ids: " . join( ", ", @$SLURM_ids ) . "\n";

$start_this = &check_time_since( $start_this, 'quantify', $SLURM_ids );

$cmd = "JoinResultTables.pl";
$cmd .= " -outfile $outpath_orig/$sname";
$cmd .= " -paths '" . join( "' '", @$paths ) . "'";
$cmd .= " -createTable";
$cmd .= " -force";               ## I want the database to be re-created!
$cmd .= " -debug" if ($debug);

print "JoinResultTables is started like that:\n" . $cmd . "\n";

if ( -f "$outpath/$sname.sqlite" ) {
	warn( "The sqlite databse exists - please remove if I should re-create it\n"
		  . "rm -f $outpath/$sname.sqlite\n" );
}

$SLURM_local->run( $cmd, "$outpath/JoinResultTables_$sname",
	"$outpath/$sname.sqlite" );

warn "please delete all files in to outpath $outpath after inspection!\n"
  if ($debug);

$start_this = &check_time_since( $start_this, 'JoinResultTables' );

$end = DateTime->now();
print "\n10xpipeline.pl Finished at" . $end->time() . "\n";
print "n10xpipeline.pl total run time: "
  . join( ":",
	$end->subtract_datetime($start)->in_units( 'days', 'hours', 'seconds' ) )
  . "\n";

sub get_files_from_path {
	my ( $path, @matches ) = @_;
	opendir( DIR, $path )
	  or die "I could not read from the directory '$path'\n$!\n";
	my @dat = grep { !/^\./ } readdir(DIR);
	closedir(DIR);
	print "I get files from path '$path'\n";
	foreach my $select (@matches) {

		print "I select all files matching '$select'\n"
		  . join( "\n", @dat ) . "\n";
		@dat = grep { /$select/ } @dat;

		print "Still in the game:" . join( "\n", @dat ) . "\n\n";
	}
	@dat = map { "$path/$_" } @dat;
	return @dat;
}

=head3 check_time_since ( $start_last, $msg )

check_time_since ( <DateTime::object>, <string> )
check_time_since ( <DateTime::object 10 min ago>, "Setup" ) would print e.g

Now: 07:09:15 Setup - processing time: 0:0:2

and return the DateTime->now() object used in the call.

=cut

sub check_time_since {
	my ( $start_last, $msg, $SLURM_ids ) = @_;

	if ( ref($SLURM_ids) eq "ARRAY" and scalar(@$SLURM_ids) > 0 ) {
		print "waiting for the batch job(s) to finish\n";

		while ( !$SLURM->pids_finished(@$SLURM_ids) ) {
			print ".";
			sleep(30);    # 30 sec waiting time (to show live)
		}
		print "\n";
	}

	$end = DateTime->now();
	my $str = "# "
	  . $end->time()
	  . " - $msg - processing time: "
	  . join( ":",
		$end->subtract_datetime($start)->in_units( 'days', 'hours', 'seconds' )
	  ) . "\n";
	print $str;
	return $end;
}

sub QuantifyBamFile_CMD {
	my $bam_file = shift;
	my $fm       = root->filemap($bam_file);
	my $chr_slice;
	if ( $bam_file =~ m/.(chr[\w\d]+:\d+-\d+)/ ) {
		$chr_slice = $1;
	}
	else {
		$chr_slice = 'chrunknown_slice';
	}
	my @tmp = split( /chr/, $chr_slice );
	$chr_slice = 'chr' . pop(@tmp);
	my $opath = "$fast_tmp/$chr_slice" . "_$sname";
	$opath =~ s/\./_/g;
	$opath .= "_$fm->{'filename_base'}";

	if ( -d $opath ) {
		system( 'rm -RF ' . $opath );
	}
	my $cmd = "QuantifyBamFile.pl";

	$cmd .= " -infile $bam_file";
	$cmd .= " -outfile $opath.sqlite";
	$cmd .= ' -options ' . $slurmOptions
	  if ( defined $options[0] );
	$cmd .= " -fastqPath $outpath";
	$cmd .= " -sampleID $sname.$fm->{'filename_base'}_$chr_slice";
	$cmd .= " -gtf_file $gtf";

	my $tmp = $fm->{'filename'};
	if ( $tmp =~ m/(.*)_$sname/ ) {

		#die "I would add the chromosome ".&getChr_name($1)."\n";
		$cmd .= " -drop_chr " . &getChr_name($1);
	}
	else {
		die
"Hey we should be able to identify a chromosome in $tmp - but we failed!\n";
	}

	$cmd .= " -debug" if ($debug);

	$cmd .= "\nmv $opath $outpath/$chr_slice" . "_$sname";
	$cmd .= "\nmv $opath/$chr_slice*.log $outpath/$chr_slice"
	  . "_$sname/";    #copy the log file, too.

	return ( $cmd, "$outpath/$chr_slice" . "_$sname/" );
}

sub byFileSize {
	-s $b <=> -s $a;
}


sub FileSizeOrder {
	my @files = @_;
	my $i = 0;
	my $order = { map { $_ => $i ++ }  @files  } ; 
	my @ret;
	foreach ( sort byFileSize @files ) {
		push( @ret, $order->{$_});
	}
	return @ret;
}

sub getChr_name {
	my $chr = shift;
	unless ( defined $avail_chr ) {
		open( IN, "<$coverage" )
		  or die
		  "getChr_name - I can not read from coverage file '$coverage'\n$!\n";
		my @line;
		while (<IN>) {
			chomp;
			@line = split( /\s+/, $_ );
			$avail_chr->{ $line[0] } = 1;
		}
		close(IN);
	}
	if ( $avail_chr->{$chr} ) {
		return $chr;
	}
	$chr =~ s/^chr//;
	if ( $avail_chr->{$chr} ) {
		return $chr;
	}
	## now it starts to get tricky - there might be a .numer in the end that I loose in the process.
	## there is only a .1 at the moment - so lets hope that stays like that and go for it.
	if ( $avail_chr->{ $chr . ".1" } ) {
		return $chr . ".1";
	}
	## so now I am not sure what to do - lets die
	die "Sorry, but the chr $chr is not in my list of chromosomes: \n\t"
	  . join( "\n\t", sort keys %$avail_chr ) . "\n";
}

sub QuantifyBamFile_local {
	my ( $SLURM, @chrBams ) = @_;
	my ( @ids, @paths );
	$SLURM->{'options'}->add( 'n', 1 );

	print "using this comuter to calculate - this can take a lot of time!\n";
  FILES:
	foreach my $bam_file ( sort byFileSize @chrBams ) {
		my $fm = root->filemap($bam_file);
		my ( $cmd, $path ) = QuantifyBamFile_CMD($bam_file);
		push( @paths, $path );
		my $pid = $pm->start and next FILES;
		print "producing files in $path\n";
	#	print "Quantify cmd:" . $cmd . "\n";
		$SLURM_local->run( $cmd,
			"$outpath/QuantifyBamFile_$fm->{'filename'}" . "_$sname",
			"$path/barcodes.tsv" );
		$pm->finish;
	}
	$pm->wait_all_children;
	return ( \@ids, \@paths );
}

sub QuantifyBamFile {
	my ( $SLURM, @chrBams ) = @_;
	my ( @ids, @paths );
	$SLURM->{'options'}->add( 'n', 1 );

	foreach my $bam_file ( sort byFileSize @chrBams ) {

		my $fm = root->filemap($bam_file);

		my ( $cmd, $path ) = QuantifyBamFile_CMD($bam_file);
		push( @paths, $path );

		push(
			@ids,
			$SLURM->run(
				$cmd, "$outpath/QuantifyBamFile_$fm->{'filename'}" . "_$sname",
				"$path/barcodes.tsv"
			)
		);

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
	my $fm               = root->filemap( $R1[0] );
	my @path             = split( "/", $fm->{'path'} );
	my @path2;
	my $problems;
	my $seen;
	for ( my $i = 0 ; $i < @R1 ; $i++ ) {
		my $fm = root->filemap( $R1[$i] );
		@path2 = split( "/", $fm->{'path'} );
		my $key;
		if ( scalar(@use_as_separator) > 0 ) {
			$key = join( "_", @path2[@use_as_separator], $fm->{'filename'} );
		}
		else {
			$key = $fm->{'filename'};
		}
		if ( $seen->{$key} ) {
			for ( my $a = 0 ; $a < @path2 ; $a++ ) {
				next
				  unless ( defined $path2[$a] )
				  ;    # the first one is alwas undefined...
				 #warn "Which one is undefined? 2 ($a)? '$path2[$a]' or 1? '$path[$a]'?\n";
				unless ( $path2[$a] eq $path[$a] ) {
					$problems->{$a}->{'orig'} = $path[$a];
					$problems->{$a}->{ $path2[$a] }++;
				}

			}
		}
		$seen->{$key} = $R1[$i];
	}
	$problems->{33} = { 'fix' => $seen };

	my @tmp = sort { $a <=> $b } keys %$problems;
	my @fix = shift(@tmp);
	my $i   = 0;

	#warn "this is not working: ".scalar(keys %$problems )."\n";
	if ( scalar( keys %$problems ) > 1 ) {

#print "problematic path parts:\n \$problems = ".root->print_perl_var_def( $problems ).";\n";
		## So now I would suppose I should try to start in the beginning to separate the files.
		$problems = &fix_path_problems(@fix);
		while ( scalar( keys %$problems ) > 1 ) {

#print "problematic path parts with adding path [".join(', ',@fix)."] to the filenames:\n \$problems = ".root->print_perl_var_def( $problems ).";\n";
			@tmp = sort { $a <=> $b } keys %$problems;
			push( @fix, shift(@tmp) );
			$problems = &fix_path_problems(@fix);
			if ( $i++ == 20 ) {

				#warn "Path problems cound not be fixed!\n";
				#last();
			}
		}
	}
	return $problems;
}

sub SpitCell_CMD {
	my ( $R1, $R2, $I1, $fast_tmp, $ofile ) = @_;
	my $cmd = 'SplitToCells.pl';
	$cmd .= " -I1 $I1";
	$cmd .= " -R1 $R1";
	$cmd .= " -R2 $R2";
	$cmd .= " -outpath $fast_tmp";
	$cmd .= " -options oname $ofile" . $slurmOptions;
	$cmd .= "\nmv $fast_tmp/$ofile.annotated.fastq.gz $outpath\n";
	$cmd .= "\nmv $fast_tmp/$ofile.per_cell_read_count.xls $outpath\n";
	$cmd .= "mv $fast_tmp/$ofile*_SplitToCells.pl.log $outpath\n";
	return $cmd;
}


sub SplitToCell_local {
	my ($SLURM) = @_;

	## now I need to start one script for every R1 R2 I1 combination
	my ( @slurmIDs, @fastqs );

	# Here I need to make sure, that the path the files are in are all the same

	my $rev = &fix_path_problems()->{33}->{'fix'};
	my $fix = { map { $rev->{$_} => $_ } keys %$rev };
	print "using this comuter to calculate - this can take a lot of time!\n";

#die "this should fix the path problems::\n \$problems = ".root->print_perl_var_def( $fix ).";\n";
  FILES:
	foreach my $i ( &FileSizeOrder(@R2) ) { ## start with the biggest
		my $pid   = $pm->start and next FILES;
		my $ofile = $fix->{ $R1[$i] };
		my $f     = root->filemap( $R1[$i] );
		my $start = DateTime->now();
		$ofile =~ s/.R1.*$//;
		print
"processing files R1 R2 and I1 like $f->{'filename'} in path $f->{'path'}\n\tUsing script $outpath/$f->{'filename'}_SplitToCell.sh\n";

		my $cmd = &SpitCell_CMD( $R1[$i], $R2[$i], $I1[$i], $fast_tmp, $ofile );
		if ($debug) {
			system( 'touch ' . "$outpath/$ofile.annotated.fastq.gz" );
		}
		else {
			unless ( -f "$outpath/$ofile.annotated.fastq.gz" ) {

				warn
"the outfile '$outpath/$ofile.annotated.fastq.gz' is missing\n";
				push(
					@slurmIDs,
					$SLURM_local->run(
						$cmd,
						"$outpath/$ofile.SplitToCell",
						"$outpath/$ofile.annotated.fastq.gz"
					)
				);
			}
		}
		push( @fastqs, "$outpath/$ofile.annotated.fastq.gz" );

		$start_this =
		  &check_time_since( $start, 'SplitToCell ' . $f->{'filename'},
			\@slurmIDs );

		$pm->finish;
	}
	$pm->wait_all_children;

	opendir( DIR, $outpath );
	@fastqs =
	  map { "$outpath/$_" } grep { /.annotated.fastq.gz$/ } readdir(DIR);
	closedir(DIR);

	return ( [], \@fastqs );
}

sub SplitToCell {
	my ($SLURM) = @_;

	## now I need to start one script for every R1 R2 I1 combination
	my ( @slurmIDs, @fastqs );
	$SLURM->options( 'n', 2 )
	  ; ## set the number of cores to two (gzip in gzip out and the main reformate)
	my $rev = &fix_path_problems()->{33}->{'fix'};
	my $fix = { map { $rev->{$_} => $_ } keys %$rev };

	foreach my $i ( &FileSizeOrder(@R2) ) { ## start with the biggest
		my $ofile = $fix->{ $R1[$i] };
		my $f     = root->filemap( $R1[$i] );

		$ofile =~ s/.R1.*$//;

		my $cmd = &SpitCell_CMD( $R1[$i], $R2[$i], $I1[$i], $fast_tmp, $ofile );
		if ( !-f "$outpath/$ofile.annotated.fastq.gz" ) {
			push(
				@slurmIDs,
				$SLURM->run(
					$cmd,
					"$outpath/$ofile.SplitToCell",
					"$outpath/$ofile.annotated.fastq.gz"
				)
			);
		}
		push( @fastqs, "$outpath/$ofile.annotated.fastq.gz" );
		if ($debug) {
			system( 'touch ' . "$outpath/$ofile.annotated.fastq.gz" );
		}

	}
	return ( \@slurmIDs, \@fastqs );
}

