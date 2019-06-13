#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-02-21 Stefan Lang

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

=head1  SYNOPSIS

    QuantifyBamFile.pl
       -infile   :one bam file containing thousands of different cells (Illumnina 10x Chromium)
       -gtf_file :the gtf file for quantification
       -outfile  :the tab separated outfile (will be overwritten)
       -drop_chr :state here the chromosome you want to read from the gtf file
                  read all if undefined
       
       -fastqPath :The path where the original fastq files were converted to the here used input files
       -sampleID  :The sampleID of the 10X sample to process
       -is10x     :10x data reads aligne on the same strand as the RNA is from.
                   this is used to identify matching genes - set to 0 is you are not using 10x data here!
                   
      -asDatabase :Export the data as database, not tables (default not used)
      
       -options  : format: key_1 value_1 key_2 value_2 ... key_n value_n

           quantify_on: which feature to select for quantification from the gtf file (default 'exon')
             report_on: which column in the gtf file object to report on (default 'gene_id')
              min_UMIs: how many UMIs have to mapp to any gene to report the sample 
                        (for one NextSeq run the default value default 100 proved usable)

       -help           :print this help
       -debug          :verbose output
       -bugfix         :extremely verbose output all data created
   
=head1 DESCRIPTION

  Take a Chromium 10x bam file and quantify it agains a gtf file.

  To get further help use 'QuantifyBamFile.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

#use POSIX;

use stefans_libs::BAMfile;
use stefans_libs::result_table;

#use stefans_libs::file_readers::gtf_file;
use Digest::MD5 qw(md5_hex);

use stefans_libs::database::Chromium_SingleCell::datavalues;
use stefans_libs::GeneModelMatcher;

use PDL;
use PDL::NiceSlice;

use Parallel::ForkManager;

$PDL::BIGPDL = 1;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

#Lpo start 87806428
#Mpo end   87804413
my (
	$help,    $debug,   $database,  $infile,   $gtf_file, $drop_chr, $is10x,
	$outfile, $options, $fastqPath, $sampleID, @options,  $asDatabase, $bugfix
);

Getopt::Long::GetOptions(
	"-infile=s"     => \$infile,
	"-gtf_file=s"   => \$gtf_file,
	"-outfile=s"    => \$outfile,
	"-options=s{,}" => \@options,
	"-drop_chr=s"   => \$drop_chr,
	"-fastqPath=s"  => \$fastqPath,
	"-sampleID=s"   => \$sampleID,
	"-asDatabase"   => \$asDatabase,
	"-bugfix"       => \$bugfix,
	"-is10x=s"		=> \$is10x,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $infile ) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $gtf_file ) {
	$error .= "the cmd line switch -gtf_file is undefined!\n";
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ( defined $options[0] ) {
	$warn .= "the cmd line switch -options is undefined!\n";
}

unless ( -d $fastqPath ) {
	$error .= "the cmd line switch -fastqPath is undefined!\n";
}
unless ( defined $sampleID ) {
	$error .= "the cmd line switch -sampleID is undefined!\n";
}
unless ( defined $is10x) {
	$warn .= "assuming 10x data\n";
	$is10x = 1;
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

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/QuantifyBamFile.pl';
$task_description .= " -infile '$infile'" if ( defined $infile );
$task_description .= " -gtf_file '$gtf_file'" if ( defined $gtf_file );
$task_description .= " -outfile '$outfile'" if ( defined $outfile );
$task_description .= ' -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );
$task_description .= ' -fastqPath ' . $fastqPath;
$task_description .= ' -drop_chr ' . $drop_chr if ( defined $drop_chr );
$task_description .= ' -sampleID ' . $sampleID;
$task_description .= ' -asDatabase' if ($asDatabase);
$task_description .= ' -bugfix' if ($bugfix);
$task_description .= ' -is10x $is10x' if (defined $is10x );

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}

### initialize default options:

$options->{'quantify_over'} ||= 'gene';
$options->{'quantify_on'} ||= 'exon';
$options->{'report_on'}   ||= 'gene_id';
unless ( defined $options->{'min_UMIs'} ){
	$options->{'min_UMIs'} = 100; ## could be set to 0
}


###

my $self = {
	'chr'         => '',
	'end'         => 0,
	'UMI'         => {},
	'mergable'    => [],
	'UMIs'        => {},
	'duplicates'  => 0,
	'total_reads' => 0,
};    # the script should store global variables here

## turn on autoflush for the process bar:
#$| = 1; ## rather not do that now that I have automated the scripts :-)

my $fm = root->filemap($outfile);
unless ( -d $fm->{'path'} ){
	mkdir( $fm->{'path'} )
}


open( LOG, ">$outfile.log" ) or die "I could not create the log file '$outfile.log'\n$!\n";
print LOG $task_description . "\n";

print "I am processing file '$infile'\n";
## Do whatever you want!

my $bam_file = stefans_libs::BAMfile->new();
if ($debug) {
	$outfile .= "_FAKE_DEBUG.sqlite"
	  unless ( $outfile =~ m/\_FAKE_DEBUG.sqlite$/ );
}
else {
	$outfile .= ".sqlite" unless ( $outfile =~ m/\.sqlite$/ );
}

## the real outfile here is 
my $tmp = $outfile;
$tmp =~ s!.sqlite$!/matrix.mtx!;
if ( -f $tmp ) {
	warn "I am going to add data to the outfile(s) '$tmp'\nremove this file if you want to recreate it.\n";
}

my $result = stefans_libs::result_table->new(
	{
		filename               => $outfile,
		'data_storage_spliced' => 0,

		#	'table_path'           => $drop_chr # not use - might crap things up
	}
);

#$fm = root->filemap($outfile);
#
#my $spliced = stefans_libs::result_table->new(
#	{
#		filename               => $fm->{'path'} ."/", $fm->{'filename_base'}."_spliced.sqlite",
#		'data_storage_spliced' => 0,
#		#	'table_path'           => $drop_chr # not use - might crap things up
#	}
#);
#
#my $unspliced = stefans_libs::result_table->new(
#	{
#		filename               => $fm->{'path'} ."/", $fm->{'filename_base'}."_unspliced.sqlite",
#		'data_storage_spliced' => 0,
#
#		#	'table_path'           => $drop_chr # not use - might crap things up
#	}
#);

my $start       = time;
my $first_start = $start;

my $GeneModelMatcher = stefans_libs::GeneModelMatcher->new( 
	{ 
		'collect_over' => $options->{'quantify_over'} , 
		'collect' => $options->{'quantify_on'},
		'id' => $options->{'report_on'}, 
		'is10x' => $is10x,
	} );

print "I read the chr $drop_chr from the gtf file\n" if ( defined $drop_chr );

$GeneModelMatcher->read_file($gtf_file, $drop_chr);

$GeneModelMatcher->{'debug'} = $debug;

$GeneModelMatcher -> {'gtf_file'}->GeneFreeSplits(); ## define the gene free splits in order to keep problems with missing genes as low as possible.

print "Available gene models:\n\t".join("\n\t", $GeneModelMatcher->Info() ).";\n" if ( $debug);

my (
	@matching_IDs, @matching_IDs2, $N,          $Seq,
	$sample_name,  $sample_table,  $sample_row, $sum,
	$max,          $min,           $rem
);

my $duration = time - $start;
print "Preparation of the gtf file took: $duration s\n";
print LOG "Preparation of the gtf file took: $duration s\n";
$start = time;

print
  "processing the initial sample summary files and detect usable cell ids.\n";
my ( $OK, @tmp, $reads );

if ( -f $fastqPath . "/.$sampleID.passing_samples.$options->{'min_UMIs'}.txt" ) {
	print "using hidden file $fastqPath/.$sampleID.passing_samples.$options->{'min_UMIs'}.txt\n";
	open( IN, "<" . $fastqPath . "/.$sampleID.passing_samples.$options->{'min_UMIs'}.txt" );
	while (<IN>) {
		chomp;
		@tmp = split( "\t", $_ );
		$OK->{ $tmp[0] } = $tmp[1];
	}
	close(IN);
}
else {
	opendir( DIR, $fastqPath )
	  or die "I could not open the path '$fastqPath'\n$!";

	my @samplefiles = grep { !/^\./ } grep { $sampleID } readdir(DIR);
	@samplefiles = map { join( "/", $fastqPath, $_ ) } @samplefiles;
	if ( @samplefiles == 0 ){
		
	}
	foreach my $file ( grep { /xls$/ } @samplefiles ) {
		print "I process file $file\n";
		open( IN, "<$file" ) or die "I could not open file '$file'\n$!\n";
		while (<IN>) {
			next if $_ =~ m/^#/;
			chomp;
			@tmp = split( "\t", $_ );
			next if ( $tmp[1] eq "count" );
			if ( $tmp[1] >= $options->{'min_UMIs'} ) {
				if ( $tmp[0] =~m/_/ ){
					$reads = $tmp[1];
					@tmp = split( "_", $tmp[0] );
					$OK->{ $tmp[3] } += $reads;
				}else {
					$OK->{ $tmp[0] } += $tmp[1];
				}
			}
		}
		close(IN);
	}
	if ( !-f . $fastqPath . "/.$sampleID.passing_samples.$options->{'min_UMIs'}.txt" ) {
		## an other process might have been faster ;-)
		open( OUT, ">" . $fastqPath . "/.$sampleID.passing_samples.$options->{'min_UMIs'}.txt" )
		  or die "I could not create the hidden samples file!\n$!\n";
		foreach ( sort keys %$OK ) {
			print OUT "$_\t$OK->{$_}\n";
		}
		close(OUT);
	}

}
$sum = 0;
map { $sum += $_ } values %$OK;

print "I got "
  . scalar( keys %$OK )
  . " passing samples with a total of $sum reads\n" if ( $options->{'min_UMIs'} > 0 );

$sample_table = data_table->new();
$sample_table->Add_2_Header(
	[ 'sample tag', 'total reads', 'mapped', 'in gene', 'unmapped' ] );
foreach my $cell ( sort keys %$OK ) {
	push( @{ $sample_table->{'data'} }, [ $cell, $OK->{$cell} ] );
}
$sample_table->createIndex('sample tag');

my $runs = 0;
my $adds = 0;
my $UMI;

sub measure_time_and_state {
	my ($msg) = @_;
	$duration = time - $start;
	$duration =
	    ( ( $duration / ( 60 * 60 ) ) % 24 )
	  . " hours, "
	  . ( ( $duration / 60 ) % 60 )
	  . " min and "
	  . ( $duration % 60 )
	  . " seconds";
	print "\n$msg took: $duration s\n";
	print LOG "$msg took: $duration s\n";
	$start = time;
}

sub table_sum {
	my ( $table, @data ) = @_;
	map { $sum += $_ } map {
		unless ($_) { 0 }
		else        { $_ }
	} @data;
}

sub table_min {
	my ( $table, @data ) = @_;
	map { $min = $_ if ( defined $_ and $min > $_ ) } @data;
}

sub table_max {
	my ( $table, @data ) = @_;
	map { $max = $_ if ( $max < $_ ) } map {
		unless ($_) { 0 }
		else        { $_ }
	} @data;
}

sub sample_and_umi_my {
	my (@bam_line) = @_;

#NB501227:57:HGJM5BGX2:1:23206:13379:6742:S_GATCTCAG_C_TATGCCCTCCATTCTA_CAGGTGAAGC
	my ( $sample_name, $UMI );
	my @matching_IDs = split( ":", $bam_line[0] );
	@matching_IDs = split( "_", pop(@matching_IDs) );

	if ( ( !defined $matching_IDs[2] ) or ( !defined $matching_IDs[3] ) ) {
		return sample_and_umi_cellranger(@bam_line);
		Carp::confess(
			"Sorry I lack the UMI information - have you used SplitToCells.pl?"
		);
	}
	$sample_name = $matching_IDs[3];

#$sample_name = $self->{'cell_id_to_name'}->{$sample_name} if ( defined $self->{'cell_id_to_name'}->{$sample_name});
	$UMI = $matching_IDs[4];
	return ( $sample_name, $UMI );
}

sub sample_and_umi_cellranger {
	my (@bam_line) = @_;
	my ( $sample_name, $UMI );
	foreach (@bam_line) {

		if ( $_ =~m/CB:Z:([ACTGN]+)-?\d*$/){
			$sample_name = $1 
		}elsif ( $_ =~ m/CR:Z:([AGCTN]+)$/ ) {
			$sample_name = $1
		}
		if ( $_ =~ m/UB:Z:([ACTGN]+)$/ ){
			$UMI  = $1 ;
		}elsif ( $_ =~ m/UR:Z:([ACTGN]+)$/ ) {
			$UMI  = $1 ;
		}
		
	}
	return ( $sample_name, $UMI );
}

sub filter {
	shift;
	my @bam_line = split( "\t", shift );
	return if ( $bam_line[0] =~ m/^\@/ );
	#die if ( $runs == 100);
	$runs++;

	if ( $runs % 1e4 == 0 ) {
		print ".";
	}
	( $sample_name, $UMI ) = sample_and_umi_my(@bam_line);
	if ( $options->{'min_UMIs'} > 0 ){
		unless ( $OK->{$sample_name} ){
			warn "sample $sample_name should not be analyzed (>$options->{'min_UMIs'} UMIs)\n" if ( $bugfix );
			return;
		}    ## new way to get rid of crap
	}
#	return  if ($GeneModelMatcher->{'match_save'}->{'this_region_start'} > $bam_line[3] );
	Carp::confess(
		"Sorry I lack the UMI information - have you used SplitToCells.pl?\n\t"
		  . join( "\t", @bam_line ) )
	  unless ( defined $UMI and defined $sample_name );
	if ( ! $self->{'UMI'}->{'chr'} eq $bam_line[2] or $self->{'UMI'}->{'start'} + 100 <  $bam_line[3]) {
		delete( $self->{'UMI'} );
		$self->{'UMI'} = { 'chr' =>  $bam_line[2], 'start' => $bam_line[3] };
	}
	$self->{'UMI'}->{$sample_name} ||= { };
	$self->{'UMI'}->{$sample_name}->{$UMI}++;
	$self->{'UMIs'}->{$UMI}++;
	
	## have we seen this UMI for this sample already?
	if ( $self->{'UMI'}->{$sample_name}->{$UMI} > 1 ) {
		$self->{'duplicates'}++;
		warn "\tNo a UMI duplicate ($UMI)\n" if ($bugfix);
		return;
	}

	#next if ( $bugfix and $bam_line[3] == $self->{'last_start'});
	warn
	  "I have got a read for sample $sample_name $bam_line[2], $bam_line[3]\n"
	  if ($bugfix);
	## start with the real matching
	#warn "Not using FLAG==16 here!}\n";
	#TODO: use @bam_line[1] == 16 AND GET THE RIGHT SEQUENCE BACK! giving an edge above cellranger!
	## I get a @bam_line[1] == 16 if the gene is in antisens (orientation -) using 10x data
	&add_to_summary( $sample_name, $GeneModelMatcher->match_cigar( @bam_line[2,3,5,1] ) ) ; ##chr start cigar 

}

#my $pm = Parallel::ForkManager->new($forks);

$runs = 0;

$self->{'end'}        = 0;
$self->{'next_start'} = 0;
$self->{'last_IDS'}   = [];
$self->{'UMI'} = { 'chr' => 'none', 'start' => 0 };

unless ($debug) {    ## debugging the 10x pipeline here
	## turn on autoflush for the process bar:
	if ( -t STDOUT ) {
		$| = 1;
	}
	$bam_file->filter_file( $infile, \&filter );
	$self->{total_reads} = $runs;
	$| = 0;
	&measure_time_and_state("Mapping the UMIs to the transcriptome");

	print
"In total I have processed $self->{'total_reads'} reads and identfied $self->{'duplicates'} UMI duplicates [6bp ("
	  . sprintf( "%.3f",
		( $self->{'duplicates'} / $self->{'total_reads'} ) * 100 )
	  . "%)\n";

##And here is where I stop the whole process this time as all other steps can be done in R.
	if ($asDatabase) {
		$result->print2file( );
#		$spliced->print2file( );
#		$unspliced->print2file( );
	}
	else {
		$result->print2table( );
#		$spliced->print2table( );
#		$unspliced->print2table( );
	}
}
else {
	## as a very short cut I will just create some fake data in the outfiles
	## but also change the outfile (done in the startup)
	my $t = 0;
	for ( my $i = 1 ; $i < 100 ; $i++ ) {
		$result->AddDataset(
			{
				'Gene_ID'   => "Gene000$i",
				'Sample_ID' => "Sample000$i",
				'value'     => $i
			},
			0
		);
		if ($t) {
			$t = 0;
			$result->AddDataset(
				{
					'Gene_ID'   => "Gene000$i",
					'Sample_ID' => "Sample000$i spliced",
					'value'     => $i
				},
				0
			);
		}
		else {
			$t = 1;
		}
	}
	print
"As we are in debug mode data was not quantified, but made up instead ;-)\n";
	$result->print2table();

}
&measure_time_and_state("Saving the outfiles");
print "The file '$outfile' now contains all results from bam file '$infile'\n";

close(LOG);

exit(1);

sub reinit_UMI {
## here I want to not only re-init the UMI has, but I would also like to identify cells that share a UMI in an exon.
## as these different cell_ids are with the highest probability from the same cell, but with a sequencing error.
## And I want to rename these cells for ever.
	$self->{'cell_id_to_name'}     ||= {};
	$self->{'cell_id_to_name_md5'} ||= {};
	$self->{'UMIs'} = {};    #the single UMI storage is not needed here
	my $tmp;
	foreach my $cell_id ( keys %{ $self->{'UMI'} } ) {
		map { $tmp->{$_} ||= {}; $tmp->{$_}->{$cell_id}++ }
		  keys %{ $self->{'UMI'}->{$cell_id} };
	}
	foreach my $UMI ( keys %$tmp ) {
		if ( scalar( keys %{ $tmp->{$UMI} } ) > 1 ) {

#warn "I have found a significant sample overlapp ($UMI) in the samples ". join(", ", keys %{$tmp->{$UMI}} )."\n";
			if ( &check_UMI_overlap( keys %{ $tmp->{$UMI} } ) ) {
				push(
					@{ $self->{'mergable'} },
					[ sort keys %{ $tmp->{$UMI} } ]
				);
			}
		}
	}
	$self->{'next_start'} = 0;
	$self->{'end'}        = 0;
	$self->{'last_IDS'}   = [];
	$self->{'UMI'}        = {};
}



sub add_to_summary {
	my ( $sampleID, $MapResult ) = @_;
	# $MapResult is a hash with keys === 'gene_id' and value IN ('exon', 'spliced' 'primary')
	return if ( scalar( keys %$MapResult) == 0);

	warn "\tI add to sample $sampleID and gene "
	  . join( ",", keys %$MapResult )
	  . "\n"
	  if ($bugfix);
	foreach my $gene_id ( keys %$MapResult) {
		#warn "Quantify - add to sample $sampleID and gene $gene_id one $MapResult->{$gene_id} read\n";
		
		if ( $MapResult->{$gene_id} eq "exon" ) {	
			&add_2( $sampleID, $gene_id, $result);
		}
		elsif ( $MapResult->{$gene_id} eq "spliced" ) {
			&add_2( $sampleID, $gene_id, $result);
			&add_2( $sampleID, $gene_id."_spliced", $result);
		}
		elsif ( $MapResult->{$gene_id} eq "primary" ) {
			&add_2( $sampleID, $gene_id."_primary", $result);
		}
	}
}

sub add_2 {
	my ( $sampleID, $geneID, $where ) = @_;
	
	$where->AddDataset(
			{ 'Gene_ID' => $geneID, 'Sample_ID' => $sampleID, 'value' => 1 },
			1 );    ## add to existing values
	
	if ( $where->Lines() > 1000 ) {
		#print "I am storing 1950 gene results in the tables\n";
		if ($asDatabase) {
			$where->print2file( undef, 700 );
		}
		else {
			$where->print2table( undef, 700 );
		}

	}
}

sub lintersect {
	my $h;

	#warn "lintersect gets ".join(", ",@_)."\n";
	map { $h->{$_}++ } @_;
	my @ret;
	map { push( @ret, $_ ) if ( $h->{$_} > 1 ) } keys %$h;

	#warn "And returns ". join(", ",@ret)."\n";
	return @ret;
}

sub lsum {
	my $sum = 0;
	foreach (@_) {
		$_ ||= 0;
		$sum += $_;
	}
	return $sum;
}

sub lmax {
	return 0 if ( @_ == 0 );
	my $max = shift;
	foreach (@_) {
		$_ ||= 0;
		$max = $_ if ( $_ > $max );
	}
	return $max;
}

sub lmin {
	return 0 if ( @_ == 0 );
	my $min = shift;

	#Carp::confess ( 'min: '.join(", ",@_)."\n");
	foreach (@_) {
		$min = $_ if ( $_ < $min );
	}
	return $min;
}

sub get_chr_end_4_ids {
	my ( $gtf, @lines ) = @_;
	return unless ( scalar(@lines) and defined $lines[0] );
	my @ret = map {
		if ( defined $_ ) {
			@{ @{ $gtf->{'data'} }[$_] }[ $gtf->Header_Position('end') ];
		}
	} @lines;

#	die "the column "
#	 . $gtf->Header_Position('end'). " using ".scalar(@lines)." interesting features "
#	 . " In the gtf file should contain the end position of the feature - right?\nlines: ". join(", ",@lines)."\nto ends:".join(", ",@ret)."\n";
	return @ret;
}

sub get_chr_start_4_ids {
	my ( $gtf, @lines ) = @_;
	return unless ( scalar(@lines) and defined $lines[0] );
	return map {
		if ( defined $_ ) {
			@{ @{ $gtf->{'data'} }[$_] }[ $gtf->Header_Position('start') ];
		}
	} @lines;
}

sub get_reporter_ids {
	my ( $gtf, @lines ) = @_;
	return unless ( scalar(@lines) and defined $lines[0] );
	return &unique(
		map {
			@{ @{ $gtf->{'data'} }[$_] }
			  [ $gtf->Header_Position( $options->{'report_on'} ) ];
		} @lines
	);
}

sub unique {
	my $h;
	map { $h->{$_}++ } @_;
	return sort keys %$h;
}

sub unique_cell_hashes {
	my $h;
	map { $h->{ $_->{'final'} } = $_ } @_;
	return sort { $a->{'final'} cmp $b->{'final'} } values %$h;
}
