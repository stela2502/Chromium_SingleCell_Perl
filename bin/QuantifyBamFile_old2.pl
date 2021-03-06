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

use stefans_libs::GeneModelMatcher;
use Digest::MD5 qw(md5_hex);

use stefans_libs::database::Chromium_SingleCell::datavalues;

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
	$help,    $debug,   $database,  $infile,   $gtf_file, $drop_chr,
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

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}

### initialize default options:

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
mkdir( $fm->{'path'} ) unless ( -d $fm->{'path'} );

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

if ( -f $outfile ) {
	warn
"I am going to add data to the outfile '$outfile'\nremove this file if you want to recreate it.\n";
}

my $result = stefans_libs::result_table->new(
	{
		filename               => $outfile,
		'data_storage_spliced' => 0,

		#	'table_path'           => $drop_chr # not use - might crap things up
	}
);

$fm = root->filemap($outfile);

my $spliced = stefans_libs::result_table->new(
	{
		filename               => $fm->{'path'} ."/", $fm->{'filename_base'}."_spliced.sqlite",
		'data_storage_spliced' => 0,
		#	'table_path'           => $drop_chr # not use - might crap things up
	}
);

my $unspliced = stefans_libs::result_table->new(
	{
		filename               => $fm->{'path'} ."/", $fm->{'filename_base'}."_unspliced.sqlite",
		'data_storage_spliced' => 0,

		#	'table_path'           => $drop_chr # not use - might crap things up
	}
);

my $start       = time;
my $first_start = $start;

my $gtf_obj = stefans_libs::file_readers::gtf_file->new();
if ( defined $drop_chr ) {
	print "I read the chr "
	  . $gtf_obj->_checkChr($drop_chr)
	  . " from the gtf file\n";
	my $filter;
	if ( $drop_chr =~ m/(chr[\w\d\.]*):(\d+)-(\d+)/ ) {
		$gtf_obj->{'tmp_start'} = $2;
		$gtf_obj->{'tmp_stop'}  = $3;
		$gtf_obj->{'tmp_chr'}   = $gtf_obj->_checkChr($1);

		warn
"I am selecting only features in $gtf_obj->{'tmp_chr'}:$gtf_obj->{'tmp_start'}-$gtf_obj->{'tmp_stop'}\n";
		$filter = sub {
			my ( $self, $array, $i ) = @_;
			unless ( defined @$array[0] and @$array[4] ) {
				warn "on line $i something is not defined: '"
				  . join( "', '", @$array[ 1, 3, 4 ] ) . "'\n";
			}
			return (
				( @$array[0] eq $self->{'tmp_chr'} )
				  and ( @$array[3] < $gtf_obj->{'tmp_stop'}
					and @$array[4] > $gtf_obj->{'tmp_start'} )
			);
		};

	}
	else {
		$filter = sub {
			my ( $self, $array, $i ) = @_;
			return @$array[0] eq $self->{'tmp'};
		};
		$gtf_obj->{'tmp'} = $gtf_obj->_checkChr($drop_chr);
	}

	$gtf_obj->read_file( $gtf_file, $filter );
}
else {
	$gtf_obj->read_file($gtf_file);

}

$gtf_obj->{'debug'} = $debug;

$gtf_obj -> GeneFreeSplits(); ## define the gene free splits in order to keep problems with missing genes as low as possible.

$gtf_obj =
  $gtf_obj->select_where( 'feature',
	sub { $_[0] eq $options->{'quantify_on'} } );
	

if ( $gtf_obj->Lines() == 0 ) {
	die "I do not have any information in the gtf file for chr '"
	  . $gtf_obj->_checkChr($drop_chr)
	  . "' - useless aprocess.\n";
}

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
  . " passing samples with a total of $sum reads\n";

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
	print "$msg took: $duration s\n";
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

sub result_areas {
        my ( $start, $cigar ) = @_;
        if ( $cigar =~m/^(\d+)S/ ) {
                $start += $1;
                $cigar =~s/^\d+S//;
        }
        my @r;
        $cigar =~ s/([SNM])(\d)/$1-$2/g;
        foreach ( split("-", $cigar) ){
                if ( $_ =~ m/(\d+)M/ ) {
                        push ( @r, [$start, ($start + $1)]);
                        $start += $1 -1;
                }elsif( $_ =~ m/(\d+)[NS]/ ) {
                        $start += $1 -1;
                }else {
                        die "I do not know what do do with that cigar part: $_\n";
                }
        }
        return @r;
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
	Carp::confess(
		"Sorry I lack the UMI information - have you used SplitToCells.pl?\n\t"
		  . join( "\t", @bam_line ) )
	  unless ( defined $UMI and defined $sample_name );
	$sample_row = undef;
	($sample_row) =
	  $sample_table->get_rowNumbers_4_columnName_and_Entry( 'sample tag',
		$sample_name );
	unless ( defined $sample_row ) {
		$sample_table->AddDataset(
			{
				'sample tag'  => $sample_name,
				'total reads' => "NA",
				'mapped'      => 0,
				'in gene'     => 0,
				'unmapped'    => 0
			}
		);
		($sample_row) =
		  $sample_table->get_rowNumbers_4_columnName_and_Entry( 'sample tag',
			$sample_name );
	}
	unless ( $bam_line[2] =~ m/^chr/ ) {
		@{ @{ $sample_table->{'data'} }[$sample_row] }[4]++;
		return;
	}
	else {
		@{ @{ $sample_table->{'data'} }[$sample_row] }[2]++;
	}
	
	#next if ( $bugfix and $bam_line[3] == $self->{'last_start'});
	warn
	  "I have got a read for sample $sample_name $bam_line[2], $bam_line[3]\n"
	  if ($bugfix);
	## start with the real matching
	
	
	my @read_areas = &result_areas( $bam_line[3], $bam_line[5] );
	my @matching_IDs =  &get_matching_ids( $gtf_obj, $bam_line[2], @{$read_areas[0]}[0], @{$read_areas[0]}[1], 1 ); ## 1 == update
	my @unspliced_IDs;
	@matching_IDs = @{$matching_IDs[0]};
	@unspliced_IDs = @{$matching_IDs[1]} if ( ref($matching_IDs[1]) eq "ARRAY");
	
	if ( @matching_IDs == 0 ) { ## the first part of the match did not match an exon - kick it?
	    if ( @unspliced_IDs > 0 ){
	    	## we have an unspliced read, that does not match another overlapping exon
	    	## or better a read from an unspliced nascent RNA - COOL!
	    	## Store that somewhere
	    	my @matching_Names = &get_reporter_ids( $gtf_obj, @unspliced_IDs );
	    	map{ &add_2( $sample_name, $_, $unspliced ) } @matching_Names;
	    }
		warn "No matching features!\n"if ( $bugfix );
		return; # yes kick it!
#		$Seq = 0;
#		map { $Seq += $_ }
#		  $bam_line[5] =~ m/(\d+)[NMS]/g
#		  ;  ## the fist match should be in the exon and overlap to at least 50%
#		     #$Seq = ceil($Seq/2);
#		@matching_IDs =
#		  &get_matching_ids( $gtf_obj, $bam_line[2], $bam_line[3] + $Seq ); ## do not update!
	#	warn
	#  "I have got a read for sample $sample_name $bam_line[2], $bam_line[3] and ".($bam_line[3] + $Seq)."\n"
	}

	warn "\tGot a matching ID? " . join( ", ", @matching_IDs ) . "\n"
	  if ($bugfix);

	## I got a match and therfore need to store the sample info and UMI info
	@{ @{ $sample_table->{'data'} }[$sample_row] }[2]++;

	$self->{'UMI'}->{$sample_name} ||= {};
	$self->{'UMI'}->{$sample_name}->{$UMI}++;
	$self->{'UMIs'}->{$UMI}++;
	
	## have we seen this UMI for this sample already?
	if ( $self->{'UMIs'}->{$UMI} > 1 ) {
		$self->{'duplicates'}++;
		warn "\tNo a UMI duplicate ($UMI)\n" if ($bugfix);
		return if $self->{'UMI'}->{$sample_name}->{$UMI};
	}

	my @matching_Names = &get_reporter_ids( $gtf_obj, @matching_IDs );

	if ( @matching_Names > 0 ){
		## do we have a spliced read?
		if ( @read_areas > 1 ) {
			## So the real match can only be on an element on both ends of the splice
			foreach my $region ( @read_areas[1..$#read_areas] ) {
			@matching_Names = &lintersect(
			@matching_Names,
			&get_reporter_ids(
				$gtf_obj,
				$gtf_obj->efficient_match_chr_position(
					$bam_line[2],
					@$region
				)
			)
			);
			}
			if ( @matching_Names == 0 ) {
				##This does mean I have a read starting in a exon, but ending in an intron
				## https://www.biorxiv.org/content/early/2017/10/19/206052
				## I probably should collect this info, too
				
				return;
			}
			else{
				&add_to_summary( $sample_name, \@matching_Names, 1, "not used anyhow" );
			}
			
		}else { ## no spiced read
			&add_to_summary( $sample_name, \@matching_Names, 0, "$bam_line[3] and $bam_line[5]" );
		}
	}
}



#my $pm = Parallel::ForkManager->new($forks);

$runs = 0;

$self->{'end'}        = 0;
$self->{'next_start'} = 0;
$self->{'last_IDS'}   = [];

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
"In total I have processed $self->{'total_reads'} and identfied $self->{'duplicates'} UMI duplicates [6bp ("
	  . sprintf( "%.3f",
		( $self->{'duplicates'} / $self->{'total_reads'} ) * 100 )
	  . "%)\n";

##And here is where I stop the whole process this time as all other steps can be done in R.
	if ($asDatabase) {
		$result->print2file( );
		$spliced->print2file( );
		$unspliced->print2file( );
	}
	else {
		$result->print2table( );
		$spliced->print2table( );
		$unspliced->print2table( );
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
&measure_time_and_state("Saving th outfiles");
print "The file '$outfile' now contains all results from bam file '$infile'\n";

close(LOG);

exit(1);

sub get_matching_ids {
	my ( $gtf, $chr, $start, $end, $update ) = @_;
	my ( @IDS, $nextID, @unspliced );
	unless ( $self->{'chr'} eq $chr ) {
		$self->{'chr'}        = $chr;
		$self->{'end'}        = 0;
		$self->{'next_start'} = 0;
		$self->{'skip'}       = '-not-a-chr-';
		$self->{'last_start'} = 0;
		&reinit_UMI();
		warn "Quant init chr $chr\n" if ($bugfix);
	}
	if ( $self->{'last_start'} > $end ) {
		## shit - that should not be possible!
		## better save than sorry:
		#warn "reinit end and next_start (".($self->{'last_start'} - $end)."bp)";
		$self->{'end'} = $self->{'next_start'} = 0;
	}
	$self->{'last_start'} = $start if ( $update );
	return ([],[]) if ( $self->{'skip'} eq $chr );
	if ( $self->{'end'} < $start and $start < $self->{'next_start'} ) {
		warn "get_matching_ids intergenic - return ($self->{'end'} < $start) & ($start < $self->{'next_start'})\n" if ($bugfix);
		return ([],[]);
	}
	#if ( $self->{'end'} <= $start ) {    # match to a new feature
	if ( $start < $self->{'end'} ) {
		warn
"\tNo match needed ( $start < $self->{'end'} ) I can use the old match\n"
		  . join( "\n", @{ $self->{'last_IDS'} } ) . "\n"
		  if ($bugfix);
		@IDS = @{ $self->{'last_IDS'} };
	}
	elsif ( $start >= $self->{'next_start'}) {
		warn
"Do a new match : get_matching_ids end < start ($self->{'end'} < $start)\n"
		  if ($bugfix);
		&reinit_UMI();
		$self->{'last_IDS'} = [];
		warn "lib change using $start -50!\n";
		@IDS = $gtf->efficient_match_chr_position_plus_one( $chr, $start-50, $end );
		if ( scalar(@IDS) == 1 and !defined $IDS[0] ) {
			## useless chr!
			$self->{'skip'} = $chr;
			return ();
		}
		my ( $match, $next, @OK, $lastOK, @starts, @ends );
		@starts = &get_chr_start_4_ids( $gtf, @IDS );
		@ends = &get_chr_end_4_ids( $gtf, @IDS );
		$lastOK = 0;
		if ($update) {
			$self->{'end'} =
			  $start;    ## we work on ordered reads hence this is OK
			$self->{'next_start'} = undef;
			for ( my $i = 0 ; $i < scalar(@IDS) ; $i++ )
			{    # @IDs also conatin the next entry after all matching ones!
				if ( $starts[$i] <= $start and $ends[$i] >= $start ) {
					## matching
					#warn "\tAnd I get a match!\n";
					push( @OK, $IDS[$i] );
					$self->{'end'} = $ends[$i];
				}
				if ( $starts[$i] > $start ) {
					$self->{'next_start'} = $starts[$i];
				}
			}
			## so if I am not totally missguided we should have some values now:
			$self->{'next_start'} ||= $starts[$#starts];
		}
		else {       ## no update
			for ( my $i = 0 ; $i < scalar(@IDS) ; $i++ )
			{    # @IDs also conatin the next entry after all matching ones!
				if ( $starts[$i] <= $start and $ends[$i] >= $start ) {
					## matching
					push( @OK, $IDS[$i] );
				}elsif ( $start < $ends[$i] and $end > $starts[$i] ) {
					push ( @unspliced, $IDS[$i] );
				}
			}
		}
		if ( !defined $self->{'next_start'} ){
			## OK that can only mean, that with this mapper there are no new entries untill the next chr mapper would be used.
			my $chr_id = $gtf->get_chr_subID_4_start($chr, $start);
			$self->{'next_start'} = @{ @{ $gtf->{'__GeneFreeSplits__'}->{'data'} }[$chr_id] }[ 2 ]; # the end of the area covered by the mapper.
		}

		@starts = @ends = undef;
		if ( $bugfix and scalar( @IDS == 0 ) ) {

			#			## seams to be fixed!
			warn
"If all is OK I now should have a $self->{'end'} that is at least the $start ("
			  . ( $self->{'end'} - $start )
			  . ") and a $self->{'next_start'} for the next possible exon.\n"
			  . "I got these matches for $start OK("
			  . join( ", ", @OK ) . "):\n"
			  . join( "\n",
				map { join( "\t", @{ @{ $gtf->{'data'} }[$_] } ) } @IDS )
			  . "\n";
		}

		@IDS = @OK;
	}
	elsif ( $start < $self->{'next_start'} ) {
		warn "\tNo match needed?!( $start < $self->{'next_start'} )\n"
		  if ($bugfix);
		@IDS = ();    ## no match needed...
	}
	else {
		die "I assume I should never get here!?\n";
	}
	$self->{'last_IDS'} = [@IDS];
	return \@IDS, \@unspliced ;
}

sub reinit_UMI {
## here I want to not only re-init the UMI has, but I would also like to identify cells that share a UMI in an exon.
## as thesedifferent cell_ids are with the highest probability from the same cell, but with a sequencing error.
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

sub check_UMI_overlap {
	my (@cells) = @_;
	my ( $m, @ids );
	$m = 1;
	map { $m = 0 unless defined $_ } @cells;
	die "are there not defined cells?: '" . join( "', '", @cells ) . "'\n"
	  unless ($m);
	for ( my $i = 0 ; $i < @cells ; $i++ ) {
		for ( my $a = 1 ; $a < @cells ; $a++ ) {
			next if ( $cells[$i] eq $cells[$a] );
			$m = join( " ", sort @cells[ $i, $a ] );
			$ids[@ids] = $self->{'cell_id_to_name_md5'}->{$m} += 1;
		}
	}
	foreach (@ids) {
		return 1 if ( $_ == 1 );
	}
	return 0;
}

sub stringdist {
	my ( $a, $b ) = @_;
	my @A = split( "", $a );
	my @B = split( "", $b );
	my $ret = 0;
	for ( my $i = 0 ; $i < @A ; $i++ ) {
		if ( !defined $A[$i] or !defined $B[$i] ) {
			Carp::confess("String problems at position '$i' a:'$a' b:'$b'\n");
		}
		$ret++ unless ( $A[$i] eq $B[$i] );
	}
	return $ret;
}

sub identify_merge_target {
	my $cname = shift;
	my $return;

	#warn "\nYou try to identify the merge target for '$cname'\n";
	eval {
		($return) =
		  $sample_table->get_value_for( 'sample tag', $cname, 'merged to' );
	};

#warn   "                 I got the return value  '$return'  stringdist = ".&stringdist($cname,$return)."\n";
	unless ( defined $return ) {

		#$sample_table->write_table("problematic_samples_missing_$cname");
		# this function is used as test too - hence I must not die here!
		#Carp::confess("Could not identify sample $cname!\n");
		return undef;
	}
	if ( $return =~ m/C_[AGTCN]*/
		and defined $result->Header_Position($return) )
	{
		#warn "identify_merge_target identified the target $cname -> $return\n";
		return $return;
	}
	elsif ( $return =~ m/C_[AGTCN]*/ )
	{    ## most likely a additional target for this sample
		 #warn "The secondary target $return has also been merged! - re-init\n";
		return &identify_merge_target($return);
	}
	Carp::confess(
"identify_merge_target($cname) -> You should not have reached this point! '$return'\n"
	);
}

sub merge_cells {
	return if ( @_ < 2 );
	if ( $runs % 1e+5 == 0 ) {
		print ".";
	}

	#	if ( $runs % 6e+6 == 0 ) {
	#		print "\n";
	#	}
	my @cells = @_;

	#warn "I am going to merge the cells:\n\t".join("\n\t", @cells)."\n";
	## first check that these columns still exist!
	my ( @tmp, @already_merged, $mtcol, $lineHash );

	foreach my $cell_name (@cells) {    #the sample has not been merged before
		if (    defined $result->Header_Position($_)
			and defined $sample_table->createIndex('sample tag')->{$cell_name} )
		{
			$tmp[@tmp] = { 'orig' => $cell_name, 'final' => $cell_name };
		}
		elsif ( defined $sample_table->createIndex('sample tag')->{$cell_name} )
		{                               ## one partner has been merged before
			## now we need to try to identify the respective target and check that it is not the other sample
			$mtcol = &identify_merge_target($cell_name);
			unless ( defined $mtcol ) {

				# the cell was simply not in the dataset - strange!
				# no action needed?!
				next;
			}
			if ( $mtcol =~ m/C_[AGTCN]*/ ) {
				$tmp[@tmp] = { 'orig' => $cell_name, 'final' => $mtcol };
			}
		}
		else {
#warn "I assume one of the partners had not a single unique mapped read and was therefore not kept:".join(", ", @cells)."\n";
			return;
		}
	}
	@tmp = &unique_cell_hashes(@tmp);
	if ( @tmp < 2 ) {

#warn "I could not merge the cells ". join( ", ", @cells ). " as I only identified ". scalar(@tmp). " rows I could merge\n";
		return;
	}
	@cells = @tmp;

#Carp::confess ( "is this a real sable cell has array?". root->print_perl_var_def( @cells)."\n");

	my @sums = map {
		$sample_table->get_value_for(
			'sample tag',
			$_->{'final'},
			'total UMI count'
		);
	} @cells;

	#warn "the sums: '".join("', '", @sums ) ."'\n";
	my $max = &lmax(@sums);
	my ( $main_cell, $position, @slave_cols, $line_id );
	for ( my $i = 0 ; $i < @sums ; $i++ ) {
		if ( $sums[$i] == $max ) {
			$position  = $i;
			$main_cell = $cells[$i];
			splice( @cells, $i, 1 );    ## kick the merge target
			last;
		}
	}
	$sample_table->set_value_for( 'sample tag', $main_cell, 'total UMI count',
		&lsum(@sums) );
	foreach my $cell_hash (@cells) {
		$lineHash = $sample_table->get_line_asHash(
			$line_id = $sample_table->get_rowNumbers_4_columnName_and_Entry(
				'sample tag', $cell_hash->{'orig'}
			)
		);
		$lineHash->{'merged to'}      = $main_cell->{'final'};
		$lineHash->{'merge evidence'} = $self->{'cell_id_to_name_md5'}
		  ->{ join( " ", sort $main_cell->{'orig'}, $cell_hash->{'orig'} ) };
		$lineHash->{'stringdist'} =
		  &stringdist( $main_cell->{'final'}, $cell_hash->{'final'} );
		$lineHash->{'mapped	in gene'} ||= 1;
		if ( $lineHash->{'merge evidence'} / $lineHash->{'mapped	in gene'} >
			0.05 )
		{
			#merge this cell
		}
		else {
			root->print_perl_var_def(
				{ 'target' => $main_cell, 'source' => $cell_hash } )
			  . "\nI assume I should not merge this cell:\n"
			  . join( "\t", @{ $sample_table->{'header'} } ) . "\n"
			  . join( "\t", @{ @{ $sample_table->{'data'} }[$line_id] } )
			  . "\nto the master cell (merge evidence = $lineHash->{'merge evidence'}; stringdist=$lineHash->{'stringdist'}):\n"
			  . join(
				"\t",
				@{
					@{ $sample_table->{'data'} }[
					  $sample_table->get_rowNumbers_4_columnName_and_Entry(
						  'sample tag', $main_cell->{'final'} )
					]
				}
			  ) . "\n";
		}
		$sample_table->set_value_for( 'sample tag', $cell_hash, 'merged to',
			$main_cell->{'final'} );
		$sample_table->set_value_for(
			'sample tag', $cell_hash,
			'merge evidence',
			$lineHash->{'merge evidence'}
		);
		$sample_table->set_value_for( 'sample tag', $cell_hash, 'stringdist',
			&stringdist( $main_cell->{'final'}, $cell_hash->{'final'} ) );
		$lineHash = $sample_table->get_line_asHash(
			$sample_table->get_rowNumbers_4_columnName_and_Entry(
				'sample tag', $cell_hash->{'final'}
			)
		);

	}

	$self->{'PDL'}
	  ->slice( ":," . $self->{'samples_in_result'}->{ $main_cell->{'final'} } )
	  += $self->{'PDL'}->slice(
		":,"
		  . join( ":",
			map { $self->{'samples_in_result'}->{ $_->{'final'} } } @cells )
	  );

	$self->{'PDL'}->slice( ":,"
		  . $self->{'samples_in_result'}->{ $main_cell->{'final'} . " spliced" }
	  ) += $self->{'PDL'}->slice(
		":,"
		  . join(
			":",
			map { $self->{'samples_in_result'}->{ $_->{'final'} . " spliced" } }
			  @cells
		  )
	  );
	foreach my $cell_hash (@cells) {

		$result->rename_column( $cell_hash->{'final'},
			"merged $cell_hash->{'final'}" );
		$result->rename_column(
			$cell_hash->{'final'} . " spliced",
			"merged $cell_hash->{'final'} spliced"
		);
	}

}

sub add_2 {
	my ( $sampleID, $geneIDs, $where ) = @_;
	foreach my $gene_id (@$geneIDs) {
		$where->AddDataset(
			{ 'Gene_ID' => $gene_id, 'Sample_ID' => $sampleID, 'value' => 1 },
			1 );    ## add to existing values
	}
	if ( $where->Lines() > 1000 ) {
		#print "I am storing 1950 gene results in the tables\n";
		if ($asDatabase) {
			$where->print2file( undef, 950 );
		}
		else {
			$where->print2table( undef, 950 );
		}

	}
}

sub add_to_summary {
	my ( $sampleID, $geneIDs, $is_spliced, $start ) = @_;
	$is_spliced = 1 if ($is_spliced);

	warn "\tI add to sample $sampleID and gene "
	  . join( ",", @$geneIDs )
	  . " is_splice==$is_spliced\n"
	  if ($bugfix);
	foreach my $gene_id (@$geneIDs) {
		&add_2( $sampleID, $geneIDs, $result);

		if ($is_spliced) {
			&add_2( $sampleID, $geneIDs, $spliced);
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
