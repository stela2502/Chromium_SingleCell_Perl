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
       -options  : format: key_1 value_1 key_2 value_2 ... key_n value_n

           quantify_on: which feature to select for quantification from the gtf file (default 'exon')
             report_on: which column in the gtf file object to report on (default 'gene_id')
              min_UMIs: how many UMIs have to mapp to any gene to report the sample (default 100)

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Take a Chromium 10x bam file and quantify it agains a gtf file.

  To get further help use 'QuantifyBamFile.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::BAMfile;
use stefans_libs::flexible_data_structures::data_table;
use stefans_libs::file_readers::gtf_file;
use Digest::MD5 qw(md5_hex);

use PDL;
use PDL::NiceSlice;
$PDL::BIGPDL = 1;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, $infile, $gtf_file, $outfile, $options,
	@options );

Getopt::Long::GetOptions(
	"-infile=s"     => \$infile,
	"-gtf_file=s"   => \$gtf_file,
	"-outfile=s"    => \$outfile,
	"-options=s{,}" => \@options,

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

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}

### initialize default options:

$options->{'quantify_on'} ||= 'exon';
$options->{'report_on'}   ||= 'gene_id';
$options->{'min_UMIs'}    ||=  100;

###

my $self = {
	'chr'      => '',
	'end'      => 0,
	'UMI'      => {},
	'mergable' => [],
	'UMIs'     => {},
};    # the script should store global variables here

## turn on autoflush for the process bar:
$| = 1;

my $fm = root->filemap($outfile);
mkdir( $fm->{'path'} ) unless ( -d $fm->{'path'} );

open( LOG, ">$outfile.log" ) or die $!;
print LOG $task_description . "\n";


## Do whatever you want!

my $bam_file = stefans_libs::BAMfile->new();
my $result   = data_table->new();
$result->Add_2_Header('Gene_ID');
$result->createIndex('Gene_ID');

my $start = time;
my $first_start = $start;

my $gtf_obj = stefans_libs::file_readers::gtf_file->new();

$gtf_obj->read_file($gtf_file);

my $quantifer =
  $gtf_obj->select_where( 'feature',
	sub { $_[0] eq $options->{'quantify_on'} } );

my (
	@matching_IDs, $N,   $Seq, $sample_name, $sample_table,
	$sample_row,   $sum, $max, $min
);

my $duration = time - $start;
print "Preparation of the gtf file took: $duration s\n";
print LOG "Preparation of the gtf file took: $duration s\n";
$start = time;

$sample_table = data_table->new();
$sample_table->Add_2_Header(
	[ 'sample tag', 'mapped', 'in gene', 'unmapped' ] );
$sample_table->createIndex('sample tag');

my $runs = 0;
my $adds = 0;
my $UMI;

sub measure_time_and_state {
	my ( $msg ) = @_;
	$duration = time - $start;
	$duration = (($duration/(60*60))%24) . " hours, ".(($duration/60)%60)." min and ".($duration%60)." seconds";
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

sub filter {
	shift;
	my @bam_line = split( "\t", shift );
	return if ( $bam_line[0] =~ m/^\@/ );
	$runs++;

	if ( $runs % 1e4 == 0 ) {
		print ".";
	}

#NB501227:57:HGJM5BGX2:1:23206:13379:6742:S_GATCTCAG_C_TATGCCCTCCATTCTA_CAGGTGAAGC
	@matching_IDs = split( ":", $bam_line[0] );
	@matching_IDs = split( "_", pop(@matching_IDs) );

	$sample_name = join( "_", @matching_IDs[ 2, 3 ] );

#$sample_name = $self->{'cell_id_to_name'}->{$sample_name} if ( defined $self->{'cell_id_to_name'}->{$sample_name});
	$UMI = $matching_IDs[4];

	$sample_row = undef;
	($sample_row) =
	  $sample_table->get_rowNumbers_4_columnName_and_Entry( 'sample tag',
		$sample_name );
	unless ( defined $sample_row ) {
		$sample_table->AddDataset(
			{
				'sample tag' => $sample_name,
				'mapped'     => 0,
				'in gene'    => 0,
				'unmapped'   => 0
			}
		);
		($sample_row) =
		  $sample_table->get_rowNumbers_4_columnName_and_Entry( 'sample tag',
			$sample_name );
	}
	unless ( $bam_line[2] =~ m/^chr/ ) {
		@{ @{ $sample_table->{'data'} }[$sample_row] }[3]++;
		return;
	}
	else {
		@{ @{ $sample_table->{'data'} }[$sample_row] }[1]++;
	}

	## start with the real matching
	@matching_IDs = &get_matching_ids( $quantifer, $bam_line[2], $bam_line[3] );

	return if ( @matching_IDs == 0 );    ## no match to any gene / exon

	@{ @{ $sample_table->{'data'} }[$sample_row] }[2]++;

	$self->{'UMI'}->{$sample_name} ||= {};
	$self->{'UMI'}->{$sample_name}->{$UMI}++;
	$self->{'UMIs'}->{$UMI}++;
	return
	  if ( $self->{'UMIs'}->{$UMI} > 1 )
	  ;    ## the sample specific UMIs will lead to a merge of the samples
	##and therefore I suold not add more than one UMI per exon

	@matching_IDs = &get_reporter_ids( $quantifer, @matching_IDs );
	$N            = 0;
	$Seq          = 0;
	map { $N   += $_ } $bam_line[4] =~ m/(\d+)N/g;
	map { $Seq += $_ } $bam_line[4] =~ m/(\d+)M/g;

	if ( $N > $Seq ) {
		## most probably spliced! (?)
		## So the real match can only be on an element on both ends of the splice
		$Seq = 0;
		map { $Seq += $_ } $bam_line[4] =~ m/(\d+)/g;
		@matching_IDs = &lintersect(
			&get_reporter_ids(
				$quantifer,
				$quantifer->efficient_match_chr_position(
					$bam_line[2], $bam_line[3] + $Seq
				)
			)
		);
		&add_to_summary( $sample_name, @matching_IDs, 1 );
	}
	else {
		&add_to_summary( $sample_name, @matching_IDs, 0 );
	}

}

$bam_file->filter_file( $infile, \&filter );

&measure_time_and_state("Mapping the UMIs to the transcriptome");

## now I can drop the quantifier and the gff

$quantifer = undef;
$gtf_obj   = undef;

### done with the quantification - now the error check starts

my ( @tmp, $sample_row_id, @export, @sampleNames, $tmp );

$| = 0;    ## buffer file output again (default)

if ( $outfile =~ m/txt$/ ) {
	$outfile =~ s/txt$/xls/;
}
unless ( $outfile =~ m/xls$/ ) {
	$outfile .= ".xls";
}
$tmp = $outfile;
$tmp =~ s/.xls$//;

if ($debug) {
	$result->write_table( $tmp . ".original.xls" );
}

$runs = 0;
print
"\nFinished with quantification (".$sample_table->Rows()." samples and ".$result->Rows()." genes)\n";

## calculate the sample UMI count and add it to the summary table


if ($debug) {
	$sample_table->write_table( $tmp . "sample_table_original" );
	&measure_time_and_state("DEBUG: Writing of the raw data files");
}


print "Syncing sample table and data table\n";
my $i = 0;
$self->{'samples_in_result'} = { map {$_ => $i ++} @{ $result->{'header'} } };

$sample_table->drop_rows( 'sample tag', sub{ ! defined $self->{'samples_in_result'}->{$_[0]} } );

if ($debug) {
	$sample_table->write_table( $tmp . "sample_table_only_used_samples" );
}

$sample_table = $sample_table->Sort_by([ [  'sample tag', $self->{'samples_in_result'} ] ] );


if ( $debug ){
	my $tmp_table = data_table->new();
	$tmp_table -> Add_2_Header( ['samples', 'results'] );
	my $OK =1;
	for ( my $i = 0; $i < $sample_table->Columns();$i ++ ){
		#print "new entry in $tmp_table->{'data'} at position $i = \n";
		#print "\t [ \n\n\t\t".@{@{$sample_table->{'data'}}[$i]}[0].", \n";
		#print "\t\t".@{$result->{'header'}}[($i*2) +1 ]." ];\n";
		@{$tmp_table->{'data'}}[$i] = [ @{@{$sample_table->{'data'}}[$i]}[0], @{$result->{'header'}}[($i*2) +1 ] ];
		$OK = 0 unless ( @{@{$sample_table->{'data'}}[$i]}[0] eq @{$result->{'header'}}[ ($i*2) +1 ] );
	}
	die "The order has not been chenged in the right way!\n" unless ( $OK );
}

$sample_table->write_table( $tmp . "sample_table_original" );

&measure_time_and_state("Syncing the sample and data tables");

delete( $result->{'index_length'}->{'Gene_ID'});
delete( $result->{'index'}->{'Gene_ID'});

my $Gene_IDs = $result->GetAsArray('Gene_ID');
$result->drop_column('Gene_ID');
foreach ( keys %{$self->{'samples_in_result'}} ){
	$self->{'samples_in_result'}->{$_} --;
}
$| = 1; 
print "Calculating colum sums for ". ($sample_table->Rows()/2)." samples:\n";
my $data_col_name;


$self->{PDL} = $result->GetAsPDL();

&measure_time_and_state("Creating the PDL");

my $sums = $self->{PDL}->sumover();
$sums = unpdl($sums);

$sample_table->add_column( 'total UMI count', @$sums[( map{$_ * 2 } 0..($sample_table->Rows()-1) )] );
$sample_table->add_column( 'total UMI count no merge', @$sums[( map{$_ * 2 } 0..($sample_table->Rows()-1) )] );
$sample_table->add_column( 'spliced UMI count', @$sums[ map{$_ * 2 +1} 0..($sample_table->Rows()-1) ] );
$sample_table->add_column( 'spliced UMI count no merge', @$sums[ map{$_ * 2 +1} 0..($sample_table->Rows()-1) ] );

&measure_time_and_state("Calculate the sum over ".$sample_table->Rows()." samples" );

## now I need to merge the samples!

print "\nStarting to merge cells that share a UMI in an exon\n";

$sample_table->Add_2_Header(['merged to','merge evidence'] );

foreach ( @{ $self->{'mergable'} } ) {
	&merge_cells(@$_);
}
delete($self->{'mergable'});

&measure_time_and_state("Mering cells where the same UMI tagged the same transcript");

## as the new merge function works on the PDL data I now need to copy the pdl data back to the result table
$result->{'data'} = unpdl($self->{'PDL'}->xchg(0,1));
$result->add_column( 'Gene_ID', @$Gene_IDs);

&measure_time_and_state("Restoring data table from merged PDL");

$sample_table->write_table( $tmp . ".samples_after_sum_add.xls" );

$sample_table->drop_rows('merged to', sub{ defined $_[0] } );


$| = 0; 

if ($debug) {
	open( my $MERGE, ">$tmp.merge.log" ) or die $!;
	print $MERGE "to\t" . join( "\t", map { "from $_" } 1 .. 40 ) . "\n";
	foreach my $target ( sort keys %{ $self->{'cell_id_to_name'} } ) {
		print $MERGE $target . "\t"
		  . join( "\t", sort keys %{ $self->{'cell_id_to_name'}->{$target} } )
		  . "\n";
	}
	close($MERGE);
}

if ($debug) {
	$result->write_table( $tmp . ".original_merged.xls" );
}

$result->drop_column( [ grep (/merged/, @{$result->{'header'}})]);

print "\ncollapse the data to the unique cells only\n";

## write the data

$result->define_subset( 'keep',
	[ grep ( !/merged/, @{ $result->{'header'} } ) ] );
@export      = $result->Header_Position('keep');
@sampleNames = @{ $result->{'header'} }[@export];

$sample_table = $sample_table->select_where('total UMI count',	sub { ( defined $_[0] and $_[0] > 0 ) } ) ;

my $columns_to_drop = $sample_table->GetAsHash('sample tag', 'total UMI count');
$sample_table = $sample_table->select_where('total UMI count',	sub { ( defined $_[0] and $_[0] > $options->{'min_UMIs'}-1 ) });


my @columns_to_drop;
foreach ( keys %$columns_to_drop) {
	if ( $columns_to_drop->{$_} <  $options->{'min_UMIs'} ) {
		push ( @columns_to_drop , $_,  "$_ spliced" ) if ( defined $result->Header_Position($_) );
	}
}
## in the test data I have identified 199748 samples with a total of 0 UMIS in a gene - totally useless to keep them!
## the data table has to be modified in place: drop_columns?
$result->drop_column( \@columns_to_drop ) if ( @columns_to_drop > 0);

&measure_time_and_state("dropping merged samples from the result and samples tables");


eval {
	open ( my $Mreport , ">$tmp.merge_report.txt" );
	foreach ( keys %{$self->{'cell_id_to_name_md5'}} ){
		$self->{'cell_id_to_name_md5'}->{$_} ||= 0;
		print $Mreport "$self->{'cell_id_to_name_md5'}->{$_} $_\n";
	}
	close ( $Mreport );
};

&measure_time_and_state("Writing the merge report to '$tmp.merge_report.txt'");


my $msg= "\nFinished with mapping: "
  . $sample_table->Rows()
  . " samples and "
  . $result->Rows()
  . " $options->{'report_on'}'s detected\n";

print $msg;

$sample_table->write_table( $tmp . ".samples.xls" );

print LOG $msg;

#previousely run:
#$result->define_subset( 'all reads', [ grep ( !/spliced/, @sampleNames )]);
#$result->define_subset( 'spliced reads', [ $sampleNames[0], grep ( /spliced/, @sampleNames )]);

@sampleNames = @{ $result->{'header'} };
$result->define_subset( 'all reads', [ grep ( !/spliced/, @sampleNames ) ] );
$result->define_subset( 'spliced reads',
	[ $sampleNames[0], grep ( /spliced/, @sampleNames ) ] );

my ( @all_reads, @spliced_reads );
@all_reads     = $result->Header_Position('all reads');
@spliced_reads = $result->Header_Position('spliced reads');

open( my $OUT, ">$outfile" )
  or die "I could not create the outfile '$outfile'\n$!\n";

open( my $SPLICED, ">$tmp.spliced.xls" )
  or die "I could not create the spliced reads data '$tmp.spliced.xls'\n$!\n";

print $OUT join( "\t", @{ $result->{'header'} }[@all_reads] ) . "\n";
print $SPLICED join( "\t", @{ $result->{'header'} }[@spliced_reads] ) . "\n";

for ( my $i = 0 ; $i < $result->Rows() ; $i++ ) {
	print $OUT join(
		"\t",
		map {
			unless ($_) { '' }
			else        { $_ }
		} @{ @{ $result->{'data'} }[$i] }[@all_reads]
	) . "\n";
	print $SPLICED join(
		"\t",
		map {
			unless ($_) { '' }
			else        { $_ }
		} @{ @{ $result->{'data'} }[$i] }[@spliced_reads]
	) . "\n";

}

close($OUT);
close($SPLICED);

&measure_time_and_state("Saving data to '$outfile' and '$tmp.spliced.xls'");

close(LOG);

print "Done\n";
$start = $first_start;
&measure_time_and_state("Total run");


sub get_matching_ids {
	my ( $gtf, $chr, $start ) = @_;
	my ( @IDS, $nextID );
	unless ( $self->{'chr'} eq $chr ) {
		$self->{'chr'} = $chr;
		$self->{'end'} = 0;

		#$self->{'skip_to_next_chr'} = 0;
		&reinit_UMI();
	}

#return if ( $self->{'skip_to_next_chr'} ); # we would not have annotations anyhow...
	if ( $self->{'end'} <= $start ) {
		&reinit_UMI();
		@IDS = $gtf->efficient_match_chr_position_plus_one( $chr, $start );

		#warn(
		#	"Some fuckup happened - I did not get any result ($chr, $start)\n'"
		#	  . join( "','", @IDS )
		#	  . "\n" )
		return
		  if ( @IDS == 0 or !defined( $IDS[0] ) )
		  ;    ## pdl has no more matching entries
		$nextID = pop(@IDS);
		unless ( defined $nextID )
		{      ## the end of the chromosome annotation data has been reached
			    #$self->{'skip_to_next_chr'} = 1;
			$self->{'next_start'} =
			  ( $gtf->get_chr_subID_4_start($start) + 1 ) *
			  $gtf->{'slice_length'};
			return;
		}
		$self->{'next_start'} =
		  @{ @{ $gtf->{'data'} }[$nextID] }[ $gtf->Header_Position('start') ];
		$self->{'end'} = &lmin( &get_chr_end_4_ids( $gtf, @IDS ) );
		$self->{'last_IDS'} = \@IDS;
	}
	elsif ( $start < $self->{'next_start'} ) {
		@IDS = ();    ## no match needed...
	}
	else {
		## the read is still in the old matching area
		@IDS = @{ $self->{'last_IDS'} };
	}
	return @IDS;
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
			my $m =  join( " ",sort keys %{ $tmp->{$UMI} } ) ;
			$self->{'cell_id_to_name_md5'}->{$m} ++;
			next if ( $self->{'cell_id_to_name_md5'}->{$m} > 1 );
			push( @{ $self->{'mergable'} }, [sort keys %{ $tmp->{$UMI} } ] );
		}
	}
	$self->{'UMI'} = {};
}

sub merge_cells {
	return if ( @_ < 2 );
	if ( $runs % 1e5 == 0 ) {
		print ".";
	}
	my @cells = @_;
	my $evidence = $self->{'cell_id_to_name_md5'}->{ join(" ", @_) } || 0;

	#warn "I am going to merge the cells:\n\t".join("\n\t", @cells)."\n";
	## first check that these columns still exist!
	my @tmp;
	my @already_merged;
	foreach (@cells) {
		if ( defined $result->Header_Position($_) and defined  $sample_table->createIndex( 'sample tag' )->{$_} ){
			$tmp[@tmp]= $_ ;
		}else {
			$already_merged[@already_merged] = $_;
		}
	}
	if ( @tmp < 2 ) {
		return;
	}
	@cells = @tmp;
	#Carp::confess ( "is this a real index?". root->print_perl_var_def( 
	#{index => $sample_table->createIndex( 'sample tag' ) } )."\n");
	
	my @sums = map{$sample_table->get_value_for( 'sample tag', $_, 'total UMI count' ) } @cells;
	#warn "the sums: '".join("', '", @sums ) ."'\n";
	my $max = &lmax(@sums);
	my ( $main_cell, $position, @slave_cols);
	for ( my $i = 0 ; $i < @sums ; $i++ ) {
		if ( $sums[$i] == $max ) {
			$position  = $i;
			$main_cell = $cells[$i];
			splice( @cells, $i ,1); ## kick the merge target
			last;
		}
	}
	$sample_table->set_value_for( 'sample tag', $main_cell, 'total UMI count',
		&lsum(@sums) );
	map { 
		$sample_table->set_value_for( 'sample tag', $_ , 'merged to', $main_cell );
		$sample_table->set_value_for( 'sample tag', $_ , 'merge evidence', $evidence );
	} @cells;
	
#	Carp::confess "\$exp = ".root->print_perl_var_def( {
#		'PDL shape' => unpdl ( $self->{'PDL'}->shape()),
#		'samples_table shape' => [ $sample_table->Columns, $sample_table->Rows()],
#		
#		'target' => unpdl( $self->{'PDL'}->slice(":,".$self->{'samples_in_result'}->{$main_cell})   ),
#		'add' =>  unpdl( $self->{'PDL'}->slice(":,(".join(",",map{$self->{'samples_in_result'}->{$_}} @cells).")") ),
#		
#		'cells' => [$main_cell, @cells],
#		'cells_using_results_location' => [@{$result->{'header'}}[map{$self->{'samples_in_result'}->{$_}}$main_cell, @cells] ],
#		'PDL_a' => ":,".$self->{'samples_in_result'}->{$main_cell},
#		'PDL_b'=> ":,(".join(",",map{$self->{'samples_in_result'}->{$_}} @cells).")",
#	} ).";\n"; 
	$self->{'PDL'}->slice(":,".$self->{'samples_in_result'}->{$main_cell}) += 
		$self->{'PDL'}->slice(":,".join(":",map{$self->{'samples_in_result'}->{$_}} @cells));
		
	$self->{'PDL'}->slice(":,".$self->{'samples_in_result'}->{$main_cell. " spliced"}) += 
		$self->{'PDL'}->slice(":,".join(":",map{$self->{'samples_in_result'}->{$_. " spliced"}} @cells));

	map { $result->rename_column( $_, "merged $_" ); $result->rename_column( $_. " spliced", "merged $_ spliced" ) } @cells;

}

sub add_to_summary {
	my ( $sampleID, $geneIDs, $is_spliced ) = @_;
	$is_spliced = 1 if ($is_spliced);

	#$sampleID = $self->{'cell_id_to_name'}->{$sampleID} || $sampleID;
	my @col_number =
	  $result->Add_2_Header( [ $sampleID, $sampleID . " spliced" ] );
	$geneIDs = [$geneIDs] unless ( ref($geneIDs) eq "ARRAY" );
	my ( @row_numbers, $row );
	foreach my $gene_id (@$geneIDs) {
		@row_numbers =
		  $result->get_rowNumbers_4_columnName_and_Entry( 'Gene_ID', $gene_id );
		if ( @row_numbers == 0 ) {    ## gene never occured
			    #print "new gene $gene_id detected for sample $sampleID\n";
			my $hash = { 'Gene_ID' => $gene_id, $sampleID => 1 };
			$hash->{ $sampleID . " spliced" } = 1 if ($is_spliced);
			$result->AddDataset($hash);
		}
		else {
			foreach my $row (@row_numbers) {
				@{ $result->{'data'} }[$row] ||= [];
				if ($is_spliced) {
					@{ @{ $result->{'data'} }[$row] }
					  [ @col_number[ 0, $is_spliced ] ]++;
				}
				else {
					@{ @{ $result->{'data'} }[$row] }[ $col_number[0] ]++;
				}
			}
		}
	}

}

sub lintersect {
	my $h;
	map { $h->{$_}++ } @_;
	my @ret;
	map { push( @ret, $_ ) if ( $h->{$_} > 1 ) } @_;
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
	return
	  map { @{ @{ $gtf->{'data'} }[$_] }[ $gtf->Header_Position('end') ]; }
	  @lines;
}

sub get_reporter_ids {
	my ( $gtf, @lines ) = @_;
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
