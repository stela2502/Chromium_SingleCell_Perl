#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-02-17 Stefan Lang

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

    SplitToCells.pl
       -R1       :the R1 fastq file
       -R2       :the R2 fastq file
       -I1       :the I1 fastq file
       -outpath  :the outpath for the sample specific fastq files
       -split    :if you supply a list of fastq files for each R1 R2 and I1
                  this data can be summed up in one fastq file (default)
                  or split into n fastq files using this option
       
       -options  :the options in format key 'value1 value2 value3 ...'
       
               oname :the outfile prefix (default 'projectX' )
 cell_barcode_length :the length of the cell barcode (should be 16 unless changed by Illumina)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  take three fastq files and plit them into single cells

  To get further help use 'SplitToCells.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::FastqFile;

use IO::File;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help, $debug,   $database, @R1,      @R2,
	@I1,   $outpath, $options,  @options, $split
);

Getopt::Long::GetOptions(
	"-R1=s{,}"      => \@R1,
	"-R2=s{,}"      => \@R2,
	"-I1=s{,}"      => \@I1,
	"-outpath=s"    => \$outpath,
	"-options=s{,}" => \@options,
	"-split"        => \$split,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( -f $R1[0] ) {
	$error .= "the cmd line switch -R1 is undefined!\n";
}
unless ( -f $R2[0] ) {
	$error .= "the cmd line switch -R2 is undefined!\n";
}
unless ( -f $I1[0] ) {
	$error .= "the cmd line switch -I1 is undefined!\n";
}
unless ( defined $outpath ) {
	$error .= "the cmd line switch -outpath is undefined!\n";
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

### initialize default options:

$options->{'cell_barcode_length'} ||= 16;
$options->{'oname'}               ||= 'projectX';
###

my ($task_description);
for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}

my $tmp_oname;

for ( my $i = 0 ; $i < @R1; $i ++ ){
	$tmp_oname = $options->{'oname'};
	$options->{'oname'} .= ".$i";
	$task_description .= 'perl ' . $plugin_path . '/SplitToCells.pl';
$task_description .= " -R1 '$R1[$i]'";
$task_description .= " -R2 '$R2[$i]'";
$task_description .= " -I1 '$I1[$i]'";
$task_description .= " -outpath '$outpath'" if ( defined $outpath );
$task_description .= ' -options \'' . join("' '", (%$options) ) . "'"
  if ( defined $options[0] );
$task_description .= ' -split' if ($split);
$task_description .= "\n";
	 $options->{'oname'} = $tmp_oname;
}


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################

mkdir($outpath) unless ( -d $outpath );
open( LOG,
	    ">$outpath/$options->{'oname'}.annotated.fastq.gz"
	  . $$
	  . "_SplitToCells.pl.log" )
  or die $!;
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!

###Where does the sample / UMI tag come from

### Main sequence:
#R2:     GCAGTGGTATCAACGCAGAGTACATGGGGGGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGCAAATCGAGGCC
#ref_R2: GGCCTCGATTTGCAGTGATGGCCGACTAGGCCATCTTTTGATACATATGCAGCTAGAGTCAAGAGCTCCCCCCATGTACTCTGCGTTGATACCACTGC
#BAM:    GGCCTCGATTTGCAGTGATGGCCGACTAGGCCATCTTTTGATACATATGCAGCTAGAGTCAAGAGCTCCCCCCATGTACTCTGCGTTGATACCACTGC
#
#
#### CR:Z Chromium cellular barcode sequence as reported by the sequencer.
#R1:          GTAACGTCAACGATCTCGATTCCGCG
#BAM:    CR:Z:GTAACGTCAACGATCT
#
#### UR:Z        Chromium molecular barcode sequence as reported by the sequencer.
#R1:          GTAACGTCAACGATCTCGATTCCGCG
#BAM:                    UB:Z:CGATTCCGCG
#
#### BC:Z Sample index read.
#I1:          TGCTCGTA
#BAM:    BC:Z:TGCTCGTA

## convert this info into a perl function...
## that get one set of fastq entries at a time
#@NB501227:57:HGJM5BGX2:1:22212:15012:5038 1:N:0:TGCTCGTA
#GTAACGTCAACGATCTCGATTCCGCG
#+
#AAAAAEEEEEEEE/EEEEEEEEEEEE

my $ofiles;
my $counter;
my $filtered_out   = 0;
my $filtered_polyT = 0;


## turn on autoflush for the process bar:
my $flush_counter = 0;
$| = 1;


sub filter_reads {
	my ( $read, $min_length ) = @_;

	$min_length ||= 10;

	## filter polyT at read end
	if ( $read->sequence() =~ m/([Aa]{9}[Aa]+)$/ ) {
		$read->trim( 'end', length( $read->sequence() ) - length($1) );
		$filtered_polyT++;
	}
	if ( $read->sequence() =~ m/([Tt]{9}[Tt]+)$/ ) {
		$read->trim( 'end', length( $read->sequence() ) - length($1) );
		$filtered_polyT++;
	}
	
	$read->filter_low_quality(20); ## throw away crap!
	
	## filter reads with high ployX (>= 50%)

	my $str = $read->sequence();
	foreach ( 'Aa', 'Cc', 'Tt', 'Gg' ) {
		foreach my $repl ( $str =~ m/[$_]{$min_length}[$_]*/g ) {
			my $by = 'N' x length($repl);
			$str =~ s/$repl/$by/;
		}
	}
	my $count = $str =~ tr/N/n/;

	if ( $count != 0 and $count / length($str) > 0.5 ) {
		return undef;
	}

	if ( ++ $flush_counter % 10000 == 0 ) {
		print '.';
	} 
	## return filtered read

	return $read;
}
my ( $ofile, $OUT );

my $func = sub {
	my ( $fastqfile, @entries ) = @_;    ## f1, f2, i1

	my $UMI_tag = filter_reads( $entries[0] );

	if ( defined $UMI_tag and length( $UMI_tag->sequence() ) > 50 ) {
		$entries[0] = $UMI_tag;
		my $cellid = "S_"
		  . $entries[2]->sequence() . "_C_"
		  . substr( $entries[1]->sequence, 0,
			$options->{'cell_barcode_length'} );
		$UMI_tag =
		  substr( $entries[1]->sequence, $options->{'cell_barcode_length'} );
		$counter->{$cellid}++;
		$entries[0]->Add_UMI_Tag( $cellid . "_" . $UMI_tag );
		$entries[0]->write($OUT);
	}
	else {
		$filtered_out++;
	}

	for ( my $i = 0 ; $i < @entries ; $i++ ) {
		$entries[$i]->clear();
	}
};

my ( $R1, $R2, $I1 );

if ($split) {
	my $worker = stefans_libs::FastqFile->new();
	print "Started a split run on ".scalar(@R2)." file sets\n";
	for ( my $i = 0 ; $i < @R2 ; $i++ ) {
		print "Process file $i R1: '$R1[$i]'\n";
		$ofile = "$outpath/$options->{'oname'}.$i.annotated.fastq.gz";
		open( $OUT, "| /bin/gzip -c > $ofile" )
		  or die
		  "I could not open the out pipe '| /bin/gzip -c > $ofile'\n$!\n";

		( $R1, $R2, $I1 ) = ( $R1[$i], $R2[$i], $I1[$i] );
		$worker->filter_multiple_files( $func, $R2, $R1, $I1 );

		close($OUT);

		open( OUT, ">$outpath/$options->{'oname'}.$i.per_cell_read_count.xls" )
		  or die
"I could not open the file '$outpath/$options->{'oname'}.$i.per_cell_read_count.xls'\n$!\n";
		print OUT "#reads containing polyA:\t$filtered_polyT\n";
		print OUT "#filtered low complexity reads:\t$filtered_out\n";

		print OUT "Cell_ID\tcount\n";
		foreach my $key ( sort keys %$counter ) {
			print OUT "$key\t$counter->{$key}\n";
		}
		close(OUT);
		$filtered_polyT = $filtered_out = 0;
		$counter = undef;
	}
}
else {
	print "Started a summary run on ".scalar(@R2)." file sets\n";
	
	$ofile = "$outpath/$options->{'oname'}.annotated.fastq.gz";
	open( $OUT, "| /bin/gzip -c > $ofile" )
	  or die "I could not open the out pipe '| /bin/gzip -c > $ofile'\n$!\n";

	my $worker = stefans_libs::FastqFile->new();

	for ( my $i = 0 ; $i < @R2 ; $i++ ) {
		print "Process file $i R1: '$R1[$i]'\n";
		
		( $R1, $R2, $I1 ) = ( $R1[$i], $R2[$i], $I1[$i] );
		$worker->filter_multiple_files( $func, $R2, $R1, $I1 );
	}

	close($OUT);

	open( OUT, ">$outpath/$options->{'oname'}.per_cell_read_count.xls" )
	  or die
"I could not open the file '$outpath/$options->{'oname'}.per_cell_read_count.xls'\n$!\n";
	print OUT "#reads containing polyA:\t$filtered_polyT\n";
	print OUT "#filtered low complexity reads:\t$filtered_out\n";

	print OUT "Cell_ID\tcount\n";
	foreach my $key ( sort keys %$counter ) {
		print OUT "$key\t$counter->{$key}\n";
	}
	close(OUT);
}

print
"A detailed per cell read count has been written to '$outpath/$options->{'oname'}.per_cell_read_count.xls'\n";

