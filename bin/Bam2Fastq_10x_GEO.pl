#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-03-26 Stefan Lang

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
   
   binCreate.pl from git@gitlab:stefanlang/Stefans_Lib_Esentials.git commit e7545a27e262074f058adf34953b79b184a19248
   

=head1  SYNOPSIS

    Bam2Fastq_10x_GEO.pl
       -infiless       :a input SRA or bam file downloaded from NCBI GEO containing 10x data
       -outpath        :the outpath

       -options:       A list of 'key' 'value' combinations
       
           sname  :the prefix for the outfile

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Convert a 10x sra/bam file downloaded from the NCBI GEO archive into a fastq file rewriting the sequence name to be compatible with QuantifBamFile.pl.

  To get further help use 'Bam2Fastq_10x_GEO.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::BAMfile;
use stefans_libs::FastqFile::FastqEntry;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @infiles, @options, $outpath, $options);

Getopt::Long::GetOptions(
	 "-infiles=s{,}"    => \@infiles,
	 "-outpath=s"       => \$outpath,
	  "-options=s{,}"   => \@options,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -f $infiles[0]) {
	$error .= "the cmd line switch -infiles is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}

unless ( defined $options[0] ) {
	$warn .= "the cmd line switch -options is undefined!\n";
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
	print "$errorMessage.\n";
	pod2usage(q(-verbose) => 1);
}


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/Bam2Fastq_10x_GEO.pl';
$task_description .= " -infiless '".join("' '",@infiles)."'" if (defined $infiles[0]);
$task_description .= " -outfile '$outpath'" if (defined $outpath);
$task_description .= ' -options \'' . join("' '", (%$options) ) . "'";


use stefans_libs::Version;
my $V = stefans_libs::Version->new();



mkdir($outpath) unless ( -d $outpath );
open( LOG,
	    ">$outpath/$options->{'oname'}.annotated.fastq.gz"
	  . $$
	  . "_SplitToCells.pl.log" )
  or die $!;
print LOG '#library version '.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG $task_description . "\n";
close(LOG);


## Do whatever you want!

my ( $counter, $OUT, $filtered_out, $total_reads);

sub sample_and_umi_cellranger {
	my (@bam_line) = @_;
	my ( $sample_name, $UMI );
	#print join("\t", @bam_line )."\n";
	foreach (@bam_line) {
		if (scalar( split(":", $_)) == 3 ) {
			#print "\t$_";
		}else {
			next;
		}
		#$sample_name = $1 if ( $_ =~m/CB:Z:([ACTGN]+)-?\d*$/);
		if ( $_ =~ m/C[BR]:Z:([AGCTN]+)-?\d*$/ ){
			#print "Matching to sample $_\n";
			$sample_name = $1 ;
		}
		#$sample_name = $1 if ( $_ =~ m/CB:Z:([AGCTN]+)$/ );
		if ( $_ =~ m/U[BR]:Z:([ACTGN]+)$/ and ! defined $UMI){
			#print "Matching UMI $_\n";
			$UMI  = $1 ;
		}
		
	}
	#print "\n";
	return ( $sample_name, $UMI );
}
my ( $sample_name, $UMI, $entry );

#	my @matching_IDs = split( ":", $bam_line[0] );
#	@matching_IDs = split( "_", pop(@matching_IDs) );
#
#	if ( ( !defined $matching_IDs[2] ) or ( !defined $matching_IDs[3] ) ){
#		$sample_name = $matching_IDs[3];
#		$UMI = $matching_IDs[4];
#	}

sub filter_read {
	my ( $read, $min_length ) = @_;

	$min_length ||= 10;

	## filter polyT at read end
	if ( $read->sequence() =~ m/([Aa]{5}[Aa]+)$/ ) {
		$read->trim( 'end', length( $read->sequence() ) - length($1) );
	}
	if ( $read->sequence() =~ m/([Tt]{5}[Tt]+)$/ ) {
		$read->trim( 'end', length( $read->sequence() ) - length($1) );
	}
		
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

	## return filtered read

	return $read;
}

sub filter {
	shift; ## get rid of the BAMfile object
	my @bam_line = split( "\t", shift );
	return if ( $bam_line[0] =~ m/^\@/ );
	$total_reads ++;
	if ( $total_reads % 1e+5 == 0 ) {
		print ".";
	}
	( $sample_name, $UMI ) = &sample_and_umi_cellranger( @bam_line );
	
	if ( defined $UMI and defined $sample_name ) {
		#print "process next line and got '$sample_name' and '$UMI'\n";
		$entry->clear();
		$entry->name ( '@'."$bam_line[0]:".join("_", 'C','ACGT','C',$sample_name,$UMI));
		$entry->sequence($bam_line[9]);
		$entry->quality($bam_line[10]);
		$entry = filter_read($entry);
		if ( length( $entry->sequence() ) > 20 ) {
			$entry->write($OUT);
			$counter->{$sample_name}++;
		}else {
			$filtered_out++;
		}
	}
	else {
		#print "No UMI and sampleID\n";
		$filtered_out++;
	}

};

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

if ( -t STDOUT ) {
	$| = 1;
}


my ($ofile);

if ( -t STDOUT ) { ## for the progress bar
	$| = 1;
}

$entry = stefans_libs::FastqFile::FastqEntry->new();
print "Started a summary run on ".scalar(@infiles)." file sets\n";
	
	$ofile = "$outpath/$options->{'oname'}.annotated.fastq.gz";
	open( $OUT, "| /bin/gzip -c > $ofile" )
	  or die "I could not open the out pipe '| /bin/gzip -c > $ofile'\n$!\n";

	my $bam_file = stefans_libs::BAMfile->new();
	$bam_file ->{'debug'} = 1 if ( $debug ); #restrict output to first 1000 reads

	for ( my $i = 0 ; $i < @infiles ; $i++ ) {
		print "Process file $i '$infiles[$i]'\n";	
		$bam_file->filter_file( $infiles[$i], \&filter );

	}

	close($OUT);

	open( OUT, ">$outpath/$options->{'oname'}.per_cell_read_count.xls" )
	  or die
"I could not open the file '$outpath/$options->{'oname'}.per_cell_read_count.xls'\n$!\n";
	print OUT "#total of $total_reads reads processed\n";
	print OUT "#no sample/UMI information in $filtered_out reads\n";
	print OUT "Cell_ID\tcount\n";
	foreach my $key ( sort keys %$counter ) {
		print OUT "$key\t$counter->{$key}\n";
	}
	close(OUT);
	print "Finished with fastq file '$ofile'\n";
	print "A total of $total_reads reads were processed\n";
	print "and $filtered_out reads had no sample/UMI information and were filtered out\n";
	print sprintf("%.3f",(($total_reads - $filtered_out)/ $total_reads ) * 100 )."% reads could be used\n";
	

