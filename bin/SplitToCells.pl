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


my ( $help, $debug, $database, $R1, $R2, $I1, $outpath, $options, @options);

Getopt::Long::GetOptions(
	 "-R1=s"    => \$R1,
	 "-R2=s"    => \$R2,
	 "-I1=s"    => \$I1,
	 "-outpath=s"    => \$outpath,
       "-options=s{,}"    => \@options,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $R1) {
	$error .= "the cmd line switch -R1 is undefined!\n";
}
unless ( defined $R2) {
	$error .= "the cmd line switch -R2 is undefined!\n";
}
unless ( defined $I1) {
	$error .= "the cmd line switch -I1 is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $options[0]) {
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

### initialize default options:

$options->{'cell_barcode_length'} ||= 16;
$options->{'oname'} ||= 'projectX';
###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/SplitToCells.pl';
$task_description .= " -R1 '$R1'" if (defined $R1);
$task_description .= " -R2 '$R2'" if (defined $R2);
$task_description .= " -I1 '$I1'" if (defined $I1);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################

mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/$outpath/$options->{'oname'}.annotated.fastq.gz".$$."_SplitToCells.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


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

my $ofile = "$outpath/$options->{'oname'}.annotated.fastq.gz";
open ( my $OUT ,"| /bin/gzip -c > $ofile" ) or die "I could not open the out pipe '| /bin/gzip -c > $ofile'\n$!\n";


my $func = sub{
	my ($fastqfile, @entries) = @_;## f1, f2, i1
	
	my $cellid = "S_". $entries[2]->sequence()."_C_". substr($entries[1]->sequence,0,$options->{'cell_barcode_length'});
	my $UMI_tag = substr($entries[1]->sequence,$options->{'cell_barcode_length'});
	$counter->{$cellid} ++;
	$entries[0]->Add_UMI_Tag($cellid."_".$UMI_tag);
	$entries[0]->write( $OUT );	
	
	for ( my $i = 0 ; $i < @entries ; $i++ ) {
		$entries[$i]->clear();
	}
};

my $worker = stefans_libs::FastqFile->new();


$worker->filter_multiple_files(
	$func, $R2, $R1, $I1
);

close ( $OUT );


open ( OUT ,">$outpath/$options->{'oname'}.per_cell_read_count.xls" ) or die "I could not open the file '$outpath/$options->{'oname'}.per_cell_read_count.xls'\n$!\n";
print OUT "Cell_ID\tcount\n";
foreach my $key ( sort  keys %$counter ) {
	print OUT "$key\t$counter->{$key}\n";
}
close ( OUT );
print "A detailed per cell read count has been written to '$outpath/$options->{'oname'}.per_cell_read_count.xls'\n";

