#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-06-26 Stefan Lang

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

    MassQuantifyBamFiles.pl
       -infiles   :the result files from ReshuffleBamFiles.pl
       -gtf_file  :the GTF file including the regions to quantify 
       -outfile   :the final outfile
       -options   : format: key_1 value_1 key_2 value_2 ... key_n value_n

           quantify_on: which feature to select for quantification from the gtf file (default 'exon')
             report_on: which column in the gtf file object to report on (default 'gene_id')
              min_UMIs: how many UMIs have to mapp to any gene to report the sample 
                        (for one NextSeq run the default value default 100 proved usable)
                        
       -slurmOptions: format: key_1 value_1 key_2 value_2 ... key_n value_n
       

       -fastqPath :The path where the original fastq files were converted to the here used input files
       -sampleID  :The sampleID of the 10X sample to process

       -bugfix    :turn on additional bugfix in the QuantifyBamFiles runs

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Quantify results from ReshuffleBamFiles.pl to speed the processing of large files up.

  To get further help use 'MassQuantifyBamFiles.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @infiles, $gtf_file, @slurmOptions, $outfile, $options, @options, $fastqPath, $sampleID, $bugfix);

Getopt::Long::GetOptions(
       "-infiles=s{,}"    => \@infiles,
	 "-gtf_file=s"    => \$gtf_file,
	 "-outfile=s"    => \$outfile,
       "-options=s{,}"    => \@options,
       "-slurmOptions=s{,}"    => \@slurmOptions,
       
	 "-fastqPath=s"    => \$fastqPath,
	 "-sampleID=s"    => \$sampleID,
       "-bugfix"    => \$bugfix,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -f $infiles[0]) {
	$error .= "the cmd line switch -infiles is undefined!\n";
}
unless ( -f $gtf_file) {
	$error .= "the cmd line switch -gtf_file is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $slurmOptions[0]) {
	$warn .= "the cmd line switch -slurmOptions is undefined!\n";
}
unless ( -d $fastqPath) {
	$error .= "the cmd line switch -fastqPath is undefined!\n";
}
unless ( defined $sampleID) {
	$error .= "the cmd line switch -sampleID is undefined!\n";
}
# bugfix - no checks necessary


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

#$options->{'n'} ||= 10;

###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/MassQuantifyBamFiles.pl';
$task_description .= ' -infiles "'.join( '" "', @infiles ).'"' if ( defined $infiles[0]);
$task_description .= " -gtf_file '$gtf_file'" if (defined $gtf_file);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= ' -slurmOptions "'.join( '" "', @slurmOptions ).'"' if ( defined $slurmOptions[0]);

$task_description .= " -fastqPath '$fastqPath'" if (defined $fastqPath);
$task_description .= " -sampleID '$sampleID'" if (defined $sampleID);
$task_description .= " -bugfix " if ( $bugfix);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################

my $slurmOptions;
for ( my $i = 0 ; $i < @slurmOptions ; $i += 2 ) {
	$slurmOptions[ $i + 1 ] =~ s/\n/ /g;
	$slurmOptions->{ $slurmOptions[$i] } = $slurmOptions[ $i + 1 ];
}
$slurmOptions->{'A'} ||= 'lsens2017-3-2';
$slurmOptions->{'n'} ||= 1;
$slurmOptions->{'N'} = 1;
$slurmOptions->{'t'} ||= '02:00:00';
$slurmOptions->{'p'} ||= 'dell';

use stefans_libs::Version;

use stefans_libs::SLURM;
my $V = stefans_libs::Version->new();
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG '#library version '.$V->version( 'Chromium_SingleCell_Perl' )."\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $worker = stefans_libs::SLURM->new($slurmOptions, 0);
$worker->{'debug'} = $debug;

my (@tmp, $ofile, $i);
mkdir ( File::Spec->catfile($fm->{'path'}, 'chrTmp' ) );
$i = 0;
foreach my $file (sort byFileSize @infiles ) {
	$i ++;
	$fm = root->filemap( $file);
	@tmp = split("_", $fm->{'filename'});
	$ofile = File::Spec->catfile($fm->{'path'}, 'chrTmp', "part_$i" );
 	my $cmd = 'QuantifyBamFile.pl' . ' -infile "'.$file.'"'
  . " -gtf_file '$gtf_file'"
  . " -outfile '$ofile'"
  . " -drop_chr $tmp[0]"
  . ' -options "'.join( '" "', @options ).'"' 
  . " -fastqPath '$fastqPath'" 
  . " -sampleID '$sampleID'";
  
   $cmd .= " -bugfix "  if ( $bugfix );
   
   $worker->run( $cmd, $ofile, File::Spec->catfile($ofile, 'matrix.mtx' ) )
   
}

print "All processes started - hopefully\n";


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

