#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-10-30 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
   

=head1  SYNOPSIS

    CellRangerAggregate.pl
       -path       :<please add some info!>
       -options     :<please add some info!> you can specify more entries to that
                         format: key_1 value_1 key_2 value_2 ... key_n value_n
       -outpath       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  start an agrgregate run

  To get further help use 'CellRangerAggregate.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use File::Spec;

use stefans_libs::SLURM;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $path, $options, @options, $outpath, $summary_ID );

Getopt::Long::GetOptions(
	 "-path=s"    => \$path,
     "-options=s{,}"    => \@options,
	 "-outpath=s"    => \$outpath,
	 "-summary_ID=s" => \$summary_ID,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $path) {
	$error .= "the cmd line switch -path is undefined!\n";
}
unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $summary_ID) {
	$error .= "the cmd line switch -summary_ID is undefined!\n";
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

#$options->{'n'} ||= 10;

###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/CellRangerAggregate.pl';
$task_description .= " -path '$path'" if (defined $path);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -outpath '$outpath'" if (defined $outpath);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################
mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_CellRangerAggregate.pl.log") or die $!;
print LOG $task_description."\n";



## Do whatever you want!

##
# first find all <sampleName>/outs/molecule_info.h5 files used the sample name to create the Aggregation CSV
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
# and I use find here...

open ( IN, "find $path -name 'molecule_info.h5' |" ) or die "could not start the 'find' fork\n$!\n";
open (OUT, ">$outpath/Aggregation.csv" ) or die "I could not create the file '$outpath/Aggregation.csv'\n$!\n";
print OUT "library_id,molecule_h5\n";
my @tmp;
while ( <IN> ) {
	chomp();
	
	@tmp = split( "/", $_ );
	print OUT "$tmp[-3],".File::Spec->rel2abs($_)."\n";
}

close ( IN );
close ( OUT );

## I will use the runCommand.pl script to actually run the command.
## This allows to adapt to different settings. e.g. runCommand uses SLURM or some other method 
## or it would simply create a bash script.

my $cmd = "cd $path\ncellranger aggr --id=$summary_ID --csv=$outpath/Aggregation.csv --normalize=none\n";
print "CellRanger cmd:\t$cmd\n";
print LOG "#CellRanger cmd:\t$cmd\n";
my $runner = "runCommand.pl -options ".join(" ", @options )
. " -loadModules 'GCCcore/6.3.0' 'bcl2fastq/2.19.1' 'cellranger/2.2.0'"
. " -cmd '$cmd'"
. " -outfile $path/CellRangerAggregate"
;
if ( $debug ) {
	$runner .= " -debug";
}

print LOG "#run:\t$runner\n";

close ( LOG );

print $runner."\n";

system( $runner );

