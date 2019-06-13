#!/usr/bin/perl
use strict; 
use warnings;

use stefans_libs::file_readers::gtf_file;

#NB501227:125:HTL2KBGX3:3:21403:14386:11684:S_TGCGCAGC_C_CATCAAGGTGCAGTAG_GCAGTATGAG	210175	32M15125N15M110462N50M
#NB501227:125:HTL2KBGX3:3:13601:5186:16519:S_CATCACTT_C_GATGAGGCAACACGCC_TACTTACAGA	222414	6S14M55166N77M
sub sample_and_umi_my {
	my (@bam_line) = @_;
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
	$UMI = $matching_IDs[4];
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
			push ( @r, "$start-".($start + $1) );
			$start += $1 -1;
		}elsif( $_ =~ m/(\d+)[NS]/ ) {
			$start += $1 -1;
		}else {
			die "I do not know what do do with that cigar part: $_\n";
		}
	}
	return @r;
}

my ( @line,$sample_name, $UMI, $seen, $duplicate , $result );
while ( <> ){
	chomp;
	@line = split("\t", $_ );
	($sample_name, $UMI) = &sample_and_umi_my(@line);
	if ( $seen->{ "$sample_name, $UMI"} ){
		$duplicate ++;
		next;
	}
	$seen->{ "$sample_name, $UMI"} =1;
	map { $result->{$_}->{$sample_name} ++ } &result_areas( @line[1,2] );
}

sub unique {
	my $h;
	map { $h->{$_}++ } @_;
	return sort keys %$h;
}

my $gtf_obj = stefans_libs::file_readers::gtf_file->new();
$gtf_obj->read_file('Spliced_Reads.gtf');
my $chr = 'chrJH792828.1';
my $rep_id = $gtf_obj->Header_Position( 'gene_id' );
my @gene_ids;
foreach ( sort keys %$result ) {
	@gene_ids = &unique(
		map {@{ @{ $gtf_obj->{'data'} }[$_] } [ $rep_id ] } 
		  $gtf_obj->efficient_match_chr_position($chr, split( "-", $_ ) )  
		);
	print join("\t", $_, scalar(keys %{$result->{$_}} ), join(" ",sort keys %{$result->{$_}}  ), join(" ",@gene_ids))."\n";
} 
