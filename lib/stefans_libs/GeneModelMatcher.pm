package stefans_libs::GeneModelMatcher;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
#created by bib_create.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
use strict;
use warnings;

use stefans_libs::GeneModelMatcher::geneModel;
use stefans_libs::GeneModelMatcher::matchSave;
use stefans_libs::file_readers::gtf_file; ## to use the efficient match without re-coding

=head1 LICENCE

  Copyright (C) 2018-11-28 Stefan Lang

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


=for comment

This document is in Pod format.  To read this, use a Pod formatter,
like 'perldoc perlpod'.

=head1 NAME

stefans_libs::gtf_file::reader

=head1 DESCRIPTION

use the gtf_file class to read a list of gene models instead of simple exon lists and make them available for matching cigar strings

=head2 depends on


=cut


=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::gtf_file::reader.
All entries of the hash will be copied into the objects hash - be careful to use that right!

=cut

sub new{

	my ( $class, $hash ) = @_;

	my ( $self );

	$self = {
		'__gname2id__' => {},
		'gtf_file' => undef,
		'data' => [],
		'id' => 'gene_id',
		'collect_over' => 'gene',
		'collect' => 'exon',
		'match_save' => stefans_libs::GeneModelMatcher::matchSave->new(),
  	};
  	foreach ( keys %{$hash} ) {
  		$self-> {$_} = $hash->{$_};
  	}

  	bless $self, $class  if ( $class eq "stefans_libs::GeneModelMatcher" );
	
	$self->check();
	if ( defined $self->{'filename'} ) {
		$self->read_file( $self->{'filename'} );
	}
	
  	return $self;

}

sub Info {
	my ( $self ) = @_;
	return map{ $_->Summary() } @{$self->{'data'}};
}

sub check{
	my ( $self ) = @_;
	my $err = "";
	foreach ( qw(id) ) {
		$err .= "key '$_' is missing in the object new(\$hash)\n" unless ( defined $self->{$_} );
	}
	Carp::confess( $err ) if ( $err =~ m/\w/ );
	return $self;
}

sub read_file {
	my ( $self, $filename, $chr ) = @_;
	$self->{'filename'} = $filename if ( -f $filename); ## open self->filename is set otherwiese ;-)
	
	my ($f, $collected, $start, $end);
		
	if ( defined $chr ) {
		if ( $chr =~ m/^([\w\\._]+)\:(\d+)-(\d+)/ ) {
			$chr = $1; $start = $2; $end = $3;
		}
	}
	$collected = 0;
	$self->{'gtf_file'} = stefans_libs::file_readers::gtf_file->new();
	@{$self->{'gtf_file'}->{'header'}}[8] = $self->{'id'};
	delete($self->{'gtf_file'}->{'header_position'}->{'annotation'});
	$self->{'gtf_file'}->{'header_position'}->{ $self->{'id'} } = 8;
	
	open ( $f, "<$self->{'filename'}" ) or Carp::confess ("Sorry I could not open the file '$self->{'filename'}'\n$!\n");
	my $line = 0;
	while( <$f> ) {
		$line ++;
		next if ( $_ =~ m/^#/ );
		chomp($_);
		my @tmp = split(/\t/, $_);
		if ( defined $chr ) {
			next unless ( $tmp[0] eq $chr );
			if ( defined $start ) {
				next unless ( $start < $tmp[4] and $end > $tmp[3] );
				#warn "Line $line made it through the $chr:$start-$end selection:\n$_\n";
			}
		}
		if ( @tmp != 9 ) {
			Carp::confess("not a gtf or gff3 line: $_\n" );
		}
		#chr1	HAVANA	gene	3073253	3074322	.	+	. "some more annotation processed by split_annotation"
		my $anno = $self->split_annotation( pop(@tmp) ); ## @tmp looses the annotation 
		if ( $tmp[2] eq $self->{'collect_over'}) {
			## that is what we need to keep here! All exons would go into the gene model anyhow.
			if ( defined $self->{'__gname2id__'}->{$anno->{$self->{'id'}}} ){
				Carp::confess ( "the same gene twice here: $_\n");
			}
			$self->{'__gname2id__'}->{$anno->{$self->{'id'}}} = scalar(@{$self->{'data'}} );
			push (@{$self->{'data'}}, stefans_libs::GeneModelMatcher::geneModel->new( @tmp, $anno ));
			push( @{$self->{'gtf_file'}->{'data'}}, [@tmp, $anno->{$self->{'id'}} ] );
		}elsif( $tmp[2] eq $self->{'collect'} ) {
			Carp::confess( "I have no '$self->{'collect_over'}' entry for this '$self->{'collect'}' feature on line $line\n". join("\t", @tmp)."\n$_\n")
				unless ( defined  $self->{'__gname2id__'}->{$anno->{$self->{'id'}}} );
			@{$self->{'data'}}[ $self->{'__gname2id__'}->{$anno->{$self->{'id'}}} ]->Add( @tmp, $anno );
			$collected++;
		}
	}
	print scalar(@{$self->{'data'}}), " $self->{'collect_over'}(s) collected storing a total of $collected $self->{'collect'}(s)\n";
	close ( $f );
	return $self;
}

=head3 genes_at_position_plus_one( $chr, $start, $end, $add )

uses gtf_file::efficient_match_chr_position_plus_one to select all matching features + the next matching feature.
Returns an array of stefans_libs::GeneModelMatcher::geneModel objects.

=cut

sub genes_at_position_plus_one {
	my ( $self, $chr, $start, $end, $add ) = @_;
	$add ||= 0;
	my @ids = $self->{'gtf_file'}->efficient_match_chr_position_plus_one( $chr, $start, $end, $add );
	return @{$self->{'data'}}[@ids];
}


=head3 match_cigar ( $chr, $start, $cigar)

Matches the cigar information of a sam/bam file to the gene models - please use sorted bam/sam files as 
the order speeds up this match enormousely!

returns the id of the matching genes id if any and an 'exon' match, 'spliced' match or 'primary' transcript match as hash.

=cut	

sub match_cigar{
	my ( $self, $chr, $start, $cigar) = @_;
	## the cigar looks like 3S8M21445N87M
	if ( $cigar =~ s/^(\d)+S// ) {
		#warn "there was a $1 bp mismatch in the start of the read - adjust $start";
		$start += $1;
	}
	my $end = $start;
	map { $end += $_ } split( /[A-Z]/, $cigar);

	unless ( $self->{'match_save'}->match( $chr, $start, $end) ){
		my @matches = $self->genes_at_position_plus_one( $chr, $start, $end  );
		#warn "rematch to genes $chr, $start, $end did produce ".scalar(@matches)." matching entries";
		$self->{'match_save'}->Update( $self->{'id'}, $chr, $start, $end, @matches );
	}
	#warn "matching the cigar $chr, $start, $cigar\n";
	return $self->{'match_save'}->match_cigar( $chr, $start, $cigar );
}


sub split_annotation {
	my ( $self, $anno ) = @_;
	my ( $ret );
	if ( $self->{'filename'} =~ m/gtf$/ ){
		#gene_id "ENSMUSG00000102693.1"; gene_type "TEC"; gene_name "4933401J01Rik"; level 2; havana_gene "OTTMUSG00000049935.1";
		$anno =~ s/;$//; 
		$anno =~ s/"//g; 
		$ret = { map{ 
			if ( $_ =~ m/([\w_]+) (.+)/) { 
				$1 => $2 
			} else { 
				Carp::confess( "wrong gft annotation format on line". scalar(@{$self->{'data'}}).": $_\n" )
			} } 
			split( /; / , $anno) 
		};		
	}elsif ( $self->{'filename'} =~ m/gff3$/ ){
		#ID=ENSMUSG00000102693.1;gene_id=ENSMUSG00000102693.1;gene_type=TEC;gene_name=4933401J01Rik;level=2;havana_gene=OTTMUSG00000049935.1
		$ret = { map{ 
			if ( $_ =~ m/([\w_]+)=(.+)/) { 
				$1 => $2 
			} else { 
				Carp::confess( "wrong gff3 annotation format on line". scalar(@{$self->{'data'}}).": $_\n" )
			} } 
			split( /;/ , $anno) 
		};
		
	}else {
		Carp::confess( "Only gtf anf gff3 files are supported here, not: '$self->{'filename'}'\n" );
	}
	return $ret;
}
1;
