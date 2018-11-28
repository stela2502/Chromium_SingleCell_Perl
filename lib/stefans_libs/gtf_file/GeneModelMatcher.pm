package stefans_libs::gtf_file::GeneModelMatcher;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
#created by bib_create.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
use strict;
use warnings;

use stefans_libs::gtf_file::geneModel;

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
		'data' => [],
		'id' => 'gene_id',
		'collect_over' => 'gene',
		'collect' => 'exon',
  	};
  	foreach ( keys %{$hash} ) {
  		$self-> {$_} = $hash->{$_};
  	}

  	bless $self, $class  if ( $class eq "stefans_libs::gtf_file::GeneModelMatcher" );
	
	$self->check();
	if ( defined $self->{'filename'} ) {
		$self->read_file( $self->{'filename'} );
	}
	
  	return $self;

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
	my ( $self, $filename ) = @_;
	$self->{'filename'} = $filename if ( -f $filename); ## ope self->filename is set otherwiese ;-)
	my ($f, $collected);
	$collected = 0;
	open ( $f, "<$self->{'filename'}" ) or Carp::confess ("Sorry I could not open the file '$self->{'filename'}'\n$!\n");
	while( <$f> ) {
		next if ( $_ =~ m/^#/ );
		chomp($_);
		my @tmp = split(/\t/, $_);
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
			push (@{$self->{'data'}}, stefans_libs::gtf_file::geneModel->new( @tmp, $anno ));
		}elsif( $tmp[2] eq $self->{'collect'} ) {
			@{$self->{'data'}}[ $self->{'__gname2id__'}->{$anno->{$self->{'id'}}} ]->Add( @tmp, $anno );
			$collected++;
		}
	}
	print scalar(@{$self->{'data'}}), " $self->{'collect_over'}(s) collected storing a total of $collected $self->{'collect'}(s)\n";
	close ( $f );
	return $self;
}


sub split_annotation {
	my ( $self, $anno ) = @_;
	my ( $ret );
	if ( $self->{'filename'} =~ m/gtf$/ ){
		#gene_id "ENSMUSG00000102693.1"; gene_type "TEC"; gene_name "4933401J01Rik"; level 2; havana_gene "OTTMUSG00000049935.1";
		$anno =~ s/;$//; 
		$ret = { map{ 
			if ( $_ =~ m/([\w_]+) "(.+)"/) { 
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
