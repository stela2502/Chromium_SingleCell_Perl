package stefans_libs::gtf_file::geneModel;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
#created by bib_create.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
use strict;
use warnings;


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

stefans_libs::gtf_file::geneModel

=head1 DESCRIPTION

the gtf_file gene model - no transcript support at the moment

=head2 depends on


=cut


=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::gtf_file::geneModel.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new{

	my ( $class, @data ) = @_;
	
	my $hash = pop(@data);
	unless ( ref($hash) eq "HASH" ) {
		Carp::confess ( "I expect a hash as last argument - not". $hash);
	}
	unless ( @data == 8 ) {
		Carp::confess ( "I expect an 9 entry array as argument" );
	}
	my ( $self );
	#chr3	HAVANA	gene	37838468	37839573	.	+	.	$hash(ed information)
	$self = {
		'chr' => $data[0],
		'start' => $data[3],
		'end' => $data[4],
		'orientation' => $data[6],
		'info' => $hash,
		'last' => undef,
		'data' => [],
  	};

  	bless $self, $class  if ( $class eq "stefans_libs::gtf_file::geneModel" );

  	return $self;

}

sub Add{
	my ( $self, @data ) = @_;
	my $added = 0;
	foreach ( @{$self->{'data'}} ) {
		if ($_ -> match ($data[3], $data[4], 0) ) {
			$_ -> Add( $data[3], $data[4] );
			$added = 1;
			last;
		}
	}
	unless ( $added ) {
		push( @{$self->{'data'}}, this::entry->new(@data[3,4]) );
	} 
}

=head3 match( start, end )

find inside exact matches for this gene + 50 bp at the end of the gene; returns 1 if match and 0 if not

=cut

sub match {
	my ( $self, $start, $end ) = @_;
	## lets do the most likely first: only match to the last one;
	unless ( defined $self->{'last'} ){
		@{$self->{'data'}} = [ sort { $a->{'start'} <=> $b->{'start'} } @{$self->{'data'}} ];
		if ( $self->{'orientation'} eq '+' ) {
			my $id = pop(@{$self->{'data'}});
			$self->{'last'} = this::entry->new( $id->{'start'}, $id->{'end'} + 50 );
		}else {
			my $id = shift( @{$self->{'data'}} );
			$self->{'last'} = this::entry->new( $id->{'start'} - 50, $id->{'end'} );
		}
	}
	if ( $self->{'last'}->exact_match( $start, $end ) ) {
		return 1;	
	}
	for( my $i = 0; $i < @{$self->{'data'}}; $i ++ ) {
		return 1 if (@{$self->{'data'}}[$i]-> exact_match( $start, $end ) )	
	}
	return 0;
}


package this::entry;

sub new{

	my ( $class, $start, $end ) = @_;
	

	my ( $self );
	#chr3	HAVANA	gene	37838468	37839573	.	+	.	$hash(ed information)
	if ( $start < $end ){
		$self = {
			'start' => $start,
			'end' => $end
  		};
	}else {
		$self = {
			'start' => $end,
			'end' => $start
  		};
	}
  	bless $self, $class  if ( $class eq "this::entry" );
	
  	return $self;

}

sub Add {
	my ( $self, $start, $end ) = @_;
	if ( $self->{'start'} > $start ) {
		$self->{'start' } = $start;
	}
	if ( $self->{'end'} < $end ) {
		$self->{'end'} = $end; 
	}
	return $self;
}

sub match{
	my ( $self, $start, $end, $add ) = @_;
	$add ||=0;
	if ( $start <= $self->{'end'} + $add and $end >= $self->{'start'} - $add ) {
		return 1;
	}
	return 0;
}

sub exact_match{
	my ( $self, $start, $end) = @_;
	if ( $start >= $self->{'start'} and $end <= $self->{'end'} ) {
		return 1;
	}
	return 0;
}

1;
