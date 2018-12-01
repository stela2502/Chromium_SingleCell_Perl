package stefans_libs::GeneModelMatcher::matchSave;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
#created by bib_create.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
use strict;
use warnings;

=head1 LICENCE

  Copyright (C) 2018-11-29 Stefan Lang

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

stefans_libs::GeneModelMatcher::matchSave

=head1 DESCRIPTION

stores information to the last match

=head2 depends on


=cut


=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::GeneModelMatcher::matchSave.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new{

	my ( $class, $hash ) = @_;

	my ( $self );

	$self = {
		'skip_rest' => "GunsetChr",
		'this_region_start' => 0,
  	};
  	foreach ( keys %{$hash} ) {
  		$self-> {$_} = $hash->{$_};
  	}

  	bless $self, $class  if ( $class eq "stefans_libs::GeneModelMatcher::matchSave" );

  	return $self;

}

=head3 match( $chr, $start, $end )

Checks if the region is still covered by the match

returns 1 => OK 0 => no

=cut

sub match{
	my ( $self, $chr, $start, $end ) = @_;
	return 1 if ( $self->{'skip_rest'} eq $chr );
	return 0 unless ( defined $self->{'next'});
	return ($self->{'next'}->{'start'} > $end and $self->{'chr'} eq $chr);	
}
=head3 Update( $nameTag, $chr, $start, $end, @geneModels )

re-inits all and stores the information

returns 1 => OK 0 => no

=cut
sub Update {
	my ( $self, $nameTag, $chr, $start, $end, @geneModels ) = @_;
	if ( scalar(@geneModels) == 0 ){
		## OK rest of the chromosome is not informative - Match to all other queries and always return a empty results hash!
		$self->{'skip_rest'} = $chr;
		return $self;
	}
	$self->{'nameTag'} = $nameTag;
	$self->{'chr'} = $chr;
	$self->{'next'} = pop(@geneModels);
	$self->{'models'} = [ @geneModels ];
	if ( scalar(@geneModels) == 0 ) {
		## that can meen two things: no match here or last match in chr
		if ( $self->{'next'}->total_match( $start, $end) ) {
			#warn "no next entry!?\n";
			$self->{'models'} = [ $self->{'next'} ];
			## make up the next one - will not be used anyhow
			$self->{'next'} = 
			   stefans_libs::GeneModelMatcher::geneModel->new(
			   $chr, "unused", 
			   "gene", @{$self->{'models'}}[0]->{'end'} + 10 , @{$self->{'models'}}[0]->{'end'} +33, 
				".",	"+", ".", {$self->{'nameTag'} => 'Fake' });
		}else {
			#warn "wanted $chr:$start-$end length ". ($end - $start)." vs match ".$self->{'next'}->Summary()." Mean distance ".( ($end + $start)/2 -($self->{'next'}->{'start'}+ $self->{'next'} ->{'end'} ) /2 )." No matching gene - wait for the next...\n";		
		}
		## should lead to a re-mapping after the actual gene
	}
	if ( scalar(  @{$self->{'models'}} ) > 0){
		$self->{'this_region_start'} = @{$self->{'models'}}[0]->{'start'}
	}
	else {
		$self->{'this_region_start'} =  $self->{'next'}->{'start'}
	}
	
	$self->{'next_model_change'} = 0;
	$self->{'last_match'} = {};
	return $self;
}

=head3 match_cigar ( $chr, $start, $cigar )

Processes the cigar string and checks if the read matches to any know gene, any exon or is a spliced exon matching read.
This function expects, that the object already contains all optional genes to ever match this read.

returns the id of the matching genes and if any an 'exon' match, 'spliced' match or 'primary' transcript match as hash.

Say genes G1 and G2 are available as they both span this area; G1 has an primary match and G2 has an spliced exon match this would retiurn:

$ret = {
	'G1' => 'primary',
	'G2' => 'spliced'
}
	 

=cut

sub match_cigar{
	my ( $self, $chr, $start, $cigar ) = @_;
	my $saveC;
	# cigar like 8M21445N87M
	my $ret = {};
	$self->{'slast_start'} ||= $start;
	Carp::confess( "This function requres to get ordered start positions and the last wone was after this one.\nThis will lead to errors!\n" )
		if ( $self->{'slast_start'} > $start);
	if ( $self->{'skip_rest'} eq $chr ){
		#warn "No more genes in chr $chr:$start and $cigar\n";
		return $ret;
	}
	my $end = $start;
	#warn "I should match the values $chr, $start, $cigar to my ".scalar(@{$self->{'models'}})." gene models";
	if ( $cigar =~ s/^(\d+)M//) {
		$end += $1;
	}
	$saveC = "$chr:$start-$end";
	#warn "I have set the start to $start and end to $end\n";
	if ( $start < $self->{'next_model_change'} ) {
		#warn ref($self).": I use the old match ($start < $self->{'next_model_change'})\n" if (1);
		$ret = $self->{'last_match'};
	}
	else {
		#warn "matching first part: $chr:$start-$end ( I have ".scalar(@{$self->{'models'}})." gene models to match with)\n";
		foreach my $geneModel ( @{$self->{'models'}}) {
			if (  $geneModel->{'start'} <  $start and $geneModel->{'end'} >= $end ) {
				#warn "I check the model $geneModel->{'info'}->{$self->{'nameTag'}}\n";
				if ( $geneModel->match( $start, $end, 10 ) ){
					#warn "ret $geneModel->{'info'}->{$self->{'nameTag'}} set to 'exon'";
					$ret->{$geneModel->{'info'}->{$self->{'nameTag'}}} = 'exon';
				}else {
					#warn "ret $geneModel->{'info'}->{$self->{'nameTag'}} set to 'primary'";
					$ret->{$geneModel->{'info'}->{$self->{'nameTag'}}} = 'primary';
				}
				my $tmp = $geneModel->NextChange( $start );
				$self->{'next_model_change'} = $tmp if ( $tmp < $self->{'next_model_change'} or $self->{'next_model_change'} ==0);
			}
		}
		$self->{'last_match'} = {map { $_ => $ret->{$_}} keys %$ret}; ## copy this hash as that could be modified later.
	}
	if ( $cigar =~ m/\d/ ) {
		#warn ( "OK there might be more to this - I could find a spliced hit" );
		$start = $end;
		if ( $cigar =~ s/(\d+)N// ){
			$start += $1 ;
		}
		$end = $start;
		if ($cigar =~ s/(\d+)M// ){
			$end += $1;
		}
		#warn "matching second part of cigar: $chr:$start-$end having initailly matched to the genes '".join(" ", keys %$ret)."'\n";
		foreach my $geneModel ( @{$self->{'models'}} ) {
			next unless ( defined $ret->{$geneModel->{'info'}->{$self->{'nameTag'}}});
			if (  $geneModel->{'start'} <  $start and $geneModel->{'end'} >= $end ) {
				if ( $geneModel->match( $start, $end ) ){
					$ret->{$geneModel->{'info'}->{$self->{'nameTag'}}} = 'spliced';
				}else {
					$ret->{$geneModel->{'info'}->{$self->{'nameTag'}}} = 'primary';
				}
			}elsif ( defined $ret->{$geneModel->{'info'}->{$self->{'nameTag'}}} ) {
				#warn "no match for $chr:$start-$end to gene $geneModel->{'info'}->{$self->{'nameTag'}}: ret $geneModel->{'info'}->{$self->{'nameTag'}} removed";
				delete( $ret->{$geneModel->{'info'}->{$self->{'nameTag'}}}  ) ;
			}
		}
		$saveC .= " + $chr:$start-$end";
	}
	#warn "$chr, $start, $saveC: return matches for ".scalar( keys %$ret)." genes :". join(" ", keys %$ret )."\n";
	return $ret;
}
1;
