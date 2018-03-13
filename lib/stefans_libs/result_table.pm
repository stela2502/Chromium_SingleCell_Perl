package stefans_libs::result_table;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
#created by bib_create.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c59a74e96bf8f751bb43a81e8e7cd1cef7c84747
use strict;
use warnings;

use stefans_libs::flexible_data_structures::data_table;

use base 'data_table';

use stefans_libs::database::Chromium_SingleCell::datavalues;
use stefans_libs::database::SQLiteBatch;
use File::Spec;

=head1 LICENCE

  Copyright (C) 2018-02-22 Stefan Lang

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

stefans_libs::result_table

=head1 DESCRIPTION

A specific results table, that stores the data in a sqlite database instead of a xls file. Specificly usable for single cell expression data.

=head2 depends on


=cut

=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::result_table.
All entries of the hash will be copied into the objects hash - be careful t use that right!

Add filename to already read in a database sample and gene information

=cut

sub new {

	my ( $class, $hash ) = @_;

	my ($self);
	my $tmp = $hash->{'filename'};
	$hash->{'filename'} = undef;
	$self = data_table->new($hash);
	foreach ('data_storage_spliced') {
		$self->{$_} = $hash->{$_};
	}
	$self->{'__gene2id__'}      = {};
	$self->{'__id2sample__'}    = {};
	$self->{'__sample2id__'}    = {};
	$self->{'__lastSampleID__'} = 0;
	$self->{'__lastGeneID__'}   = 0;
	$self->{'__lastDataID__'}   = 0;

	bless $self, $class if ( $class eq "stefans_libs::result_table" );

	$self->read_file($tmp);

	return $self;

}

## This is the most basic print function replace that and write_table and write_file as well as print2file will do the same thing.

=head2 print2file ( <filename>, <lines> )

The basic write to file function from data_table overloaded to store the data in a Chromium_SingleCell::datavalues sqlite database.

Be aware, that the software removes lines from the table! So please make sure, that you finished all a line before you put it out to the database.

=cut

sub info {
	my $self = shift;
	return
	    "object of class "
	  . ref($self)
	  . "\nwith "
	  . $self->Lines()
	  . " rows and "
	  . scalar( @{ $self->{'header'} } )
	  . " columns\n"
	  . "\nThe database contains $self->{'__lastSampleID__'} samples, "
	  . "$self->{'__lastGeneID__'} genes "
	  . "and $self->{'__lastDataID__'} data points\n";

}

sub gene2id {
	my ( $self, $gene_name, $id ) = @_;
	return undef unless ( defined $gene_name );
	$id ||= $self->{'__gene2id__'}->{$gene_name};
	$id = scalar( keys %{ $self->{'__gene2id__'} } ) + 1
	  unless ( defined $id );
	if ( defined $id ) {
		$self->{'__gene2id__'}->{$gene_name} = $id;
		$self->{'__id2gene__'}->{$id}        = $gene_name;
	}
	return $self->{'__gene2id__'}->{$gene_name};
}

sub id2gene {
	my ( $self, $id, $gene_name ) = @_;
	if ( defined $gene_name ) {
		$self->{'__id2gene__'}->{$id}        = $gene_name;
		$self->{'__gene2id__'}->{$gene_name} = $id;
	}
	return $self->{'__id2gene__'}->{$id};
}

=head3 AddDataset ( {Gene_ID => 'Gene name i', Sample_ID => 'Sample name n', value => 'value i n' }, add )

Adds (+) the value to the table or if add==0 replaces the value in the table! 

=cut

sub AddDataset {
	my ( $self, $hash, $add ) = @_;
	unless ( defined $hash->{'Gene_ID'} and defined $hash->{'Sample_ID'} ) {
		return $self;
	}
	$self->Add_unique_key( 'Gene_ID', 'Gene_ID' )
	  unless ( $self->{'uniques'}->{'Gene_ID'} );
	my $key = $self->{'uniques'}->{'Gene_ID'};
	unless ( defined $key->{ $hash->{'Gene_ID'} } ) {
		push( @{ $self->{'data'} }, [ $hash->{'Gene_ID'} ] );
		$key->{ $hash->{'Gene_ID'} } = scalar( @{ $self->{'data'} } - 1 );
	}
	if ($add) {
		@{ @{ $self->{'data'} }[ $key->{ $hash->{'Gene_ID'} } ] }
		  [ $self->Add_2_Header( $hash->{'Sample_ID'} ) ] += $hash->{'value'};
	}
	else {
		@{ @{ $self->{'data'} }[ $key->{ $hash->{'Gene_ID'} } ] }
		  [ $self->Add_2_Header( $hash->{'Sample_ID'} ) ] = $hash->{'value'};
	}
	return $self;
}

sub sample2id {
	my ( $self, $sample_name, $id ) = @_;
	$sample_name =~ s/ spliced//;
	$id ||= $self->{'__sample2id__'}->{$sample_name};
	unless ( defined $id ) {
		$id = scalar( keys %{ $self->{'__sample2id__'} } ) + 1;
	}

	if ( defined $id ) {
		$self->{'__id2sample__'}->{$id}          = $sample_name;
		$self->{'__sample2id__'}->{$sample_name} = $id;
	}
	return $self->{'__sample2id__'}->{$sample_name};
}

sub id2sample {
	my ( $self, $id, $sample_name ) = @_;
	if ( defined $sample_name ) {
		$self->{'__id2sample__'}->{$id}          = $sample_name;
		$self->{'__sample2id__'}->{$sample_name} = $id;
	}
	return $self->{'__id2sample__'}->{$id};
}

=head3 import_database ( filename )

Please make sure you actually store the results of this object as I re-initialize the object - you will get a new one!

=cut

sub import_database {
	my ( $self, $filename ) = @_;
	## reset all internal counters
	$self = ref($self)
	  ->new( { 'data_storage_spliced' => $self->{'data_storage_spliced'} } );
	$self->{'__gene2id__'}          = {};
	$self->{'__id2sample__'}        = {};
	$self->{'__sample2id__'}        = {};
	$self->{'__lastSampleID__'}     = 0;
	$self->{'__lastGeneID__'}       = 0;
	$self->{'__lastDataID__'}       = 0;
	$self->{'uniques'}->{'Gene_ID'} = undef;
	$self->{'data'}                 = [];
	## read in gene an sample information:
	$self->read_file($filename);
	my $obj = stefans_libs::database::Chromium_SingleCell::datavalues->new(
		{
			file                   => $filename,
			'data_storage_spliced' => $self->{'data_storage_spliced'}
		}
	);
	my $sth = $obj->execute_for_search(
		{
			'search_columns' => [
				'datavalues.gene_id', 'datavalues.sample_id',
				'datavalues.value'
			],
			'where' => [],
		}
	);
	my $dat;

	while ( $dat = $sth->fetchrow_hashref() ) {
		$self->AddDataset(
			{
				Gene_ID   => $self->id2gene( $dat->{'gene_id'} ),
				Sample_ID => $self->id2sample( $dat->{'sample_id'} ),
				value     => $dat->{'value'}
			},
			0
		);
	}
	if ( $self->{'data_storage_spliced'} ) {
		my $sth = $obj->execute_for_search(
			{
				'search_columns' => [
					'datavalues.gene_id', 'datavalues.sample_id',
					'datavalues.spliced'
				],
				'where' => [ [ 'datavalues.spliced', '>', 'my_value' ] ],
			},
			[0]
		);
		while ( $dat = $sth->fetchrow_hashref() ) {
			$self->AddDataset(
				{
					Gene_ID   => $self->id2gene( $dat->{'gene_id'} ),
					Sample_ID => $self->id2sample( $dat->{'sample_id'} )
					  . " spliced",
					value => $dat->{'spliced'}
				},
				0
			);
		}
	}
	## now I need to update the internal objects too!
	
	$obj->{'use_this_sql'} = "select max(id) from datavalues";
	my $tmp;
	$tmp = $obj->get_data_table_4_search(
		{
			'search_columns' => ['datavalues.id'],
			'where'          => [],
		}
	);
	$self->{'__lastDataID__'} = @{ @{ $tmp->{'data'} }[0] }[0];

	$obj->{'use_this_sql'} = "select max(id) from samples";
	$tmp = $obj->get_data_table_4_search(
		{
			'search_columns' => ['datavalues.id'],
			'where'          => [],
		}
	);
	$self->{'__lastSampleID__'} = @{ @{ $tmp->{'data'} }[0] }[0];

	$obj->{'use_this_sql'} = "select max(id) from genes";
	$tmp = $obj->get_data_table_4_search(
		{
			'search_columns' => ['datavalues.id'],
			'where'          => [],
		}
	);
	$self->{'__lastGeneID__'} = @{ @{ $tmp->{'data'} }[0] }[0];

	return $self;
}

sub read_file {
	my ( $self, $file ) = @_;
	$self->Add_2_Header('Gene_ID');
	$self->{'filename'} = $file;
	if ( defined $file and -f $file ) {		
		my $obj = stefans_libs::database::Chromium_SingleCell::datavalues->new(
			{
				file                   => $file,
				'data_storage_spliced' => $self->{'data_storage_spliced'}
			}
		);
		my $tmp = $obj->{'data_handler'}->{'genes'}->get_data_table_4_search(
			{
				'search_columns' => [ 'genes.id', 'gname' ],
				'where'          => [],
			}
		);

	  #warn $obj->{'data_handler'}->{'genes'}->{'complex_search'}." failing?\n";
		for ( my $i = 0 ; $i < $tmp->Lines() ; $i++ ) {
			$self->gene2id( reverse @{ @{ $tmp->{'data'} }[$i] } );
			$self->{'__lastGeneID__'} = @{ @{ $tmp->{'data'} }[$i] }[0];
		}

		$tmp = $obj->{'data_handler'}->{'samples'}->get_data_table_4_search(
			{
				'search_columns' => [ 'samples.id', 'sname' ],
				'where'          => [],
			}
		);

	#warn $obj->{'data_handler'}->{'samples'}->{'complex_search'}." failing?\n";
		for ( my $i = 0 ; $i < $tmp->Lines() ; $i++ ) {
			$self->sample2id( reverse @{ @{ $tmp->{'data'} }[$i] } );
			$self->Add_2_Header( @{ @{ $tmp->{'data'} }[$i] }[1] );
			$self->{'__lastSampleID__'} = @{ @{ $tmp->{'data'} }[$i] }[0];
		}
		$obj->{'use_this_sql'} = "select max(id) from datavalues";

		$tmp = $obj->get_data_table_4_search(
			{
				'search_columns' => ['datavalues.id'],
				'where'          => [],
			}
		);

#warn $obj->{'complex_search'}." failing (@{ @{ $tmp->{'data'} }[0] }[0]) ?". $tmp->Lines()."\n";
		$self->{'__lastDataID__'} = @{ @{ $tmp->{'data'} }[0] }[0]
		  if ( $tmp->Lines() > 0 );
	}
	return $self;
}

=head3 import_tables ( @path )

import_tables can import a set of tables creates by the print2table() function in order to collide them into a databse.
The function can also be used to read in cellRangers matrix files.

=cut

sub import_tables {
	my ($self, @path ) = @_;
	## first check if all files are existing:
	my ($ifile, @line);
	foreach my $path ( @path ) { 
		my $err = '';
		foreach ( 'barcodes.tsv', 'genes.tsv' , 'matrix.mtx') {
			unless ( -f File::Spec->catfile($path, $_ ) ){
				$err .= "file $_ missing in path $path\n";
			}
		}
		if ( $err =~ m/\w/ ){
			warn $err;
			return 0;
		}
		$ifile = File::Spec->catfile($path, 'barcodes.tsv'); 
		open ( IN, "<$ifile" ) or die "I could not open the file '$ifile'\n$!";
		while ( <IN> ) {
			chomp();
			$self->sample2id( $_ );
			$self->Add_2_Header( $_ );
		}
		close ( IN );
		
		$ifile = File::Spec->catfile($path, 'genes.tsv'); 
		open ( IN, "<$ifile" ) or die "I could not open the file '$ifile'\n$!";
		while ( <IN> ) {
			chomp();
			@line = split(" ",$_);
			$self->gene2id( $line[0] );
		}
		close ( IN );
		
		$ifile = File::Spec->catfile($path, 'matrix.mtx'); 
		open ( IN, "<$ifile" ) or die "I could not open the file '$ifile'\n$!";
		my $i = 0;
		while ( <IN> ) {
			$i ++;
			if ($_ =~ m/^%/) {
				while ( <IN> ) { ## skip the 3 first lines
					$i ++;
					last if ( $i == 3);
				}
				## set mode to 
				$self->{'data_storage_spliced'} = 0;
			}else {
				chomp();
				@line = split(" ",$_);
				$self->AddDataset(
					{
						Gene_ID   => $self->id2gene( $line[0] ),
						Sample_ID => $self->id2sample( $line[1] ),
						value     => $line[2]
					} , 0
				);
				if ( @line == 4) {
					$self->AddDataset(
						{
							Gene_ID   => $self->id2gene( $line[0] ),
							Sample_ID => $self->id2sample( $line[1] ). " spliced",
							value     => $line[3]
						} , 0
					) if ( $line[3] > 0);
					$self->{'data_storage_spliced'} = 1;
					warn "Data storeage type set to 'splced'\n";
				}
					
			}
			last;
		}
		while ( <IN> ) {
			chomp();
			@line = split(" ",$_);
			$self->AddDataset(
				{
					Gene_ID   => $self->id2gene( $line[0] ),
					Sample_ID => $self->id2sample( $line[1] ),
					value     => $line[2]
				} , 0
			);
			if ( $self->{'data_storage_spliced'} ) {
				$self->AddDataset(
					{
						Gene_ID   => $self->id2gene( $line[0] ),
						Sample_ID => $self->id2sample( $line[1] ). " spliced",
						value     => $line[3]
					} , 0
				) if ( $line[3] > 0);
			}
		}
		close ( IN );
	}
	return $self;
}

=head3 print2table ($filename, $lines, $type)

Prints the data to a barcodes.tsv,  genes.tsv and  matrix.mtx using almost the cellranger format:

barcode.tsv contains 'only' the cell barcodes
genes.tsv contains the geneID and gene symbol
matrix.mtx the gene_id, sample_id, value and spliced (spliced is optional)

In the matrix.mtx the annotation and total count of values is missing (lines 1-3)
These files are meant to be intermediate results, for each chromosome one.

The easiest way to combine them is to use the object again - read all files and export as databse. 

=cut

sub export2file {
	my ($self, $obj, $ofile ) = @_;
	my @order;
	if ( $ofile =~m/matrix.mtx$/ ){
		@order = ( 2,1,3 );
		if ( $self->{'data_storage_spliced'}){
			push( @order,4);
		}
	}else{
		@order= (1);
	}
	if ( $obj->Lines() > 0 ) {
		if ( -f  $ofile) {
			open ( OUT, ">>$ofile") or die "could not add to file '$ofile'\n$!";
		}else {
			open ( OUT, ">$ofile") or die "could not create file '$ofile'\n$!";
		}
		while( my $line = shift(@{$obj->{'data'}}) ) {
			#warn "export to $ofile in order".join(" ",@order)." the values".join(" ", @$line[@order])."\n";
			print OUT join(" ", @$line[@order])."\n";
		}
		close ( OUT );
	}
	return $self;
}

sub print2table {
	my ( $self, $filename, $lines ) = @_;
	
	$lines ||= $self->Lines();
	$filename ||= $self->{'filename'};
	
	my $fm = root->filemap($filename);
	my $outpath = File::Spec->catfile($fm->{'path'}, $fm->{'filename_base'} );
	mkdir( $outpath ) unless ( -d $outpath);
	my $ofile;	
	
	## genes
	my $genes = $self->prepare_gene_table();
	$ofile = File::Spec->catfile($outpath, 'genes.tsv');
	$self->export2file( $genes,  $ofile);
	
	## samples
	my $samples = $self->prepare_sample_table();
	$ofile = File::Spec->catfile($outpath, 'barcodes.tsv');
	$self->export2file( $samples,  $ofile);

	## data
	my $data = $self->prepare_data_table($lines);
	$ofile = File::Spec->catfile($outpath, 'matrix.mtx');
	$self->export2file($data, $ofile);
	
	return $self,
}

sub print2file {
	my ( $self, $filename, $lines, $type ) = @_;
	$type ||= 'database';
	if ( $type eq "table" ) {
		return $self->print2table( $filename, $lines  );
	} 
	$filename ||= $self->{'filename'};
	
	$lines ||= $self->Lines();
	my $obj = stefans_libs::database::Chromium_SingleCell::datavalues->new(
		{
			'data_storage_spliced' => $self->{'data_storage_spliced'},
			file                   => $filename
		}
	);

	my $batch = stefans_libs::database::SQLiteBatch->new();

	my ( $genes, $samples, $data, $gene_id );

	#	warn $self->info();
	## save the new gene info
	$genes = $self->prepare_gene_table();

	if ( $genes->Lines() > 0 ) {
		$batch->batch_import( $obj->{'genes'}, $genes );
	}

	$samples = $self->prepare_sample_table();
	
	if ( $samples->Lines() > 0 ) {
		$batch->batch_import( $obj->{'samples'}, $samples );
	}

	$data = $self->prepare_data_table($lines);

	if ( $data->Lines() > 0 ) {
		$batch->batch_import( $obj, $data );
	}
	return $self,

}

sub prepare_gene_table {
	my $self = shift;
	
	my $genes = data_table->new();
	$genes->Add_2_Header( [ 'id', 'gname' ] );
	$self->{'as_array'} = {};
	map {
		#		warn "Processing gene $_\n";
		if ( $self->gene2id($_) > $self->{'__lastGeneID__'} ) {
			push( @{ $genes->{'data'} }, [ $self->gene2id($_), $_ ] );

			#	warn "Adding gene $_\n";
		}

	  } @{ $self->GetAsArray('Gene_ID') }
	  ;    ## only if the data is not already in the database
	if ( $genes->Lines() > 0 ) {
		$self->{'__lastGeneID__'} = scalar( keys %{ $self->{'__gene2id__'} } );
	}
	return $genes;
}

sub prepare_sample_table {
	my $self = shift;
	
	my $samples = data_table->new();
	$samples->Add_2_Header( [ 'id', 'sname' ] );
	$self->{'as_array'} = {};
	my $ids;
	map {
		if ( $self->sample2id($_) > $self->{'__lastSampleID__'} ) {
			unless ( $ids->{$self->sample2id($_)} ){
				push( @{ $samples->{'data'} }, [ $self->sample2id($_), $_ ] );
				$ids->{$self->sample2id($_)} = 1;
			#	warn "Adding sample $_ with id ".$self->sample2id($_)."\n";
			}	
		}

	  } grep { !/Gene_ID/ } @{ $self->{'header'} };
	  ;    ## only if the data is not already in the database
	if ( $samples->Lines() > 0 ) {
		$self->{'__lastSampleID__'} =
		  scalar( keys %{ $self->{'__sample2id__'} } );
	}
	return $samples;
}

sub prepare_data_table {
	my ($self, $lines)  = @_;
	my $data = data_table->new();

	my ( $gene_id );

	my ( $hash, $data_id, $sname, $key, $prog );
	## get last ID
	$data_id = $self->{'__lastDataID__'};
	$prog = 0;
	local $| = 1; ## local progress bar
	if ( $self->{'data_storage_spliced'} ) {
		$data->Add_2_Header(
			[ 'id', 'sample_id', 'gene_id', 'value', 'spliced' ] );
		$data->Add_unique_key( 'data_sg', [ 'sample_id', 'gene_id' ] );
		while ( my $line = shift( @{ $self->{'data'} } ) ) {
			@$hash{ @{ $self->{'header'} } } = @$line;

			#print "one line = \$exp = ".root->print_perl_var_def($hash ).";\n";
			$gene_id = $self->gene2id( $hash->{'Gene_ID'} );
			if ( $prog++ % 100 == 0) {
				print "-$prog-";
			}
			map {
				if ( defined $hash->{$_} )
				{    ## simple only the sample id no splicing
					unless ( $_ =~ m/spliced/ ) {
						push(
							@{ $data->{'data'} },
							[
								++$data_id, $self->sample2id($_),
								$gene_id,   $hash->{$_},
								0
							]    ## default spliced to 0
						) if ( defined $hash->{$_} );
						$key = [ $self->sample2id($_), $gene_id];
						$key = "@$key";
						$data->{'uniques'}->{'data_sg'}->{$key} = $data->Lines()-1;
					}
					else {
						## find the right ID!
						$sname = $_;
						$sname =~ s/ spliced//;
						$sname = $data->{'uniques'}->{'data_sg'}
						  ->{ join( " ", $self->sample2id($sname), $gene_id ) };

						@{ @{ $data->{'data'} }[$sname] }[4] = $hash->{$_};
					}
				}
			} grep { !/Gene_ID/ } @{ $self->{'header'} };
			last if ( ( --$lines ) == 0 );
		}
	}
	else {
		$data->Add_2_Header( [ 'id', 'sample_id', 'gene_id', 'value' ] );
		$data->Add_unique_key( 'data_sg', [ 'sample_id', 'gene_id' ] );
		while ( my $line = shift( @{ $self->{'data'} } ) ) {
			@$hash{ @{ $self->{'header'} } } = @$line;

			#print "one line = \$exp = ".root->print_perl_var_def($hash ).";\n";
			$gene_id = $self->gene2id( $hash->{'Gene_ID'} );

			map {
				if ( defined $hash->{$_} )
				{    ## simple only the sample id no splicing

					push(
						@{ $data->{'data'} },
						[
							++$data_id, $self->sample2id($_),
							$gene_id,   $hash->{$_},
							0
						]    ## default spliced to 0
					) if ( defined $hash->{$_} );
					$key = [ $self->sample2id($_), $gene_id];
					$key = "@$key";
					$data->{'uniques'}->{'data_sg'}->{$key} = $data->Lines()-1;
				}
			} grep { !/Gene_ID/ } @{ $self->{'header'} };
			last if ( ( --$lines ) == 0 );
		}
	}
	if ( $data->Lines() > 0 ) {
		$self->{'__lastDataID__'} = $data_id;
	}
	$| = 0; ## local progress bar off.
	
	$self->{'uniques'}->{'Gene_ID'} =
	  undef; 
	return $data;
}


1;
