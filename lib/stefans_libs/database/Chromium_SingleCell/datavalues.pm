package stefans_libs::database::Chromium_SingleCell::datavalues;


#  Copyright (C) 2010 Stefan Lang

#  This program is free software; you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation;
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

use stefans_libs::database::variable_table;
use base variable_table;

use stefans_libs::database::Chromium_SingleCell::genes;
use stefans_libs::database::Chromium_SingleCell::samples;
use stefans_libs::database::SQLiteBatch;


##use some_other_table_class;
use POSIX; ## we need floor

use strict;
use warnings;


sub new {

    my ( $class, $file, $debug ) = @_;
    
#    $file ||= "Chromium_SingleCell_Perl.db";
#    my $driver   = "SQLite";
#	#my $database = "test.db";
#	my $dsn = "DBI:$driver:dbname=$file";
#	my $userid = "";
#	my $password = "";
#	my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 })
#                      or die $DBI::errstr;    

    my ($self);

    $self = {
        debug => $debug,
#        dbh   => $dbh,
		'data_storage_spliced' => 0,
        connection => {
        	'driver' => "SQLite",
        	'filename' => $file,
        },
    };

    bless $self, $class if ( $class eq "stefans_libs::database::Chromium_SingleCell::datavalues" );
    
    $self->{'dbh'} = $self->getDBH();
    
    Carp::confess ("$class : new -> we need a acitve database handle at startup!, not "
	  . ref($self->{'dbh'}))
	  unless ( ref($self->{'dbh'}) =~ m/::db$/ );
	  
    $self->init_tableStructure();
    
    return $self;

}

sub insert {
	my ( $self, $hash ) = @_;
	my $err = '';
	map {$err .= "no $_ info!\n"  unless ( defined $hash->{$_} )} 'sample', 'gene', 'value';
	Carp::confess ( $err ) if ( $err =~ m/\w/ );
	my $gene_id = $self->{data_handler}->{'genes'} -> AddDataset( { 'gname' => $hash->{'gene'}});
	my $sample_id = $self->{data_handler}->{'samples'} -> AddDataset( { 'sname' => $hash->{'sample'}});
	return $self->AddDataset( { value => $hash->{'value'}, gene_id => $gene_id, sample_id => $sample_id});
}

sub DESTROY {
      my $self = shift;
      $self->dbh()->disconnect();
	  print "db saved as $self->{'connection'}->{'filename'}\n";
}

sub  init_tableStructure {
     my ($self, $dataset) = @_;
     my $hash;
     $hash->{'INDICES'}   = [];
     $hash->{'UNIQUES'}   = [];
     $hash->{'variables'} = [];
     $hash->{'table_name'} = "datavalues";

	 push(
		@{ $hash->{'variables'} },
		{
			'name'         => 'sample_id',
			'type'         => 'INTEGER UNSIGNED',
			'NULL'         => '0',
			'description'  => '',
			'data_handler' => 'samples',
			'needed'       => ''
		}
	 );
	 push(
		@{ $hash->{'variables'} },
		{
			'name'         => 'gene_id',
			'type'         => 'INTEGER',
			'NULL'         => '0',
			'description'  => '',
			'data_handler' => 'genes',
			'needed'       => ''
		}
	 );
	 
	 push(
		@{ $hash->{'variables'} },
		{
			'name'         => 'value',
			'type'         => 'INTEGER',
			'NULL'         => '0',
			'description'  => '',
			'needed'       => ''
		}
	 );
	 push(
		@{ $hash->{'variables'} },
		{
			'name'         => 'spliced',
			'type'         => 'INTEGER',
			'NULL'         => '0',
			'description'  => '',
			'needed'       => ''
		}
	 );
	 
     $self->{'table_definition'} = $hash;
     $self->{'UNIQUE_KEY'} = [ 'sample_id','gene_id' ];
	
     $self->{'table_definition'} = $hash;

     $self->{'_tableName'} = $hash->{'table_name'}  if ( defined  $hash->{'table_name'} ); # that is helpful, if you want to use this class without any variable tables

     ##now we need to check if the table already exists. remove that for the variable tables!
     unless ( $self->tableExists( $self->TableName() ) ) {
     	$self->create();
     }
     ## Table classes, that are linked to this class have to be added as 'data_handler',
     ## both in the variable definition and here to the 'data_handler' hash.
     ## take care, that you use the same key for both entries, that the right data_handler can be identified.
     $self->{'samples'} = $self->{'data_handler'}->{'samples'} = 
     	stefans_libs::database::Chromium_SingleCell::samples->new($self->{'dbh'}, $self->{'connection'}->{'filename'}, $self->{'debug'});
    $self->{'genes'} = $self->{'data_handler'}->{'genes'} = 
     	stefans_libs::database::Chromium_SingleCell::genes->new($self->{'dbh'}, $self->{'connection'}->{'filename'}, $self->{'debug'});
     	
     return $dataset;
}

=head3 store_data_table ( $self, $data_table, $gene_colname)

Uses the samples == columns and genes == rows table and stores it as SQlite database.

=cut
 
sub store_data_table {
	my ( $self, $data_table, $gene_colname ) = @_;
	my ( $samples, $genes, $data, $v ) = @_;
	my $batch = stefans_libs::database::SQLiteBatch->new();
	
	$samples = data_table->new();
	my $i = 1;
	$samples ->Add_db_result( ['id', 'sname'], [map { [ $i++, $_] } @{$data_table->{'header'}} ] );
	$batch -> batch_import ( $self->{'samples'}, $samples );
	undef($samples);
	
	$genes = data_table->new();
	$i = 1;
	$genes ->Add_db_result( ['id', 'gname'], [map { [ $i++, $_] } @{$data_table->GetAsArray($gene_colname)} ] );
	$batch -> batch_import ( $self->{'genes'}, $genes );
	undef($genes);
	($samples) = $data_table->Header_Position($gene_colname);
	
	$data = data_table->new();
	$data-> Add_2_Header(['id', 'sample_id', 'gene_id', 'value', 'spliced']);
	$i = 1;
	my $sample_ID = 0;
	if ( $samples == 0 ){
		$sample_ID =1;
	}
	if ( $self->{'data_storage_spliced'}) {
		for (my $gene_id = 0; $gene_id < $data_table->Rows(); $gene_id ++ ){
				for (my $sample_id = $sample_ID; $sample_id < $data_table->Columns();$sample_id +=2 ){
					next if ( $sample_id == $samples);
					$v = @{@{$data_table->{'data'}}[$gene_id]}[$sample_id];
					push(@{$data->{'data'}}, [$i ++, floor($sample_id/2)+1, $gene_id+1,$v, @{@{$data_table->{'data'}}[$gene_id]}[$sample_id+1]  ]) 
						if ( defined $v and $v != 0 );
			}
		}
	}
	else {
			for (my $gene_id = 0; $gene_id < $data_table->Rows(); $gene_id ++ ){
				for (my $sample_id = $sample_ID; $sample_id< $data_table->Columns();$sample_id ++ ){
					next if ( $sample_id == $samples);
					$v = @{@{$data_table->{'data'}}[$gene_id]}[$sample_id];
					push(@{$data->{'data'}}, [$i ++, $sample_id+1, $gene_id+1, $v, 0 ]) 
						if ( defined $v and $v != 0 );
				}
			}
	}
	$batch -> batch_import ($self, $data );
	undef($data);
	return $self;
}


sub expected_dbh_type {
	return 'dbh';
	#return 'database_name';
}


1;
