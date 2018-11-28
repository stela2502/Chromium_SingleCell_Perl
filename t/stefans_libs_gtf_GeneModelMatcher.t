#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 2;
BEGIN { use_ok 'stefans_libs::gtf_file::GeneModelMatcher' }

use File::Spec::Functions;

use stefans_libs::root;
use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );

my $filename = File::Spec->catfile($plugin_path, 'data','GOI.annotation.genecode.v12.gff3');

my $OBJ = stefans_libs::gtf_file::GeneModelMatcher -> new({'debug' => 1, 'filename' => $filename});
is_deeply ( ref($OBJ) , 'stefans_libs::gtf_file::GeneModelMatcher', 'simple test of function stefans_libs::gtf_file::GeneModelMatcher -> new() ');

#print "\$exp = ".root->print_perl_var_def( $OBJ->{'data'}  ).";\n";
is_deeply( 5, scalar( @{$OBJ->{'data'}}), "5 gene models? (".scalar( @{$OBJ->{'data'}}).")" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";


