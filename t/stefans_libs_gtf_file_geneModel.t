#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 2;
BEGIN { use_ok 'stefans_libs::gtf_file::geneModel' }

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );
my $OBJ = stefans_libs::gtf_file::geneModel -> new({'debug' => 1});
is_deeply ( ref($OBJ) , 'stefans_libs::gtf_file::geneModel', 'simple test of function stefans_libs::gtf_file::geneModel -> new() ');

#print "\$exp = ".root->print_perl_var_def($value ).";\n";


