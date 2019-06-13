#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 2;
BEGIN { use_ok 'stefans_libs::database::Chromium_SingleCell::samples' }

my ( $value, @values, $exp );
my $obj = stefans_libs::database::Chromium_SingleCell::samples -> new();
is_deeply ( ref($obj) , 'stefans_libs::database::Chromium_SingleCell::samples', 'simple test of function stefans_libs::database::Chromium_SingleCell::samples -> new()' );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";


