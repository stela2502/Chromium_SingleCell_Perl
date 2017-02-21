#!/usr/bin/env perl
# IMPORTANT: if you delete this file your app will not work as
# expected.  You have been warned.

use inc::Module::Install;

name "Chromium_SingleCell_Perl";
version '0.01';
author 'Stefan Lang <Stefan.Lang@med.lu.se>';

requires 'Text::Levenshtein' => 0;

sub install_scripts_in_path {
        my ( $dir ) = @_;
        my $PATH;
        opendir( $PATH, "$dir/" ) or die "I could not open the bin folder\n$!\n";
        my @files =  readdir($PATH);
        closedir ( $PATH );
        warn "processing bin folder $dir\n";
        map { if ( $_ =~ m/^\./ ){
                } elsif ( -f "$dir/$_" and $_ =~ m/.pl$/ ){
                        &install_script("$dir/$_");
                } elsif ( -d "$dir/$_" ) {
                        &install_scripts_in_path("$dir/$_")
                } } @files;

}

&install_scripts_in_path('bin');

auto_install();
WriteAll();
