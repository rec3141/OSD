#!/usr/bin/perl -w
use strict;

use Bio::DB::Taxonomy;
#use Data::Dumper;

 use warnings;
 use diagnostics;

#print("\nthis program updates the taxonomy database\n\n");

my $idx_dir = '/Volumes/ramdisk/taxonomy/';

my ( $nodefile, $namesfile ) = ( "nodes.dmp", "names.dmp" );

my $db = new Bio::DB::Taxonomy( -source    => 'flatfile',
                                -nodesfile => $idx_dir . $nodefile,
                                -namesfile => $idx_dir . $namesfile,
                                -directory => $idx_dir
);
#print("done");
