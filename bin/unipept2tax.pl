#!/usr/bin/perl
use strict;
use warnings;
no warnings qw(uninitialized);
#use diagnostics;

use Bio::DB::Taxonomy;
#require "/work/rec3141/OSD-analysis/flatfile.pm";
use Data::Dumper;

#BEGIN { *Bio::DB::Taxonomy::flatfile::DESTROY = sub {} }
#sub Bio::DB::Taxonomy::flatfile::DESTROY {};


die
"\nthis program takes as inputs:
1) a taxon name  -- e.g. 'Deltaproteobacteria' will find all Deltaproteobacteria
  OR a filename containing a list of Taxon Ids or Taxon Names or a mixture of both
  (taxa to exclude may be prefixed by a '-')

e.g. unipept2tax.pl Deltaproteobacteria

and outputs the full taxonomy in the following format:
taxid\tk__XYZ|p__XYZ|c__XYZ|o__XYZ|f__XYZ|g__XYZ|s__XYZ
\n\n" if ( scalar(@ARGV) != 1 );


### SET UP TAXONOMIC DATABASE

my $idx_dir = '/Volumes/ramdisk/taxonomy/';

my ( $nodefile, $namesfile ) = ( 'nodes.dmp', 'names.dmp' );

my $db = new Bio::DB::Taxonomy( -source    => 'flatfile',
                                -nodesfile => $idx_dir . $nodefile,
                                -namesfile => $idx_dir . $namesfile,
                                -directory => $idx_dir
);

die "Couldn't load taxonomic database\n\n" unless $db;



### FIND WHICH TAXONOMIC GROUPS WERE REQUESTED

my $printname = $ARGV[0]; $printname =~ s#[ /']#-#g;
my %taxahash;

my @ancestors; #taxon ids of ANCESTOR
my @requestednames;

# IS IT A FILE?
if ( open( IF, "<$ARGV[0]" ) ) {
    while (<IF>) {
		chomp;
    	my @line = split(/\t/,$_);
    	# list can be taxonids or taxa names or mixture
		   $taxahash{$line[0]}{'myparent_id'} = $line[2];
		   $taxahash{$line[0]}{$line[1]} = $line[0];
		   push(@ancestors,$line[2]);
		   push(@requestednames,$line[0]);
	}
    close(IF);
}

die "couldn't find requested taxonid\n" unless (@requestednames);
#die "unequal data input\n" unless length(keys %taxahash)==scalar(@requesteids);

#print Dumper(%taxahash);

#### DO WORK

### Get all ancestors of the requested taxa

foreach my $name (keys %taxahash) {
#	print "$name\t" . $taxahash{$name}{'myparent_id'} . "\n";
	my $node = $db->get_taxon(-taxonid => $taxahash{$name}{'myparent_id'});
	next unless ($node);
	while ( ! exists($taxahash{$name}{'superkingdom'}) ) {
#		 print $node->scientific_name . "(" . $node->rank . ") --> ";
		if ( ! exists($taxahash{$name}{$node->rank}) ) { $taxahash{$name}{$node->rank} = $node->scientific_name }
		$node = $db->ancestor($node);
		last unless defined($node);
	}
 #print "\n";

}

### OUTPUT NECESSARY FILES
my @okranks = qw/superkingdom phylum class order family genus species/;
my @rankabb = qw/k__ p__ c__ o__ f__ g__ s__/;

open( OUTTAXA,  ">" . $printname . ".tax" );

my @printout;
#print Dumper(%taxahash);

for(my $j=0;$j<scalar(@requestednames);$j++) {

	my $k1 = $requestednames[$j];
	
	my @outline;
	for(my $i=0;$i<scalar(@okranks);$i++) {
#		print "$i $okranks[$i] ";
		if ( exists( $taxahash{ $k1 }{ $okranks[$i] } ) ) {
			push(@outline, $rankabb[$i] . $taxahash{ $k1 }{ $okranks[$i] });
		} else {
			push(@outline, $rankabb[$i]);
		}
	}

	push(@printout, $k1 . "\t" . join("|",@outline));
}

print OUTTAXA join("\n",@printout),"\n";

close(OUTTAXA);

print "DONE $printname\n";
exit 0;
