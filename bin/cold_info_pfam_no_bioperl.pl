#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use File::Slurp;
use File::Basename;
use Data::Dumper;

# to run in parallel
# ls step1-pfam.tax/out.*.tax | parallel --gnu "./cold_info_pfam_no_bioperl.pl {}"


# read in amino acid information
my %aminos;
my @amino_lines = read_file( 'amino-acid-info.txt' ) ;
chomp(@amino_lines);
my @amino_cols = split("\t",shift(@amino_lines));
foreach my $line (@amino_lines) {
	chomp($line);
	my @cols = split("\t",$line);
	for (my $i=1; $i < scalar(@amino_cols); $i++) {
		$aminos{ $cols[0] }{ $amino_cols[$i] } = $cols[$i];
	}
}

#foreach my $aa (sort keys %aminos) {print $aa . "\t"}; print "\n";

my %aa_props; #global protein properties

# Amino acid analysis of Methylobacterium and Rhodococcus sp. has revealed 
# significant genome-wide protein adaptations leading to a general reduction 
# in the number of acidic residues, the arginine to lysine ratio, aliphaticity 
# and aromaticity. As expected, these changes favor an increase in flexibility 
# at lower temperatures.

my $file = $ARGV[0];
my ($name,$path,$suffix) = fileparse($file, ".pfam.tax");
$name =~ s/^out\.//;

	%aa_props = ();

	print "$file ... ";
	#(ugly) output from post-kraken.sh
	#HWI-M02024:121:000000000-ADJ6A:1:2111:24359:11410-1:N:0:ACGGAACA	_1_415_-	13	78	137	9.1E-5	PF00004	ATPase family associated with various cellular activities (AAA)	IPR003959	ATPase, AAA-type, core	GO:0005524		AKFRGQFEERLRSVLEEVSRPDSGVVLFVDELHTVVGSERSSTDAGSLLKPALARGDLRCIGATTP	LISRDLGALIAGAKFRGQFEERLRSVLEEVSRPDSGVVLFVDELHTVVGSERSSTDAGSLLKPALARGDLRCIGATTPEDYRRTVEKDPALNRRFQQVLIREPDLPLSLVILRGLKERYELHHGVSITDEALQAANR	d__Bacteria|p__Cyanobacteria|o__Chroococcales|g__Synechococcus|s__Synechococcus_sp._CC9902
	my @krakenlines = read_file($file);
	chomp(@krakenlines);

	my @faaline;
	foreach my $line (@krakenlines) {
		next if $line =~ m/^\s*$/;
		my @splitline = split("\t",$line);
		push(@faaline, join("",$splitline[0],$splitline[1],"\t",$splitline[13]));
	}

# get protein sequences
	print "got sequences ... ";

	# calculate parameters of interest
	my %aa_count = count_aa(\@faaline);
	
	arg_lys_ratio(\%aa_count); # R/(R+K) -- cold = decrease, slope = +
	acidic_fraction(\%aa_count); # (D+K)/all -- cold = decrease, slope = +
	aliphaticity(\%aa_count); # cold = decrease, slope = +
	aliphatic_index(\%aa_count); # cold = decrease, slope = +
	gravy(\%aa_count); # cold = , slope = 
#	proline_fraction(\%aa_count); # P/all -- cold = decrease, slope = +
	weigh_prot(\%aa_count); # average molecular weight of proteins
	count_n(\%aa_count); #average nitrogens per aa
	each_aa_fraction(\%aa_count);
	basic_fraction(\%aa_count);
	aromatic_fraction(\%aa_count);
	hbond_fraction(\%aa_count);
	sulfur_fraction(\%aa_count);
	polar_fraction(\%aa_count);
	nonpolar_fraction(\%aa_count);

	print "calculated params ... ";

	#need to reimport file to reset counters

	open(OUT, ">", "out.$name.protparam.csv");
#	method	read_id	length	acidic	aliphaticity	aliphatic_index	arg_lys	gravy	polar	nonpolar	sulfur	hbond	basic	aromatic	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	nitrogen	mweight
	
	foreach my $key (sort { $aa_count{$a}{'id'} <=> $aa_count{$b}{'id'} } keys %aa_props) {
		print OUT "rec\t" . $key . "\t" .
			$aa_count{$key}{'length'} . "\t" .
			sprintf("%.5f", $aa_props{$key}{'acidic_fraction'}) . "\t" . #%.5f for more precision
			sprintf("%.5f", $aa_props{$key}{'aliphaticity'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'aliphatic_index'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'arg_lys_ratio'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'gravy'}) . "\t" . 
			sprintf("%.5f", $aa_props{$key}{'polar'}) . "\t" . 
			sprintf("%.5f", $aa_props{$key}{'nonpolar'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'sulfur'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'hbond'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'basic'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'aromatic'}) . "\t";
		
		foreach my $aa (sort keys %aminos) {
			print OUT
			sprintf("%.5f", $aa_props{$key}{$aa . '_fraction'}) . "\t";
		}
		
		print OUT 
			sprintf("%.5f", $aa_props{$key}{'nitrogen'}) . "\t" .
			sprintf("%.5f", $aa_props{$key}{'mweight'}) . "\n";
		}

	print "output in out.$name.protparam.csv \n";

	close(OUT);





#------------------------
# SUB count_aa
# hash all amino acids
#------------------------

sub count_aa {
	my $_faaline_ref = shift;
	my @_faaline = @{$_faaline_ref};
	my %_aa_count;
	my $id = 1;
	my $inseq;

	foreach my $_seq_in (@_faaline) {
		if (length($_seq_in) > 0) {
			my @seqin = split("\t",$_seq_in);
			$_aa_count{$seqin[0]}{'id'} = $id++;
			$_aa_count{$seqin[0]}{'length'} = length($seqin[1]);
			my @characters = split(//, $seqin[1]);
			for my $char (keys %aminos) {
				$_aa_count{$seqin[0]}{$char} = 0;
			}
		for my $char (@characters) {
			$_aa_count{$seqin[0]}{$char}++;
		}
	 } else {print "error in sequence " . $_seq_in . "\n";}
	}
	
#	print Dumper(%_aa_count);
	return(%_aa_count);
}




#------------------------
# SUB acidic_fraction
# acidic residue fraction (aspartate+glutamate)/all
#------------------------

sub acidic_fraction {
	my $_aa_count_ref = shift;
	my $_asp = 'D';
	my $_glu = 'Q';
	foreach my $key (keys %{$_aa_count_ref}) {	
		$aa_props{$key}{'acidic_fraction'} = 'NaN';
		$aa_props{$key}{'acidic_fraction'} = ( $_aa_count_ref->{$key}->{$_asp} + $_aa_count_ref->{$key}->{$_glu} ) / $_aa_count_ref->{$key}->{'length'};
	}

	return(0);
}



#------------------------
# SUB basic_fraction
# basic residue fraction
#------------------------

sub basic_fraction {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'basic'} = 'NaN';
		my $tmp = 0;
		foreach my $_aa (keys %aminos) {
			$tmp += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'basic'} ;
		}
		$aa_props{$key}{'basic'} = $tmp / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}


#------------------------
# SUB aromatic_fraction
# aromatic residue fraction
#------------------------

sub aromatic_fraction {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'aromatic'} = 'NaN';
		my $tmp = 0;
		foreach my $_aa (keys %aminos) {
			$tmp += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'aromatic'} ;
		}
		$aa_props{$key}{'aromatic'} = $tmp / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}


#------------------------
# SUB hbond_fraction
# h-bonding residue fraction
#------------------------

sub hbond_fraction {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'hbond'} = 'NaN';
		my $tmp = 0;
		foreach my $_aa (keys %aminos) {
			$tmp += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'hbond'} ;
		}
		$aa_props{$key}{'hbond'} = $tmp / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}



#------------------------
# SUB sulfur_fraction
# sulfur residue fraction
#------------------------

sub sulfur_fraction {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'sulfur'} = 'NaN';
		my $tmp = 0;
		foreach my $_aa (keys %aminos) {
			$tmp += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'sulfur'} ;
		}
		$aa_props{$key}{'sulfur'} = $tmp / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}


#------------------------
# SUB polar_fraction
# polar residue fraction
#------------------------

sub polar_fraction {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'polar'} = 'NaN';
		my $tmp = 0;
		foreach my $_aa (keys %aminos) {
			$tmp += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'polar'} ;
		}
		$aa_props{$key}{'polar'} = $tmp / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}




#------------------------
# SUB nonpolar_fraction
# polar residue fraction
#------------------------

sub nonpolar_fraction {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'nonpolar'} = 'NaN';
		my $tmp = 0;
		foreach my $_aa (keys %aminos) {
			$tmp += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'nonpolar'} ;
		}
		$aa_props{$key}{'nonpolar'} = $tmp / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}



#------------------------
# SUB arg_lys_ration
# arginine/(arginine + lysine)
#------------------------

sub arg_lys_ratio {
	my $_aa_count_ref = shift;
	my $_arg = 'R';
	my $_lys = 'K';
	foreach my $key (keys %{$_aa_count_ref}) {
		$aa_props{$key}{'arg_lys_ratio'} = 'NaN';
		if ( ( $_aa_count_ref->{$key}->{$_lys} + $_aa_count_ref->{$key}->{$_arg} ) > 0) {
			$aa_props{$key}{'arg_lys_ratio'} = $_aa_count_ref->{$key}->{$_arg} / ( $_aa_count_ref->{$key}->{$_lys} + $_aa_count_ref->{$key}->{$_arg} );		
		}
	}

	return(0);
}


#------------------------
# SUB aliphaticity
# aliphatic residue fraction
#------------------------

sub aliphaticity {
	my $_aa_count_ref = shift;

	my @_aliphatics;
	foreach my $aa (keys %aminos) {
		push(@_aliphatics,$aa) if $aminos{$aa}{'aliphatic'} > 0;
	}

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'aliphaticity'} = 'NaN';
		my $tmpali = 0;
		foreach my $_ali (@_aliphatics) {
			$tmpali += $_aa_count_ref->{$key}->{$_ali} ;
		}
		$aa_props{$key}{'aliphaticity'} = $tmpali / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}


#------------------------
# SUB aliphatic_index
# volumetric aliphatic index
#------------------------

sub aliphatic_index {
	my $_aa_count_ref = shift;

	my %_aliphatics = (
        'A' => 1,
        'V' => 2.9,
        'I' => 3.9,
        'L' => 3.9
    );

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'aliphatic_index'} = 'NaN';
		my $tmpali = 0;
		foreach my $_ali (keys %_aliphatics) {
			$tmpali += $_aa_count_ref->{$key}->{$_ali} * $_aliphatics{$_ali} ;
		}
		$aa_props{$key}{'aliphatic_index'} = $tmpali / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}




#------------------------
# SUB gravy
# grand average of hydropathicity
#------------------------

sub gravy {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'gravy'} = 'NaN';
		my $tmpgravy = 0;
		foreach my $_aa (keys %aminos) {
			$tmpgravy += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'hydropathicity'} ;
		}
		$aa_props{$key}{'gravy'} = $tmpgravy / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}



#------------------------
# SUB proline_fraction
# proline fraction (proline/all)
#------------------------

sub proline_fraction {
	my $_aa_count_ref = shift;
	my $_pro = 'P';
	foreach my $key (keys %{$_aa_count_ref}) {	
		$aa_props{$key}{'proline_fraction'} = 'NaN';
		$aa_props{$key}{'proline_fraction'} = ( $_aa_count_ref->{$key}->{$_pro} ) / $_aa_count_ref->{$key}->{'length'};
	}

	return(0);
}


#------------------------
# SUB each_aa_fraction
# each amino acid fraction (aa/all)
#------------------------

sub each_aa_fraction {
	my $_aa_count_ref = shift;
	foreach my $_aa (keys %aminos) {
		foreach my $key (keys %{$_aa_count_ref}) {	
			$aa_props{$key}{$_aa . '_fraction'} = 'NaN';
			$aa_props{$key}{$_aa . '_fraction'} = ( $_aa_count_ref->{$key}->{$_aa} ) / $_aa_count_ref->{$key}->{'length'};
		}
	}
	return(0);
}


#------------------------
# SUB count_aas
# hash all amino acids using Protparam (Expasy cgi)
#------------------------

sub count_aas {
	my $_seq_in = shift;
	my %_aas_count;
	while (my $inseq = $_seq_in->next_seq) {
		my $pp = Bio::Tools::Protparam->new(seq=>$inseq->seq);
		$_aas_count{$inseq->primary_id}{'length'} = $inseq->length;
		$_aas_count{$inseq->primary_id}{'aliphatic_index'} = $pp->aliphatic_index();
		$_aas_count{$inseq->primary_id}{'gravy'} = $pp->gravy();
		$_aas_count{$inseq->primary_id}{'half_life'} = $pp->half_life();
		$_aas_count{$inseq->primary_id}{'acidic_fraction'} = $pp->num_neg()/$inseq->length;
		$_aas_count{$inseq->primary_id}{'basic'} = $pp->num_pos()/$inseq->length;
		$_aas_count{$inseq->primary_id}{'theoretical_pI'} = $pp->theoretical_pI();
		if ($pp->AA_comp('K') == 0) {
			$_aas_count{$inseq->primary_id}{'arg_lys_ratio'} = 'NaN';
		} else {
			$_aas_count{$inseq->primary_id}{'arg_lys_ratio'} = $pp->AA_comp('R')/$pp->AA_comp('K');
		}
		$_aas_count{$inseq->primary_id}{'proline_fraction'} = $pp->AA_comp('P');
		

	}
	return(%_aas_count);
}

#------------------------
# SUB weigh_prot
# average molecular weight
#------------------------

sub weigh_prot {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'mw'} = 'NaN';
		my $tmpweight = 0;
		foreach my $_aa (keys %aminos) {
			$tmpweight += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'mw'} ;
		}
		$aa_props{$key}{'mweight'} = $tmpweight / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}


#------------------------
# SUB count_n
# average nitrogens per aa
#------------------------

sub count_n {

	my $_aa_count_ref = shift;

	foreach my $key ( keys %{$_aa_count_ref} ) {
		$aa_props{$key}{'nitrogen'} = 'NaN';
		my $tmpn = 0;
		foreach my $_aa (keys %aminos) {
			$tmpn += $_aa_count_ref->{$key}->{$_aa} * $aminos{$_aa}{'nitrogen'} ;
		}
		$aa_props{$key}{'nitrogen'} = $tmpn / $_aa_count_ref->{$key}->{'length'}
	}

	return(0);
}




#amino-acid-info.txt
# abb1	abb3	amino acid	mw	acidic	basic	polar	non-polar	sulfur	aromatic	aliphatic	H-bonding	hydropathicity	nitrogen
# G	Gly	Glycine	57.021464	0	0	0	0.5	0	0	1	0	-0.4	1
# A	Ala	Alanine	71.037114	0	0	0	1	0	0	1	0	1.8	1
# S	Ser	Serine	87.032029	0	0	1	0	0	0	0	1	-0.8	1
# P	Pro	Proline	97.052764	0	0	0	1	0	0	1	0	-1.6	1
# V	Val	Valine	99.068414	0	0	0	1	0	0	1	0	4.2	1
# T	Thr	Threonine	101.04768	0	0	1	0	0	0	0	1	-0.7	1
# C	Cys	Cysteine	103.00919	0	0	0.5	1	1	0	0	1	2.5	1
# L	Leu	Leucine	113.08406	0	0	0	1	0	0	1	0	3.8	1
# I	Ile	Isoleucine	113.08406	0	0	0	1	0	0	1	0	4.5	1
# N	Asn	Asparagine	114.04293	0	0	1	0	0	0	0	1	-3.5	2
# D	Asp	Aspartate	115.02694	1	0	1	0	0	0	0	1	-3.5	1
# Q	Gln	Glutamine	128.05858	0	0	1	0	0	0	0	1	-3.5	2
# K	Lys	Lysine	128.09496	0	1	1	0	0	0	0	1	-3.9	2
# E	Glu	Glutamate	129.04259	1	0	1	0	0	0	0	1	-3.5	1
# M	Met	Methionine	131.04048	0	0	0	1	1	0	0	0	1.9	1
# H	His	Histidine	137.05891	0	0.5	1	0	0	1	0	1	-3.2	3
# F	Phe	Phenylalanine	147.06841	0	0	0	1	0	1	0	0	2.8	1
# R	Arg	Arginine	156.10111	0	1	1	0	0	0	0	1	-4.5	4
# Y	Tyr	Tyrosine	163.06333	0	0	0.5	1	0	1	0	1	-1.3	1
# W	Trp	Tryptophan	186.07931	0	0	0	1	0	1	0	1	-0.9	2