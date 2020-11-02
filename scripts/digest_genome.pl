#!/bin/perl

@arg=split/\//, $0; $scriptname=$arg[$#arg];
use lib '/home/gthansen/bin/RestrictionDigest/lib';
use RestrictionDigest;
use Getopt::Long;

######### Author: Grace Hansen ###########
# This script uses RestrictionDigest to digest a given genome according to a given restriction enzyme.
# It only works with a single restriction enzyme (although RestrictionDigest can process double enzymes)
# The genome must be a fa with one line per chromosome

# OPTIONS

# get option
%options=(
		'f|fasta=s'    => \$fasta,
		'r|enzyme=s'    => \$enzyme,
		'o|outdir=s'	=> \$outdir,
         );

GetOptions(%options) || die;


my $single_digest=RestrictionDigest::SingleItem::Single->new();
	$single_digest->add_ref(-reference=>"$fasta");
	$single_digest->add_single_enzyme(-enzyme =>"$enzyme");
	$single_digest->change_range(-start =>150,-end =>800);
	$single_digest->change_lengths_distribution_parameters(-front =>200,-behind =>800,-step =>25);
	$single_digest->add_output_dir(-output_dir=>"$outdir");
	$single_digest->single_digest();