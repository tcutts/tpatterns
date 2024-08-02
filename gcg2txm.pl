#!/bin/perl

#######################################################
#                                                     #
#  This code is copyright (C) T J R Cutts 1998        #
#  tjrc1@mole.bio.cam.ac.uk                           #
#                                                     #
#  It may be redistributed and modified, as long as   #
#  this copyright notice is retained.                 #
#                                                     #
#######################################################

# This script converts a GCG codon translation table into
# the format expected by tpatterns.

# Usage: gcg2txm.pl gcgtable.txt > tpattable.txm

undef $started;
undef %table;

foreach (<>)
{
    if (defined $started)
    {
		tr/acgt/ACGT/;
		if (/^\s*([A-Z*])\s+\S\S\S\s+(([ACGT]{3}\s+)+)\s+!/)
		{
	    	$aa = $1;
	    	@f = split(/\s+/, $2);
	    	foreach $codon (@f)
	    	{
				$table{$codon} = $aa;
	    	}
		}
    }
    else
    {
		if (/\.\.\s*$/)
		{
	    	$started = 1;
		}
    }
}

foreach $a ("A","C","G","T")
{
    foreach $b ("A","C","G","T")
    {
		foreach $c ("A","C","G","T")
		{	
	    	$codon = "$a$b$c";
	    	print "$codon\t";
	    	print defined $table{$codon} ? $table{$codon} : "x";
	    	print "\n";
		}
    }
}

