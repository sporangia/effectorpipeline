#!/usr/bin/perl

use strict;
use warnings;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_v $opt_u $opt_o);

# Usage
my $usage = '
trimfasta - a program to trim translated amino acid sequence in fasta
format, with the output to be on a single line with M signal peptide 
at the start position.

		      by
		Sydney E. Everhart
		  April 2014

Copyright (c) 2014 Sydney E. Everhart.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (canopyecology@gmail.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl bamcounter.pl options
 required:
  -a	fasta file.
 optional:
  -u    minimum length of polypeptide [default is 10].
  -o    output file name [default is outfile.fasta].
  -v	verbose mode [optional T/F, default is F].

';

# command line processing.
getopts('a:v:u:o:');
die $usage unless ($opt_a);

my ($inf, $verb, $minl, $outf);

$inf	= $opt_a if $opt_a;
$verb	= $opt_v ? $opt_v : "F";
$minl   = $opt_u ? $opt_u : 10;
$outf   = $opt_o ? $opt_o : "outfile.fasta";

##### ##### ##### ##### #####
# Main.
open(my $in,  "<",  $inf)  or die "Can't open $inf: $!";
open(my $out, ">", $outf)  or die "Can't open $outf: $!";

chomp(my $header = <$in>);
chomp(my $seq = <$in>);

$header =~ s/\s+$//;  #remove trailing whitespace
$seq =~ s/\s+$//;

while (<$in>) { # assigns each line in turn to $_
    if(/^>/) {
        $seq =~ s/^[^M]*//;
        if(length $seq >= $minl){
            print $out "$header\n";
            print $out "$seq\n";
        }
        chomp($header = $_);
        $header =~ s/\s+$//;
        $seq = "";
#       print "Just read in this line: $_";
    
    } else {
        $_ =~ s/\s+$//;
        $seq .= $_;    #adds sequence to previous line
    }
}

$seq =~ s/^[^M]*//;
if(length $seq >= $minl){
    print $out "$header\n";
    print $out "$seq\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####

# EOF.
