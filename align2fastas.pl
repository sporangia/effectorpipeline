#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_m $opt_v);

# Usage
my $usage = "
align2fastas.pl v0.1-1 
seperates a file into sequences.
		      by
		Brian J. Knaus
		  August 2011

Copyright (c) 2011 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl align2fastas.pl options
 required:
  -a	fasta file.
 optional:
  -m	use SGE [T/F, default is T].
  -v	verbose mode [optional T/F, default is F].

";

# command line processing.
getopts('a:m:v:');
die $usage unless ($opt_a);

my ($inf, $sge, $verb);

$inf	= $opt_a if $opt_a;
$sge	= $opt_m ? $opt_m : "T";
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Globals.

my ($in, $outf, $out);
my ($id, $seq);
my $i;

#my @temp;

##### ##### ##### ##### #####
# Main.

$inf = abs_path($inf);

mkdir "fastas";
chdir "fastas";

$i=1;
open($in,  "<",  $inf)  or die "Can't open $inf: $!";

chomp($id = <$in>);
$outf = join "", substr($id, 1), ".fa";
#$outf = join "", $i, ".fa";
$i++;
open($out, ">", $outf) or die "Can't open $outf: $!";
print $out "$id\n";

while(<$in>){
  chomp;
  if (/^>/){
    # New record.
    close $out or die "$out: $!";
    $outf = join "", substr($_, 1), ".fa";
#    $outf = join "", $i, ".fa";
    $i++;
    open($out, ">", $outf) or die "Can't open $outf: $!";
    print $out "$_\n";
  } else {
    print $out "$_";
  }
}

close $in or die "$in: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# EOF.
