#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_v);

# Usage
my $usage = "
vcfsummary.pl - summarizes vcf files.
		      by
		Brian J. Knaus
		  March 2013

Copyright (c) 2013 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl vcfsummary.pl options
 required:
  -a	vcf file.
 optional:
  -v	verbose mode [optional T/F, default is F].

";

# command line processing.
getopts('a:v:');
die $usage unless ($opt_a);

my ($inf, $verb);

$inf	= $opt_a if $opt_a;
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Globals.

my ($in, $outf, $out, $feat);
my @temp;
my @vkeys;
my %vars;
my %contigs;
my $i;

##### ##### ##### ##### #####
# Main.

# Manage outfile name.
$outf = outname($inf);


# Manage metadata.

open($in, "<", $inf) or die "Can't open $inf: $!";

while(<$in>){
	# Manage Metadata.
	if (/^##ALT=<ID=(.+),Description.+/){
		$vars{$1}=0;
	} elsif (/^#/) {
	
	} else {
		# Manage data.
		@temp = split (/\t/,$_);
		$temp[7]=~/.+SVTYPE=(.+)/;
		$feat = $1;
#		print $feat."\n";
		if (exists $contigs{$temp[0]}){
			# Contig exists.
			$contigs{$temp[0]}{$feat}=$contigs{$temp[0]}{$feat}+1;	
		} else {
			# New contig
#			print "$feat\t$temp[0]\n";
			$contigs{$temp[0]}={%vars};
			$contigs{$temp[0]}{$feat}=$contigs{$temp[0]}{$feat}+1;
		}
	}
}

close $in or die "$in: $!";


# Write hash of hashes to file.
@temp = keys(%contigs);
@temp = sort(@temp);
@vkeys = keys(%vars);
open($out, ">", $outf) or die "Can't open $outf: $!";

print $out "Contig\t".join("\t",@vkeys)."\n";

for $i	(0 .. $#temp){
  print $out "$temp[$i]\t";
  foreach(@vkeys){
#    print "$_\t";
    print $out $contigs{$temp[$i]}{$_}."\t";
  }
  print $out "\n";
}

close $out or die "$out: $!";

##### ##### ##### ##### #####
# Subroutines.

sub outname {
  my $inf = shift;
  my @temp = split("/", $inf);
  $inf = $temp[$#temp];
  @temp = split(/\./, $inf);
  my $outf = $temp[0];
  $outf = join("_", $outf, ".txt");

  return($outf);
}

##### ##### ##### ##### #####
# EOF.
