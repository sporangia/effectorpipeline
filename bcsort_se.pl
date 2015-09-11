#!/usr/bin/perl

use strict;
use warnings;
#use Cwd;
use Cwd 'abs_path';
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use File::Spec;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_c $opt_d $opt_f $opt_g $opt_s $opt_m $opt_v);

# Usage
my $usage = "
bcsort_se.pl v0.1-2 
Sort barcodes from a pair of fastq files.
		      by
		Brian J. Knaus
		 October 2011

Copyright (c) 2010, 2011 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl bcsort_pe.pl options
 required:
  -a	fastq input file a.
 optional:
  -c	barcode file name, no quotes, one barcode per line
	[optional, default = barcodes.txt].
  -d	output directory [optional, default is bcsort_out].
  -f	fastq file is gzipped [default is T].
  -g	filter reads [default is T].
  -s	output fasta file [default is F].
  -m	output fastq file [default is T].
  -v	verbose mode [optional T/F, default is F].

";

# command line processing.
getopts('a:c:d:f:g:s:m:n:v:');
die $usage unless ($opt_a);

my ($infa, $bcfile, $outdir, $zip, $filter, $fa, $fq, $verb);

$infa	= $opt_a if $opt_a;

$bcfile	= $opt_c ? $opt_c : "barcodes.txt"; 
#$outdir	= $opt_d ? "bcsort_out".$opt_d :"bcsort_out";
$outdir	= $opt_d ? $opt_d :"bcsort_out";
$zip	= $opt_f ? $opt_f : "T";
$filter	= $opt_g ? $opt_g : "T";
$fa	= $opt_s ? $opt_s : "F";
$fq	= $opt_m ? $opt_m : "T";
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Globals.

my $bclog = "00_bcsort_log.txt";
my ($in, $out, $out2, @temp, @temp2);
my (%barcodesH, %bcnames, @barcodesA);
#my ($anom1, $anom2);

my $ina;
my $anoma = "anomalous_reads1.fq";
my $anom1;

my ($ida, $bca, $seqa, $prba);
my ($ofile1, $ofile2);

my $filtered = 0;
my $anom_num = 0;
my $bc_num = 0;
my $ntot = 0;
my %samps;

##### ##### ##### ##### #####
# Main.

# Get absolute path of infiles.
$infa = abs_path($infa);
$bcfile = abs_path($bcfile);
$outdir = abs_path($outdir);

# Create logfile name.
$bclog = File::Spec->catfile($outdir,$bclog);

# Create anomaly file names.
$anoma = File::Spec->catfile($outdir,$anoma);

# Check infile is a valid fastq.
fastqchk($infa, $zip);

# Initialize log.
if(!-e $outdir ){mkdir $outdir;}
open($out, ">", $bclog) or die "Can't open $bclog: $!";
print $out $usage;
print $out "\n      ----- ----- --*-- ----- -----\n\n";
print $out "Infile a:\t$infa\n";
print $out "Barcodes file:\t$bcfile\n";
print $out "Out directory:\t$outdir\n";
print $out "fastq file is gzipped:\t$zip\n";
print $out "filter reads\t$filter\n";
print $out "Print fasta:\t$fa\n";
print $out "Print fastq:\t$fq\n";
print $out "Verbose:\t$verb\n";
print $out "\n      ----- ----- --*-- ----- -----\n\n";
close $out or die "$out: $!";

# Initial timestmp.
my @stime = localtime(time);
timestamp($bclog, 'Process started');

##### ##### ##### ##### #####
# Input barcodes.
timestamp($bclog, 'Processing barcodes');
open($in,  "<",  $bcfile)  or die "Can't open $bcfile: $!";

while (<$in>){
  chomp;
  @temp = split("\t", $_);
  $barcodesH{uc($temp[0])} = 0;
  $bcnames{uc($temp[0])} = $temp[1];
  push @barcodesA, uc($temp[0]);
}
@barcodesA = sort(@barcodesA);

open($out, ">>", $bclog) or die "Can't open $bclog: $!";
print $out "Barcodes:\n\n";
foreach(@barcodesA){
  print $out "$_\t=>\t";
  print $out $bcnames{$_}."\n";
}
print $out "\n      ----- ----- --*-- ----- -----\n\n";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# Input data and sort
# into a hash of arrays.
#
# Keys are barcodes.
# Values are arrays of samples.

timestamp($bclog, 'Reading in data');

if($zip eq 'T'){
  $ina = new IO::Uncompress::Gunzip ($infa, MultiStream=>1) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  $anoma = $anoma.".gz";
  $anom1 = new IO::Compress::Gzip $anoma or die "IO::Compress::Gzip failed: $GzipError\n";
} else {
#    open($in, "<", $inf) or die "Can't open $inf: $!";
  open($ina,  "<",  $infa)  or die "Can't open $infa: $!";
  open($anom1,  ">",  $anoma)  or die "Can't open $anoma: $!";
}

# Make regex.
my $test = join("|",@barcodesA);

# Read in.
while (<$ina>){
  $ntot = $ntot + 1;

  # Read a.
  chomp($ida	= $_);
  $ida = substr $ida,1;
  chomp($seqa	= <$ina>);
  <$ina>;
  chomp($prba	= <$ina>);

  # Filter for quality.
  if($filter eq 'T'){
    @temp = split(" ",$ida);
    @temp = split(":",$temp[1]);
    if($temp[1] eq 'Y'){
      $filtered++;
      # Read a.
      print $anom1 "\@$ida\n";
      print $anom1 "$seqa\n";
      print $anom1 "+$ida\n";
      print $anom1 "$prba\n";
      next;
    }
  }

  my @matches1 = ($seqa =~ /^($test)/o);
  if ($#matches1 == 0){
    # Sequence one matches a barcode.
    @temp = debarcode($matches1[0],$ida,$seqa,$prba);
    $ida = $temp[0];
    $seqa = $temp[1];
    $prba = $temp[2];
    push @{$samps{$matches1[0]}}, join("\t", $ida, $seqa, $prba); # Hash of arrays.
    $barcodesH{$matches1[0]} = $barcodesH{$matches1[0]} + 1;
    $bc_num = $bc_num + 1;
  } else {
    # Anomalous reads.
    $anom_num = $anom_num + 1;
    # Read a.
    print $anom1 "\@$ida\n";
    print $anom1 "$seqa\n";
    print $anom1 "+$ida\n";
    print $anom1 "$prba\n";
  }
  if ($ntot % 1000000 == 0){
    open($out, ">>", $bclog) or die "Can't open $bclog: $!";
    print $out $ntot/1000000, "\tmillion records sorted.\n";
    close $out or die "$out: $!";
  }
}

if($zip eq 'T'){
  close $ina or die "$ina: $GunzipError\n";
  close $anom1 or die "$anom1: $GzipError\n";
} else {
  close $ina or die "$ina: $!";
  close $anom1 or die "$anom1: $!";
}

# Report to log.
open($out, ">>", $bclog) or die "Can't open $bclog: $!";
print $out "\n";
close $out or die "$out: $!";
timestamp($bclog, 'Data read and sorted');

##### ##### ##### ##### #####
# Write data to file.

# Report to log.
timestamp($bclog, 'Writing data to file');

# fasta files.
if ($fa eq "T"){
  foreach $bca (keys %samps) {
    $ofile1 = $bca."_seq.fa";
    $ofile1 = File::Spec->catfile($outdir,$ofile1);

    open($out, ">",  $ofile1) or die "Can't open $ofile1: $!";
    foreach ( @{$samps{$bca}}){
      @temp = split("\t",$_);
      print $out ">$temp[0]\n";
      print $out "$temp[1]\n";
    }
    close $out or die "$out: $!";
  }
  # Report to log.
  timestamp($bclog, 'Fasta files written');
}

# Fastq files.
if ($fq eq "T"){
  foreach $bca (keys %samps) {
    if($zip eq 'T'){
      $ofile1 = $bca."_seq1.fastq.gz";
      $ofile1 = File::Spec->catfile($outdir,$ofile1);
      $out  = new IO::Compress::Gzip $ofile1 or die "IO::Compress::Gzip failed: $GzipError\n";
    } else {
      $ofile1 = $bca."_seq1.fastq";
      $ofile1 = File::Spec->catfile($outdir,$ofile1);
      open($out, ">",  $ofile1) or die "Can't open $ofile1: $!";
    }

    foreach ( @{$samps{$bca}}){
      @temp = split("\t",$_);
      print $out "\@$temp[0]\n";
      print $out $temp[1], "\n";
#      print $out "+$temp[0]\n";
      print $out "+\n";
      print $out "$temp[2]\n";
    }
    if($zip eq 'T'){
      close $out or die "$ina: $GunzipError\n";
    } else {
      close $out or die "out: $!";
    }
  }
  # Report to log.
  timestamp($bclog, 'Fastq files written');
}

# Report to log.
timestamp($bclog, 'Data written to file');

##### ##### ##### ##### #####
# Summary.

open($out, ">>", $bclog) or die "Can't open $bclog: $!";

#print $out "\n      ----- ----- --*-- ----- -----\n\n";
print $out "\nTotal reads:\t\t";
print $out $ntot / 1000000, "\tmillion";
print $out "\nFiltered reads:\t\t";
print $out $filtered / 1000000, "\tmillion";
print $out "\nBarcoded reads:\t\t";
print $out $bc_num / 1000000, "\tmillion";
print $out "\nNon-barcoded reads:\t";
#print $out ($ntot - $bc_num) / 1000000, "\tmillion (total-barcoded)";
print $out "\nAnomalous reads:\t";
print $out $anom_num / 1000000, "\tmillion";

#print $out "\n\n";
#print $out "\nTotal reads:\t\t";
#print $out $ntot , "\treads";
#print $out "\nFiltered reads:\t\t";
#print $out $filtered , "\treads";
#print $out "\nBarcoded reads:\t\t";
#print $out $bc_num , "\treads";
#print $out "\nNon-barcoded reads:\t";
#print $out "\nAnomalous reads:\t";
#print $out $anom_num , "\treads";

print $out "\n\n";

print $out "Barcode Summary:\n\n";
#while (my($key, $value) = each %barcodes){
#	print $log "$key \t=>\t", $value/1000000, "\tmillion", "\n";
#}
foreach (@barcodesA){
  print $out "$_\t=>\t".$barcodesH{$_}."\t=>\t".$bcnames{$_}."\n";
}
#print $out "\n      ----- ----- --*-- ----- -----\n\n";

close $out or die "$out: $!";

##### ##### ##### ##### #####
# Finish.

open($out, ">>", $bclog) or die "Can't open $bclog: $!";
print $out "\n      ----- ----- --*-- ----- -----\n\n";
print $out "Process initiated:\t";
print $out $stime[4]+1;
print $out "-$stime[3]-";
print $out $stime[5]+1900;
#print $out $year+1900;
print $out " $stime[2]:$stime[1]:$stime[0]\n";
close $out or die "$out: $!";

# Log finish time.
timestamp($bclog, 'Process finished');

open($out, ">>", $bclog) or die "Can't open $bclog: $!";
print $out "bcsort process complete.\n\n";
close $out or die "$out: $!";

my $ffile = File::Spec->catfile($outdir,"bcsort.finished");
open($out, ">>", $ffile) or die "Can't open $ffile: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# Subfunctions.             #
##### ##### ##### ##### #####

sub idbc {
	# id tab bc as input.
	my @id = split(/\#|\//, shift);
	my $bca = shift;
	$id[1] = $bca;
	return join("", $id[0], "\#", $id[1], "/", $id[2]);
}

# Check fastq file.
sub fastqchk{
  my $inf = shift;
  my $zip = shift;
  my @temp;
  my $in;
  # Check that file exists.
  if(-e $inf){} else {
    die "Error: $inf does not exist!\n";
  }

  if($zip eq 'T'){
    $in = new IO::Uncompress::Gunzip $inf or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  } else {
    open($in, "<", $inf) or die "Can't open $inf: $!";
  }
#  open($in, "<", $inf) or die "Can't open $inf: $!";
  chomp($temp[0] = <$in>);	# First line is id.
  chomp($temp[1] = <$in>);	# Second line is sequence.
  chomp($temp[2] = <$in>);	# Third line is id.
  chomp($temp[3] = <$in>);	# Fourth line is quality.
  # Check identifier line.
  if($temp[0] !~ /^@/){
    print "Unexpected file format!\n";
    print "Line1: $temp[0]\n";
    print "Line2: $temp[1]\n";
    print "Line3: $temp[2]\n";
    print "Line4: $temp[3]\n";
    die "Error: Fastq file expected.\n";
  }
  # Check sequence.
  $_ = $temp[1];
  s/A//gi;
  s/C//gi;
  s/T//gi;
  s/G//gi;
  s/N//gi;
  if(length($_)>0){
    print "Unexpected file format!\n";
    print "Line1: $temp[0]\n";
    print "Line2: $temp[1]\n";
    print "Line3: $temp[2]\n";
    print "Line4: $temp[3]\n";
    die "Invalid character in sequence: $_.\n";
  }
  # Check identifier line.
  if($temp[2] !~ /^\+/){
    print "Unexpected file format!\n";
    print "Line1: $temp[0]\n";
    print "Line2: $temp[1]\n";
    print "Line3: $temp[2]\n";
    print "Line4: $temp[3]\n";
    die "Error: Fastq file expected.\n";
  }
  if($zip eq 'T'){
    close $in or die "$in: $GunzipError\n";
  } else {
    close $in or die "$in: $!";
  }
}

sub timestamp{
  my $logf = shift;
  my $msg = shift;
  my $out;

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

  open($out, ">>", $logf) or die "Can't open $logf: $!";
  print $out $msg.":\t";
  print $out $mon+1;
  print $out "-".$mday."-";
  print $out $year+1900;
  print $out " $hour:$min:$sec\n";
  print $out "\n      ----- ----- --*-- ----- -----\n\n";
  close $out or die "$out: $!";
}

sub debarcode{
  my $bc = shift;
  my $ida = shift;
  my $seqa = shift;
  my $prba = shift;

  my $q = substr($prba, 0, length($bc));
  $seqa = substr($seqa, length($bc));
  $prba = substr($prba, length($bc));
#  $ida = $ida.$bc.':'.$q; # Not sure this is standard.
  $ida = $ida.$bc;
  return($ida,$seqa,$prba);
}

sub debarcode_v1{
  my $bcode = shift;
  my $id = shift;
  my $seq = shift;
  my $prb = shift;
  my $mach;
  my @temp;

  my $bclength = length($bcode);
  $seq = substr $seq, $bclength;
  $prb = substr $prb, $bclength;
  @temp = split /[\/#]/, $id;
  $id = "$temp[0]#$bcode/$temp[2]";

  return($id, $seq, $prb);
#      @temp = debarcode($_,$ida,$seqa,$prba);

}


##### ##### ##### ##### #####
# EOF.
