#!/usr/bin/perl -w

use strict;
use FindBin;

if ($#ARGV != 1 and $#ARGV !=2 and $#ARGV !=3) {
  print "transOutput.pl <output file name> <first result> [last result] [only ligand flag]\n";
  exit;
}

my $ligandPdb = "";
my $receptorPdb = "";

my $resFileName = $ARGV[0];
my $first = $ARGV[1];
my $last = $first;

if ($#ARGV > 1) {
  $last = $ARGV[2];
}

open(DATA, $resFileName);

my $home = "$FindBin::Bin";

while(<DATA>) {
  chomp;
  my @tmp=split('\|',$_);
  if($#tmp>0 and $tmp[0] =~/\d/) {
    my $transNum=int $tmp[0];
    if($transNum >= $first and $transNum <= $last and length $ligandPdb > 0 and length $receptorPdb > 0) {
      # apply transformation on the ligand molecule
      my $currResFile = "$resFileName.$transNum.pdb";
      unlink $currResFile;
      if ($#ARGV < 3) { #add receptor too
	`cat $receptorPdb | grep -v '^END' > $currResFile`;
      }
      print "$home/pdb_trans $tmp[$#tmp] < $ligandPdb >> $currResFile\n";
      `$home/pdb_trans $tmp[$#tmp] < $ligandPdb >> $currResFile`;
    }
  } else {
    # find the filenames in the output file
    @tmp=split(' ',$_);
    if($#tmp>0) {
      if($tmp[0] =~ /ligandPdb/) {
	$ligandPdb = $tmp[2];
	print "Ligand PDB: $ligandPdb\n";
      }
      if($tmp[0] =~ /receptorPdb/) {
	$receptorPdb = $tmp[2];
	print "Receptor PDB: $receptorPdb\n";
      }
    }
  }
}

