#!/usr/bin/perl -w

use strict;
use FindBin;

if ($#ARGV != 0)
{
  print "runCDR.pl <antibody_PDB_file>\n" ;
  exit ;
}

my $home = "$FindBin::Bin";

system("$home/cdr/cdr.Linux $ARGV[0] $home/cdr/ > cdrs");
system("$home/cdr/cdr3.Linux $ARGV[0] $home/cdr/ > cdrs3");

print "CDRs are computed in the 'cdrs' and 'cdrs3' files\n";



