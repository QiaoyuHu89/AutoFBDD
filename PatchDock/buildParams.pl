#!/usr/bin/perl -w

use strict;
use FindBin;

if ($#ARGV != 3 and $#ARGV != 2 and $#ARGV != 1) {
  print "buildParams.pl <receptorPDB> <ligandPDB> [rmsd default=4.0] [moleculesTypes]\n";
  print "The types are: EI:enzyme/inhibitor, AA:antibody/antigen, drug:protein/drug, default\n" ;
  exit;
}

my $recPdb = $ARGV[0];
my $ligPdb = $ARGV[1];
my $rmsd = -1.0;
if($#ARGV > 1) {
  $rmsd = $ARGV[2];
  if($rmsd <= 0.0 or $rmsd >= 10.0) {
    print "Invalid RMSD value, using default value!\n";
    $rmsd = -1.0;
  }
}
my $complexType= "Default";
if($#ARGV == 3) {
  $complexType = $ARGV[3];
  if($complexType ne "EI" and $complexType ne "AA" and $complexType ne "drug") {
    print "Invalid molecules type, using default!\n";
    $complexType= "Default";
  }
}

my $home = "$FindBin::Bin";

open OUT, ">params.txt";

print OUT "########################################################################\n";
print OUT "#    PatchDock Parameter File for $recPdb-$ligPdb";
print OUT "\n########################################################################\n";

print OUT "\n#    File Names:\n";

if($complexType eq "AA") {
  print OUT "receptorPdb $ligPdb\n";
  print OUT "ligandPdb $recPdb\n";
  print OUT "#receptorActiveSite epitopes.txt\n";
  print OUT "ligandActiveSite cdrs3\n";
  system("$home/cdr/cdr3.Linux $recPdb $home/cdr/ > cdrs3");
} else {
  print OUT "receptorPdb $recPdb\n";
  print OUT "ligandPdb $ligPdb\n";
  print OUT "#receptorActiveSite site1.txt\n";
  print OUT "#ligandActiveSite site2.txt\n";
  print OUT "#receptorBlockingSite site1.txt\n";
  print OUT "#ligandBlockingSite site2.txt\n";
}
print OUT "protLib $home/chem.lib\n";
print OUT "log-file patch_dock.log\n";
print OUT "log-level 2\n";

# distance constraint parameters
print OUT "\n#   Distance constraint parameters:\n";
print OUT "#       distanceConstraints <rec_atom_index> <lig_atom_index> <dist_thr>\n";
print OUT "# <rec_atom_index> - receptor atom used for constraint\n";
print OUT "# <lig_atom_index> - ligand atom used for constraint\n";
print OUT "# <dist_thr> - maximum allowed distance between the specified atom centers\n";
print OUT "#distanceConstraints rec_atom_index lig_atom_index dist_thr\n";

print OUT "#pointDistanceConstraints rec_coord lig_coord min_dist max_dist\n";
print OUT "#distanceConstraintsFile file_name\n";
print OUT "#membraneLine rec_atom_index lig_atom_index\n";

#segmentation params
print OUT "\n#    Surface Segmentation Parameters:\n";
print OUT "#       receptorSeg <low_patch_thr> <high_patch_thr> <prune_thr>\n";
print OUT "#                   <knob> <flat> <hole>\n";
print OUT "#                   <hot spot filter type>\n";
print OUT "#    <low_patch_thr>,<high_patch_thr> - min and max patch diameter\n";
print OUT "#    <prune_thr> - minimal distance between points inside the patch\n";
print OUT "#    <knob> <flat> <hole> - types of patches to dock (1-use, 0-do not use) (may need tuning)\n";
print OUT "#    <hot spot filter type> :None - 0, Antibody - 1, Antigen - 2\n";
print OUT "#                             Protease - 3, Inhibitor - 4, Drug - 5\n";
if($complexType =~ /drug/) {
  print OUT "receptorSeg 10.0 20.0 0.5 0 0 1 5\n";
  print OUT "ligandSeg 5.0 15.0 0.1 1 1 1 5\n";
} else {
  if($complexType eq "EI") {
    print OUT "receptorSeg 10.0 20.0 1.5 0 0 1 3\n";
    print OUT "ligandSeg 10.0 20.0 1.5 1 0 0 4\n";
  } else {
    if($complexType eq "AA") {
      print OUT "receptorSeg 10.0 20.0 1.5 1 0 1 2\n";
      print OUT "ligandSeg 10.0 20.0 1.5 1 0 1 1\n";
    } else {
      if(countAtoms($recPdb) > 5000 and countAtoms($ligPdb) > 5000) {
	print OUT "receptorSeg 10.0 20.0 2.0 1 0 1 0\n";
	print OUT "ligandSeg 10.0 20.0 2.0 1 0 1 0\n";
      } else {
	print OUT "receptorSeg 10.0 20.0 1.5 1 0 1 0\n";
	print OUT "ligandSeg 10.0 20.0 1.5 1 0 1 0\n";
      }
    }
  }
}

#scoring params
print OUT "\n#    Scoring Parameters:\n";
print OUT "#        scoreParams <small_interfaces_ratio> <max_penetration> <ns_thr>\n";
print OUT "#                    <rec_as_thr> <lig_as_thr> <patch_res_num> <w1 w2 w3 w4 w5>\n";
print OUT "#    <small_interfaces_ratio> - the ratio of the low scoring transforms to be removed\n";
print OUT "#    <max_penetration> - maximal allowed penetration between molecules surfaces\n";
print OUT "#    <ns_thr> - normal score threshold\n";
print OUT "#    <rec_as_thr> <lig_as_thr> - the minimal ratio of the active site area in the solutions\n";
print OUT "#    <patch_res_num> - number of results to consider in each patch\n";
print OUT "#    <w1 w2 w3 w4 w5> - scoring weights for ranges:\n";
print OUT "#                [-5.0,-3.6],[-3.6,-2.2],[-2.2,-1.0],[-1.0,1.0],[1.0-up]\n";

SWITCH: {
  if ($complexType eq "AA")
    { print OUT "scoreParams 0.3 -4.0 0.5 0.0 0.3 1500 -8 -4 0 1 0\n"; last SWITCH; }
  if ($complexType eq "EI")
    { print OUT "scoreParams 0.4 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0\n"; last SWITCH; }
  if ($complexType eq "Default")
    { print OUT "scoreParams 0.3 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0\n"; last SWITCH; }
  if ($complexType =~ /drug/)
    { print OUT "scoreParams 0.3 -4.0 0.4 0.0 0.0 200 -4 -2 0 1 0\n"; last SWITCH; }
}

#desolvation params
print OUT "\n#    Desolvation Scoring Parameters:\n";
print OUT "#        desolvationParams <energy_thr> <cut_off_ratio>\n";
print OUT "#    <energy_thr> - remove all results with desolvation energy higher than threshold\n";
print OUT "#    <cut_off_ratio> - the ratio of low energy results to be kept\n";
print OUT "#    First filtering with energy_thr is applied and the remaining results\n";
print OUT "#    can be further filtered with cut_off_ratio.\n";
print OUT "desolvationParams 500.0 1.0\n";

print OUT "\n########################################################################\n";
print OUT "#   Advanced Parameters";
print OUT "\n########################################################################\n\n";

#clustering params
print OUT "#    Clustering Parameters:\n";
print OUT "#    clusterParams < rotationVoxelSize > < discardClustersSmaller > < rmsd > < final clustering rmsd >\n";
if( $complexType =~ /drug/) {
  if( $rmsd == -1.0) {
    print OUT "clusterParams 0.05 2 1.0 1.5\n";
  } else {
    print OUT "clusterParams 0.05 2 1.0 $rmsd\n";
  }
} else {
  if($complexType eq "AA") {
    if( $rmsd == -1.0) {
      print OUT "clusterParams 0.1 3 2.0 4.0\n";
    } else {
      print OUT "clusterParams 0.1 3 2.0 $rmsd\n";
    }
  } else {
    if( $rmsd == -1.0) {
      print OUT "clusterParams 0.1 4 2.0 4.0\n";
    } else {
      print OUT "clusterParams 0.1 4 2.0 $rmsd\n";
    }
  }
}

#base params
print OUT "\n#    Base Parameters:\n";
print OUT "#    baseParams <min_base_dist> <max_base_dist> <# of patches for base: 1 or 2>\n";
if ( $complexType =~ /drug/ ) {
  print OUT "baseParams 1.0 10.0 1\n";
} else {
  print OUT "baseParams 4.0 13.0 2\n";
}

#matching params
print OUT "\n#    Matching Parameters:\n";
print OUT "#  matchingParams <geo_dist_thr> <dist_thr> <angle_thr> <torsion_thr> \n";
print OUT "#     <angle_sum_thr>\n";
if($complexType =~ /drug/ ) {
  print OUT "matchingParams 2.0 1.0 0.4 0.5 0.9\n";
} else {
  if($complexType eq "AA" ) {
    print OUT "matchingParams 1.5 1.5 0.4 0.5 0.9\n";
  } else {
    print OUT "matchingParams 1.5 1.5 0.4 0.5 0.9\n";
  }
}

print OUT "# 1 - PoseClustering (default), 2 - Geometring Hashing\n";
print OUT "matchAlgorithm 1\n";

#grid params
print OUT "\n#    Grid Parameters:\n";
print OUT "#      receptorGrid <gridStep> <maxDistInDistFunction> <vol_func_radius>\n";
print OUT "#      NOTE: the vol_func_radius of small molecules and peptides should be 3A!\n";
print OUT "receptorGrid 0.5 6.0 6.0\n";
if ( $complexType =~ /drug/ ) {
  print OUT "ligandGrid 0.35 6.0 3.0\n";
} else {
  print OUT "ligandGrid 0.5 6.0 6.0\n";
}

#ms params
print OUT "\n#    Connolly Surface Parameters:\n";
print OUT "#      receptorMs <density> <probe_radius>\n";
print OUT "receptorMs 10.0 1.8\n";
if ( $complexType =~ /drug/ ) {
  print OUT "ligandMs 20.0 1.4\n";
} else {
  print OUT "ligandMs 10.0 1.8\n";
}


close (OUT);

sub countAtoms {
  my $filename = shift;
  open FILE, "<$filename" or die "Can't open file: $filename";
  my @atomLines =  grep /(^ATOM|^HETATM)/, <FILE>;
  my $atomsNumber = grep ! /HOH/, @atomLines;
  close FILE;
  return $atomsNumber;
}
