#!usr/bin/perl
print "Thanks for using our SA_grid generation tools."."\n";
print "USAGE: perl mapping_SA_grid_generator.pl PROTEIN N CENTER"."\n";
print "PROTEIN is the prefix of your SA grids for protein and the protein name, protein with the partial hydrogen"."\n";
print "Present version require the prefixes of nrg file and pdb file are same, and the nrg file was named as PROTEIN_1.nrg, PROTEIN_2.nrg..."."\n";
print "N denotes the number of parts of the protein you have divided for mapping"."\n";
print "CENTER denotes the file name for the center list file for each part of the divided protein"."\n";
print "OK, here we go~ Hope you enjoy!";
print "***********^_^************"."\n"."\n";

 
system"rm OUT*";
system"rm temp*";

my $center_name = $ARGV[2];
open(CENTER,"$center_name")||die"can't open the CENTER file list:$!";
my @center_list = <CENTER>;
my $length_center = @center_list;
close(CENTER);
for($j =0 ;$j < $ARGV[1]; $j++)
{ 
	my $k = $j +1;
  my $gbsa_pdb = $ARGV[0]."_gbsa.pdb"."\n";
  my $box = $ARGV[0]."_grid_gbsa_box_".$k.".pdb"."\n";
  my $sa =$ARGV[0]."_SA_".$k."\n";
   	 
  open(GBFILE,"INCHEM")||die "can't open the out box file:$!";
  my @gb_content = <GBFILE>;
  close(GBFILE);
   	
  $gb_content[0] = $gbsa_pdb ;
  $gb_content[3] = $box;
  $gb_content[8] = $sa;
  open(GRD,">INCHEM")||die "can't open the out grdin_vdw file:$!";
  print GRD @gb_content;
  close(GRD);
    
  system("nchemgrid_SA");
  } 
  print "SA Grids are generated successfully."."\n";

