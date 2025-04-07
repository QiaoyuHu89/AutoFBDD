#!usr/bin/perl
print "Thanks for using our energy grid generation tools."."\n";
print "USAGE: perl mapping_nrg_grid_generator.pl PROTEIN N CENTER"."\n";
print "PROTEIN is the prefix of your .nrg file for protein and the protein name."."\n";
print "Present version require the prefixes of nrg file and pdb file are same, and the nrg file was named as PROTEIN_1.nrg, PROTEIN_2.nrg..."."\n";
print "N denotes the number of parts of the protein you have divided for mapping"."\n";
print "CENTER denotes prefix for the center list file for each part of the divided protein, such as prefix 'Center'in our example files named as Center_1.txt, Center_2.txt...; suffix fixed as '.txt' at present version"."\n";
print "OK, here we go~, this process is relatively time consuming. Being patient. Hope you enjoy!";
print "***********^_^************"."\n"."\n";

system"rm OUT*";
system"rm temp*";

   open(CENTER,"$ARGV[2]")||die"can't open the CENTER file list:$!";
   my @center_list = <CENTER>;
   my $length_center = @center_list;
   close(CENTER);
   for($j =0 ,$k = 1;$j < $ARGV[1]; $j++,$k++)
   {
   	chomp($center_list[$j]);
   	my $center = $center_list[$j]."\n";
   	open(BOX,"box.in")||die "can't open the out box file:$!";
   	my @box_content = <BOX>;
   	close(BOX);
   	my $vdw = "_grid_vdw_box"."_".$k.".pdb"."\n";
   	$box_content[2] = $center;
   	$box_content[3] ="30.0  30.0  30.0"."\n";
   	$box_content[4] =$ARGV[0].$vdw;
    open(BOX,">box.in")||die "can't open the out boxin_vdw file:$!";
    print BOX @box_content;
    close(BOX); 
    system"showbox < box.in";
    
    open(BOX,"box.in")||die "can't open the out box file:$!";
   	my @box_content = <BOX>;
   	close(BOX);
   	my $gbsa = "_grid_gbsa_box"."_".$k.".pdb"."\n";
   	chomp($box_content[3]);
   #	$box_content[3]=();
   	chomp($box_content[4]);
   	#$box_content[4]=();
   	$box_content[3] ="20.0  20.0  20.0"."\n";
   	$box_content[4] =$ARGV[0].$gbsa;
    open(BOX,">box.in")||die "can't open the out boxin_vdw file:$!";
    print BOX @box_content;
    close(BOX);
    system"showbox < box.in";
    print"BOX is generated"."\n";
    
    open(GRD,"grid.in")||die "can't open the out grid file:$!";
   	my @grd_content = <GRD>;
   	close(GRD);
   	my $mol_name = $ARGV[0].".mol2"."\n";
   	my $grd_name = $ARGV[0]."_grid_".$k."\n";
   	chomp($grd_content[14]);
   	$grd_content[14]=();
   	$grd_content[15]=();
   	$grd_content[17]=();
   	
   	$grd_content[14] = "receptor_file        ".$mol_name;
   	$grd_content[15] = "box_file       ".$ARGV[0].$vdw;
   	$grd_content[17] = "score_grid_prefix        ".$ARGV[0]."_".$k;
   	open(GRD,">grid.in")||die "can't open the out grdin_vdw file:$!";
    print GRD @grd_content;
    close(GRD);
    system"grid -i grid.in";
   	
   }
   print"Energy Grids are generated"."\n"; 
   system"rm temp*";
   system"rm OUT*";


