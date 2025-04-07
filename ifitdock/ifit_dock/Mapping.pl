#!usr/bin/perl
print "Thank you for using the druggble interface identifier develped by Fang Bai.et al."."\n";
print "USAGE: perl Mapping.pl PROTEIN N Probes CENTER"."\n";
print "PROTEIN is the prefix of your .nrg file for protein and the protein name."."\n";
print "Present version require the prefixes of nrg file and pdb file are same, and the nrg file was named as PROTEIN_1.nrg, PROTEIN_2.nrg..."."\n";
print "N denotes the number of parts of the protein you have divided for mapping"."\n";
print "Probes denotes the probe name list files, such as our case file---probes_file.txt"."\n";
print "CENTER denotes prefix for the center list file for each part of the divided protein, such as prefix 'Center'in our example files named as Center_1.txt, Center_2.txt...; suffix fixed as '.txt' at present version"."\n";
print "OK, here we go~ Hope you enjoy!";
print "***********^_^************"."\n"."\n";


open(FILE,"$ARGV[2]")||die "can't open the probe name list file.";#probes_file.txt
my @probes = <FILE>;
my $length_probe = @probes;
close(FILE);


for(my $j = 0; $j < $length_probe; $j++)
  { 
  	print "We are mapping the probe ".$probes[$j]." to the each part of the protein one by one."."\n";     
  	chomp($probes[$j]);
    my $lig = $probes[$j];
    
    for(my $k = 0; $k < $ARGV[1]; $k++)
    {
      my $nrg =  $ARGV[0];
     # print "Computating the ";
     # print $k+1;
     # print "part of the protein"."\n";
      my $lable = $k+1;
      open(DOCK,"dock.in")||die "can't open the out dock.in file:$!";
      my @dock_content = <DOCK>;
      close(DOCK);
      
      $dock_content[0] = "ligand_name"."     "."./file/".$lig."\n";
      $dock_content[5] = "grid_score_grid_prefix"."     ".$nrg."_".$lable."\n";
      $dock_content[10] ="gbsa_zou_gb_grid_prefix"."     ".$nrg."_GB_".$lable."\n";
      $dock_content[11] ="gbsa_zou_sa_grid_prefix"."     ".$nrg."_SA_".$lable."\n";
      $dock_content[12] ="gbsa_zou_vdw_grid_prefix"."     ".$nrg."_GB_".$lable."\n";
      open(DOUK,">dock.in")||die "can't open the out dock.in file:$!";
      print DOUK @dock_content;
      close(DOUK);
      
      
      open(DOC,"NSGA2.param")||die "can't open the NSGA_param file:$!";
      my @doc_content = <DOC>;
      close(DOC);
      $doc_content[1] ="RUN_MODE      0"."\n"; 
      $doc_content[3] ="NUM_OBJ_FUNC    2"."\n";
      $doc_content[4] ="STRATEGY        2"."\n";
      $doc_content[40] = "OUTPUT       ".$nrg."_".$lig."\n";
      $doc_content[41] = "RECEPTOR      ".$nrg."_noH.pdb"."\n";
      open(DOU,">NSGA2.param")||die "can't open the out NSGA_param file:$!";
      print DOU @doc_content;
      close(DOU);
    
      my $center_file = $ARGV[3]."_".$lable.".txt";
       print "#######################################2\n";
      print $center_file;
      open(FILE,"$center_file");
      my @list = <FILE>;
      my $length_file = @list;
      close(FILE);
      print "#######################################1\n";
      print $length_file;
      for(my $i = 0; $i < $length_file; $i++)
      {

      	my $line = $list[$i];
	      chomp($line);
        print $line."\n";
        my $h =$i + 1; 
        open(DOC,"NSGA2.param")||die "can't open the out box file:$!";
        my @doc_content = <DOC>;
        close(DOC);
        $doc_content[37] = "ACTIVE_SITE               ".$line."\n";
        open(DOU,">NSGA2.param")||die "can't open the out grdin_vdw file:$!";
        print DOU @doc_content;
        close(DOU);
        
        print ("#######################################1\n");
        system("./iFitDock -i dock.in ");
      }
    }
	}
	print "Yep! We finished the mapping process and let's move to clustering process."."\n";
	print "Please use script cluster.pl and follow the intruction of the cluster.pl"."\n";
	





