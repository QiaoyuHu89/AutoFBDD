#!usr/bin/perl
print "Thanks for using the chemical probes clustering tools developed by Fang Bai.et al."."\n";
print "USAGE: perl clustering.pl PROTEIN N Probes X1 Y1 X2 Y2 OUT2"."\n";
print "PROTEIN is the prefix of your .nrg file for protein and the protein name."."\n";
print "Present version require the prefixes of nrg file and pdb file are same, and the nrg file was named as PROTEIN_1.nrg, PROTEIN_2.nrg..."."\n";
print "N denotes the number of parts of the protein you have divided for mapping"."\n";
print "Probes denotes the probe name list files, such as our case file---probes_file.txt"."\n";
print "X1 and Y1 denote the required cluster density and the cluster radius for each probe respectively in the first round of clustering process."."\n";
print "Name of the output clusters obtained in the first round are automatically desinged by the script, you needn't indicate it."."\n";
print "X2 and Y2 denote the required cluster density and the cluster radius for each probe respectively in the first round of clustering process."."\n";
print "OUT2 denotes the prefix for the cluster name of obtained in the first round of clustering process."."\n";

print "OK, here we go~ Hope you enjoy!";
print "***********^_^************"."\n"."\n";



open(FILE,"$ARGV[2]");
my @probes = <FILE>;
my $length_probe = @probes;
close(FILE);

my $nrg;
print "######First Round Mapping#########"."\n";

open(DOC,"NSGA2.param")||die "can't open the NSGA_param file:$!";
my @doc_content = <DOC>;
close(DOC);

$doc_content[1] ="RUN_MODE      13"."\n"; 
$doc_content[51] ="CLUSTER_PARAMETER  ".$ARGV[3]."  ".$ARGV[4].".\n";
open(DOU,">NSGA2.param")||die "can't open the out NSGA_param file:$!";
print DOU @doc_content;
close(DOU);

for(my $j = 0; $j < $length_probe; $j++)
  { 
  	print $probes[$j]."\n";     
  	chomp($probes[$j]);
    my $lig = $probes[$j];

       $nrg = $ARGV[0]."_";
    
    open(DOCK,"dock.in")||die "can't open the out dock.in file:$!";
    my @dock_content = <DOCK>;
    close(DOCK);
      
    $dock_content[0] = "ligand_name"."     ".$nrg.$lig.".mol2"."\n";
    open(DOUK,">dock.in")||die "can't open the out dock.in file:$!";
    print DOUK @dock_content;
    close(DOUK);
    
    open(DOC,"NSGA2.param")||die "can't open the NSGA_param file:$!";
    my @doc_content = <DOC>;
    close(DOC);
    $doc_content[40] = "OUTPUT       ".$nrg."_cluster"."\n";
    open(DOU,">NSGA2.param")||die "can't open the out NSGA_param file:$!";
    print DOU @doc_content;
    close(DOU);
    
    system("./iFitDock -i dock.in ");
    }
	

print "The first round of clustering is finished."."\n";

print "######Second Round Mapping#########"."\n";

open(DOC,"NSGA2.param")||die "can't open the NSGA_param file:$!";
my @doc_content = <DOC>;
close(DOC);

$doc_content[1] ="RUN_MODE      14"."\n"; 
$doc_content[40] ="OUTPUT    ".$ARGV[7]."\n";
$doc_content[51] ="CLUSTER_PARAMETER  ".$ARGV[5]."  ".$ARGV[6].".\n";
open(DOU,">NSGA2.param")||die "can't open the out NSGA_param file:$!";
print DOU @doc_content;
close(DOU);

open(DOCK,"dock.in")||die "can't open the out dock.in file:$!";
my @dock_content = <DOCK>;
close(DOCK);
$dock_content[0] = "ligand_name"."     ".$nrg."_cluster_core".".mol2"."\n";
open(DOUK,">dock.in")||die "can't open the out dock.in file:$!";
print DOUK @dock_content;
close(DOUK);

system("./iFitDock -i dock.in ");

print "The second round of clustering is finished."."\n";
