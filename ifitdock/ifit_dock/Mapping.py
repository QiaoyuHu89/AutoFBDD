# usage : python Mapping.py PROTEIN N Probes CENTER
import sys
import os

dock_content = [
  "ligand_name                      1OWE__cluster_core.mol2\n",
  "atom_model                       all\n",
  "grid_score_rep_rad_scale         1\n",
  "grid_score_vdw_scale             1\n",
  "grid_score_es_scale              1\n",
  "grid_score_grid_prefix           target_1\n",
  "vdw_defn_file                    vdw.defn\n",
  "flex_defn_file                   flex.defn\n",
  "flex_drive_file                  flex_drive.tbl\n",
  "chem_defn_file                   chem.defn\n",
  "gbsa_zou_gb_grid_prefix          target_GB_1\n",
  "gbsa_zou_sa_grid_prefix          target_SA_1\n",
  "gbsa_zou_vdw_grid_prefix         target_GB_1\n",
  "gbsa_zou_screen_file             screen.in\n",
  "gbsa_zou_solvent_dielectric      78.300003"
]
doc_content = [
  "#Running mode choice:(0~12)\n",
  "RUN_MODE      14\n",
  "#Number of objective functions\n",
  "NUM_OBJ_FUNC    2\n",
  "STRATEGY        2\n",
  "#Choice for optimization method % NSGA2 or NSGA2&MOPSO\n",
  "OPTIMIZATION_METHOD   NSGA2\n",
  "#Number of constraint solutions\n",
  "NUM_CONS    0\n",
  "#Population size\n",
  "POP_SIZE     2000\n",
  "#Possibility of the cross over\n",
  "PCROSS     0.9\n",
  "#Possibility of the mutation\n",
  "PMUT      0.3\n",
  "#SBX varible\n",
  "ETA_C  10\n",
  "#SBX varible\n",
  "ETA_M   10\n",
  "#Generation of NSGA\n",
  "GEN    200\n",
  "#Starting mutation point of variables\n",
  "NREAL_MUTE  0.2\n",
  "#Staring cross point of variables\n",
  "NREAL_CROSS  0.8\n",
  "#Lower limit for transition\n",
  "MIN_REAL_TRANS        -2.0  -2.0  -2.0  -3.142  -3.142   -3.142\n",
  "#Upper limitfor transition\n",
  "MAX_REAL_TRANS           2.0   2.0  2.0   3.142  3.142  3.142\n",
  "#Range for rotation\n",
  "RANGE_OF_ROTATE         -3.142     3.142\n",
  "#Maximum size of the non-dominated soultions\n",
  "MAX_FEASIBLE        1000\n",
  "#Maximum speed of the partical\n",
  "MAX_VELOCITY           0.2\n",
  "#Minimum speed of the partical\n",
  "MIN_VELOCITY     1.0\n",
  "ACTIVE_SITE               8.929000 -2.000000 0.731000\n",
  "EPSILON      0.1   0.1    0.1    0.1\n",
  "ENERGY_WEIGHTS      1.0  1.0  1.0  1.0\n",
  "OUTPUT    cluster_out\n",
  "RECEPTOR      target_noH.pdb\n",
  "#Number of conformations per ligand required to output\n",
  "NUM_CONF_PER_LIG  1\n",
  "#Score cutoff for determing the output of docked ligands\n",
  "SCORE_CUTOFF   -1.0\n",
  "INDUCEFIT    1\n",
  "METAL_SWITCH   0\n",
  "METAL_FILE   metal.txt\n",
  "LIG_ELE_CHOICE   1\n",
  "CHOICE_FOR_SOLVENT 0 1B_H2O_opti.mol2\n",
  "CLUSTER_PARAMETER  3  3."
]

with open(sys.argv[3]) as f:
  probes = f.readlines()
probes = [i.strip() for i in probes]

len_probes = len(probes)
for j in range(len_probes):
  lig = probes[j]
  for k in range(int(sys.argv[2])):
    nrg = sys.argv[1]
    label = k+1
    
    dock_content[0] = f"ligand_name        ./file/{lig}\n"
    dock_content[5] = f"grid_score_grid_prefix    {nrg}_{label}\n"
    dock_content[10] = f"gbsa_zou_gb_grid_prefix     {nrg}_GB_{label}\n"
    dock_content[11] = f"gbsa_zou_sa_grid_prefix     {nrg}_SA_{label}\n"
    dock_content[12] =f"gbsa_zou_vdw_grid_prefix     {nrg}_GB_{label}\n"
    with open("dock.in",'w') as f:
      f.writelines(dock_content)
        
    doc_content[1] ="RUN_MODE      0\n"
    doc_content[3] ="NUM_OBJ_FUNC    2\n"
    doc_content[4] ="STRATEGY        2\n"
    doc_content[40] = f"OUTPUT       {nrg}_{lig}\n"
    doc_content[41] = f"RECEPTOR      {nrg}_noH.pdb\n"
      
    with open("NSGA2.param","w") as f:
      f.writelines(doc_content)
    
    with open(sys.argv[4]+"_"+str(label)+".txt") as f:
      centers = f.readlines()
    len_center = len(centers)
    
    for i in range(len_center):
        line = centers[i].strip()
        h = i + 1 
        doc_content[37] = f"ACTIVE_SITE               {line}\n"
        with open("NSGA2.param","w") as f:
          f.writelines(doc_content)
        for c in doc_content:
          print(c)
        os.system("./iFitDock -i dock.in")

print("Yep! We finished the mapping process and let's move to clustering process.\n")
print("Please use script cluster.pl and follow the intruction of the cluster.pl\n")
