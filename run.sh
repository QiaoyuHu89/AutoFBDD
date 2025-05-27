#!/bin/bash

conda activate AutoFBDD
python AutoFBDD.py --mode linking --brick_mode brick --folder WDR5 --input_pdb WDR5.pdb --center center.txt --brickfolder all --brickfile_list bricks_file.txt --num_cpu 40 --sep_bricks 100 --top_clusters 10 --top_bricks 5 --dis_val 6 --poses 3

#python AutoPROTACs.py --folder WDR5 --target WDR5.pdb --E3_ligase CRBN.pdb --E3_ligand CRBN_ligand.mol2 --min_dist 5 --max_dist 15 --num_cpu 40 --top_docking 50

#conda activate DeepPROTACs
#python DeepPROTACs.py --folder WDR5 --E3_ligase CRBN.pdb --top_docking 50
