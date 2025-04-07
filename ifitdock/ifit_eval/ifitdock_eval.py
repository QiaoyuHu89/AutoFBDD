import sys
import os
import time
import shutil
import argparse
from ifitdock_func import *
from multiprocessing import Pool

def main(pdbfile, centerfile, brickfolder, parallel, n):
    log_file = open('log.txt', 'w')
    sys.stdout = log_file
    pdbname = pdbfile.split('.')[0]
    
    start = time.time()
    # Prepare protein: add hydrogens, add charges and save it as protein.mol2 (with all H), protein_gbsa.pdb(with polar H) and protein_noH.pdb(without H)
    Pre_protein(pdbfile)
    
    # Prepare brick folder: add hydrogens, add charges and save them as mol2
    if str(parallel) != 'True':
       filter_brickfolder = filter_brick_folder(brickfolder)
       new_brickfolder = Pre_brick_folder(filter_brickfolder)
    else:
       filter_brickfolder = filter_brick_folder(brickfolder)
       new_brickfolder = Pre_brick_folder_parallel(filter_brickfolder, n)
    # new_brickfolder = 'new_' + brickfolder
    
    # Generate energy grids and GBSA grids for protein
    mol2_file = pdbname + '.mol2'
    modify_mol2(mol2_file)
    Protein(pdbfile, centerfile)
    
    # Merge brick files into whole file
    os.chdir(new_brickfolder)
    os.system('cat *.mol2 > gen_mols.mol2')
    shutil.copy('gen_mols.mol2', '../')
    os.chdir('../')
    
    # Mapping the whole file onto the receptor
    os.system('python Mapping_eval.py ' + pdbname + ' 1 gen_mols.mol2')
    print('Mapping is finished!')
    
    end2 = time.time()
    print('Running time for evaluation: %s Seconds' % (end2-start))
    log_file.close()

if __name__=="__main__":
    # Usage: python ifitdock_eval.py -i target.pdb -c center.txt -f brickfolder -p 20 --parallel
    # Usage: python ifitdock_eval.py -i target.pdb -c center.txt -f brickfolder
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-i', type=str, help='input pdb file')
    parser.add_argument('-c', type=str, help='pocket center file')
    parser.add_argument('-f', type=str, help='brick folder')
    parser.add_argument('--parallel', default='False', action='store_true', help='prepare bricks in a parallel mode or not')
    parser.add_argument('-p', type=int, help='the number of cpu used in preparing the bricks')
    args = parser.parse_args()
    
    pdbfile = args.i
    centerfile = args.c
    brickfolder = args.f
    parallel = args.parallel
    n = args.p

    main(pdbfile, centerfile, brickfolder, parallel, n)
