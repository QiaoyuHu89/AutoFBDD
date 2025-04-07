import sys
import os
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")
import time
import argparse
from ifitdock_func import *
from multiprocessing import Pool


def mapping(pdbname, sep_brickfolder, sep_bricklist, centername):
    os.chdir(sep_brickfolder)
    os.system('python Mapping.py ' + pdbname + ' 1 ' + sep_bricklist + ' ' + centername)
    os.chdir('../')

def main(pdbfile, centerfile, brickfolder, bricklist, parallel, n, sep_num):
    log_file = open('log.txt', 'w')
    sys.stdout = log_file
    pdbname = pdbfile.split('.')[0]
    centername = centerfile.split('.')[0]
    
    start = time.time()
    # Prepare protein: add hydrogens, add charges and save it as protein.mol2 (with all H), protein_gbsa.pdb(with polar H) and protein_noH.pdb(without H)
    Pre_protein(pdbfile)
    
    # Prepare brick folder: add hydrogens, add charges and save them as mol2
    if brickfolder in ['zinc', 'chembl', 'pdbbind', 'all', 'all_0.8', 'all_0.6', 'all_0.4', 'all_adme', 'all_admet']:
        new_brickfolder = AutoFBDD_FOL + '/brick_library/' + brickfolder
    else:
        if str(parallel) != 'True':
            filter_brickfolder = filter_brick_folder(brickfolder)
            new_brickfolder = Pre_brick_folder(filter_brickfolder)
        else:
            filter_brickfolder = filter_brick_folder(brickfolder)
            new_brickfolder = Pre_brick_folder_parallel(filter_brickfolder, n)
    
    # Generate energy grids and GBSA grids for protein
    mol2_file = pdbname + '.mol2'
    modify_mol2(mol2_file)
    Protein(pdbfile, centerfile)
    
    # Separate brick file list into small brick file lists and create brick folders
    bricks = os.listdir(new_brickfolder)
    with open(bricklist, 'a') as fw:
        for brick in bricks:
            fw.write(brick + '\n')
            
    if len(bricks)%sep_num == 0:
        num_brick_list = len(bricks)//sep_num
        for i in range(num_brick_list):
            with open(os.path.splitext(bricklist)[0] + '_%s.txt' % (str(i)), 'a') as fw:
                for j in range(i*sep_num, (i+1)*sep_num):
                    fw.write(bricks[j] + '\n')
    else:
        num_brick_list = len(bricks)//sep_num + 1
        for i in range(num_brick_list):
            with open(os.path.splitext(bricklist)[0] + '_%s.txt' % (str(i)), 'a') as fw:
                if (i+1)*sep_num < len(bricks):
                    for j in range(i*sep_num, (i+1)*sep_num):
                        fw.write(bricks[j] + '\n')
                else:
                    for j in range(i*sep_num, len(bricks)):
                        fw.write(bricks[j] + '\n')

    p = Pool(n)
    for i in range(num_brick_list):
        sep_bricklist = os.path.splitext(bricklist)[0] + '_%s.txt' % (str(i))
        sep_brickfolder = os.path.splitext(sep_bricklist)[0]
        os.mkdir(sep_brickfolder)
        os.mkdir(sep_brickfolder + '/file')
        os.system('cp * ' + sep_brickfolder)
        with open(sep_bricklist, 'r') as fr:
            brick_lines = fr.readlines()
        for brick_line in brick_lines:
            os.system('cp ' + new_brickfolder + '/' + brick_line.strip() + ' ' + sep_brickfolder + '/file')
        
        p.apply_async(mapping, args=(pdbname, sep_brickfolder, sep_bricklist, centername, ))
    print('waiting for all processes')
    p.close()
    p.join()
    print('Mapping is finished!')
    
    # Copy mapping output files from separate folder to current folder and clustering 
    for i in range(num_brick_list):
        sep_bricklist = os.path.splitext(bricklist)[0] + '_%s.txt' % (str(i))
        sep_brickfolder = os.path.splitext(sep_bricklist)[0]
        os.system('cp ' + sep_brickfolder + '/*.mol2.mol2 .')
        os.system('cp ' + sep_brickfolder + '/*scored.list .')
    Clustering(pdbfile, bricklist)
    
    print('Clustering is finished!')
    
    os.system('rm *.mol2.mol2')
    os.system('rm *scored.list')
    # os.system('rm -r ' + bricklist + '_*')
    
    end2 = time.time()
    print('Running time for docking and clustering: %s Seconds' % (end2-start))
    log_file.close()

if __name__=="__main__":
    # Usage: python ifitdock.py -i target.pdb -c center.txt -f brickfolder -l bricks_file.txt -n 20 -s 5 --parallel
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-i', type=str, help='input pdb file')
    parser.add_argument('-c', type=str, help='pocket center file')
    parser.add_argument('-f', type=str, help='brick folder')
    parser.add_argument('-l', type=str, help='brick file list')
    parser.add_argument('--parallel', default='False', action='store_true', help='prepare bricks in a parallel mode or not')
    parser.add_argument('-n', type=int, help='the number of cpu used in parallel run')
    parser.add_argument('-s', type=int, help='the number of bricks in each brick list file')
    args = parser.parse_args()
    
    pdbfile = args.i
    centerfile = args.c
    brickfolder = args.f
    bricklist = args.l
    parallel = args.parallel
    n = args.n
    sep_num = args.s

    main(pdbfile, centerfile, brickfolder, bricklist, parallel, n, sep_num)
