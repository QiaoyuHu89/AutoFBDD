from pymol import cmd
import sys
import os
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")
import re
import glob
import time
import shutil
from rdkit import Chem
from multiprocessing import Pool
from collections import Counter
from utils import dataset_info


# Delete hydrogens of brick sdf file and save as sdf file
def delete_H(brick):
    mol_sdf = Chem.SDMolSupplier(brick)[0]
    brick_noH_sdf = os.path.splitext(brick)[0] + '_noH.sdf'
    if mol_sdf is not None:
        mol_sdf = Chem.RemoveHs(mol_sdf)
        w = Chem.SDWriter(brick_noH_sdf)
        w.write(mol_sdf)
        w.close()
        print('The hydrogens of ' + brick + ' is deleted.')
        return brick_noH_sdf
    else:
        print('The hydrogens of ' + brick + ' can not be deleted.')
        return None

def Pre_protein(pdbfile):
    pdbname = pdbfile.split('.')[0]
    
    with open('protein_template.py', 'r') as fr:
        lines = fr.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace('target', pdbname)
    with open('ifitpre_protein.py', 'w') as fw:
        fw.writelines(lines)
    os.system('$chimera/chimera --nogui --script ./ifitpre_protein.py')
    os.system('rm ifitpre_protein.py')
    print('Preparing protein done!')

def filter_brick_folder(brickfolder):
    brickfiles = os.listdir(brickfolder)
    filter_brickfolder = 'filter_' + brickfolder
    if not os.path.exists(filter_brickfolder):
        os.mkdir(filter_brickfolder)
    for brickfile in brickfiles:
        shutil.copy(brickfolder + '/' + brickfile, filter_brickfolder + '/' + brickfile)
        brick_noH_sdf = delete_H(filter_brickfolder + '/' + brickfile)
        if brick_noH_sdf is not None:
            mol = Chem.SDMolSupplier(brick_noH_sdf)[0]
            for atom in mol.GetAtoms():
                if "%s%i(%i)" % (atom.GetSymbol(), atom.GetTotalValence(), atom.GetFormalCharge()) not in dataset_info('zinc')['atom_types']:
                    os.system('rm ' + filter_brickfolder + '/' + brickfile)
                    print(filter_brickfolder + '/' + brickfile + ' is deleted.')
                    break
                elif atom.GetSymbol() == 'S':
                    neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
                    counter = Counter(neighbor_list)
                    if counter['O'] >= 3:
                        os.system('rm ' + filter_brickfolder + '/' + brickfile)
                        print(filter_brickfolder + '/' + brickfile + ' is deleted.')
                        break
        else:
            os.system('rm ' + filter_brickfolder + '/' + brickfile)
            print(filter_brickfolder + '/' + brickfile + ' is deleted.')
    os.system('rm ' +  filter_brickfolder + '/*_noH.sdf')
    return filter_brickfolder

def Pre_brick_folder(brickfolder):
    brickfiles = os.listdir(brickfolder)
    new_brickfolder = 'new_' + brickfolder
    if not os.path.exists(new_brickfolder):
        os.mkdir(new_brickfolder)
    for brickfile in brickfiles:
        shutil.copy(brickfolder + '/' + brickfile, './')
        brickname = os.path.splitext(brickfile)[0]
        os.system('babel ' + brickfile + ' ' + brickname + '.mol2 -h')
        with open('ligand_template.py', 'r') as fr:
            lines = fr.readlines()
        for i in range(len(lines)):
            lines[i] = lines[i].replace('ligand', brickname)
        with open('ifitpre_ligand.py', 'w') as fw:
            fw.writelines(lines)
        try:
            os.system('$chimera/chimera --nogui --script ./ifitpre_ligand.py')
            os.system('mv ' + brickname + '.mol2 ' + new_brickfolder)
        except:
            pass

        os.system('rm ifitpre_ligand.py')
        os.system('rm ' + brickfile)
    print('Preparing brick folder done!')
    return new_brickfolder

def Pre_brick_folder_new(brickfolder):
    brickfiles = os.listdir(brickfolder)
    new_brickfolder = 'new_' + brickfolder
    if not os.path.exists(new_brickfolder):
        os.mkdir(new_brickfolder)
    for brickfile in brickfiles:
        brickname = os.path.splitext(brickfile)[0]
        os.system('babel ' + brickfolder + '/' + brickfile + ' ' + brickfolder + '/' + brickname + '.mol2 -h')
        os.system('mv ' + brickfolder + '/' + brickname + '.mol2 ' + new_brickfolder)
    print('Preparing brick folder done!')
    return new_brickfolder

def chimera(folder, file):
    os.chdir(folder)
    os.system('$chimera/chimera --nogui --script ./' + file)
    os.chdir('../')

def Pre_brick_folder_parallel(brickfolder, n):
    brickfiles = os.listdir(brickfolder)
    new_brickfolder = 'new_' + brickfolder
    if not os.path.exists(new_brickfolder):
        os.mkdir(new_brickfolder)
    
    start = time.time()
    p = Pool(n)
    for i in range(len(brickfiles)):
        try:
            ifitpre_file = 'ifitpre_ligand_%s.py' % (str(i))
            brickname = os.path.splitext(brickfiles[i])[0]
            if not os.path.exists(brickname):
                os.mkdir(brickname)
            shutil.copy(brickfolder + '/' + brickfiles[i], brickname)
            os.system('babel ' + brickname + '/' + brickfiles[i] + ' ' + brickname + '/' + brickname + '.mol2 -h')
            with open('ligand_template.py', 'r') as fr:
                lines = fr.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].replace('ligand', brickname)
            with open(ifitpre_file, 'w') as fw:
                fw.writelines(lines)
            shutil.move(ifitpre_file, brickname)
            p.apply_async(chimera, args=(brickname, ifitpre_file, ))
        except:
            continue
    print('waiting for all processes')
    p.close()
    p.join()
    end = time.time()
    
    for brickfile in brickfiles:
        brickname = os.path.splitext(brickfile)[0]
        os.system('mv ' + brickname + '/' + brickname + '.mol2 ' + new_brickfolder)
        os.system('rm -r ' + brickname)
    print("Total time of preparing all bricks {} s".format((end - start)))
    print('=' * 100)
    print('Preparing brick folder done!')
    return new_brickfolder

def openbabel(folder, file):
    name = os.path.splitext(file)[0]
    new_folder = 'new_' + folder
    os.system('babel ' + folder + '/' + file + ' ' + folder + '/' + name + '.mol2 -h')
    os.system('mv ' + folder + '/' + name + '.mol2 ' + new_folder)

def Pre_brick_folder_parallel_new(brickfolder, n):
    brickfiles = os.listdir(brickfolder)
    new_brickfolder = 'new_' + brickfolder
    if not os.path.exists(new_brickfolder):
        os.mkdir(new_brickfolder)
    p = Pool(n)
    for brickfile in brickfiles:
        p.apply_async(openbabel, args=(brickfolder, brickfile, ))
    print('waiting for all processes')
    p.close()
    p.join()
    print('Preparing brick folder done!')
    return new_brickfolder

def modify_mol2(mol2_file):
    with open(mol2_file, 'r') as fr:
        lines = fr.readlines()
    start_i = lines.index("@<TRIPOS>BOND\n")
    end_i= lines.index("@<TRIPOS>SUBSTRUCTURE\n")
    for i in range(start_i+1, end_i):
        lines[i] = lines[i].rstrip() + "    DICT\n"
    with open(mol2_file, 'w') as fw:
        fw.writelines(lines)

def is_in_interval(x, minx, maxx):
    """ x: [i,j,k] in (minx, maxx]
    """
    if x[0]> minx[0] and x[0]<=maxx[0]:
        if x[1]> minx[1] and x[1]<=maxx[1]:
            if x[2]> minx[2] and x[2]<=maxx[2]:
                return True
    return False

def Protein(pdbfile, centerfile):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()

    # get the pdb and process
    cmd.load(pdbfile)
    pdbname = pdbfile.split('.')[0]

    # boundary and center
    boundary=cmd.get_extent()
    boundary[0]=[int(x-3) for x in boundary[0]]
    boundary[1]=[int(x+3) for x in boundary[1]]
    d0 = (boundary[1][0] - boundary[0][0])//4
    d1 = (boundary[1][1] - boundary[0][1])//4
    d2 = (boundary[1][2] - boundary[0][2])//4

    # write into the file
    with open('lattice.in','w') as fw:
        fw.write("BOX_DIMENSIONS_LOWER %.0f %.0f %.0f\n" % (boundary[0][0], boundary[0][1],boundary[0][2]))
        fw.write("BOX_DIMENSIONS_HIGHER %.0f %.0f %.0f\n" % (boundary[1][0], boundary[1][1],boundary[1][2]))
        fw.write("GRID_DIMENSIONS %d %d %d\n" % (d0, d1, d2))

    # get box grid, outfile = outlattice_center.pdb
    os.system('./lattice 1 lattice.in out')

    # get the grid file "modify_grid.pdb"
    cmd.load('outlattice_center.pdb')
    cmd.select('outer__',pdbname+' around 3')
    cmd.select('inner__',pdbname+' around 1')
    cmd.remove('outlattice_center and not outer__')
    cmd.remove('outlattice_center and inner__')
    cmd.save('modify_grid.pdb',"outlattice_center")
    cmd.delete('all')

    with open(centerfile, 'r') as fr:
        center_coor = fr.readline().split()
        center_i, center_j, center_k = float(center_coor[0]), float(center_coor[1]), float(center_coor[2])
    centername = centerfile.split('.')[0]
    center_file = centername + '_1.txt'

    ## get coor
    with open('modify_grid.pdb') as fr:
        modify_grid = fr.readlines()
    modify_grid_1 = [[float(x[30:38]),float(x[38:46]),float(x[46:54])] for x in modify_grid[:-1]]

    center_file_content = []
    ## if grid_i in each small grid
    for grid_i in range(len(modify_grid_1)):
        x0, y0, z0 = modify_grid_1[grid_i]
        if abs(x0-center_i)<5 and abs(y0-center_j)<5 and abs(z0-center_k)<5:
            center_file_content.append(modify_grid_1[grid_i])

    ## assign modify_grid, delete less than 30
    with open(center_file, 'w') as fw:
        for line in center_file_content:
            fw.write("%f %f %f\n"%(line[0],line[1],line[2]))

    # Energy grid generation;
    os.system("perl mapping_nrg_grid_generator.pl %s 1 %s" % (pdbname, centerfile))

    # GB Grid Generation;("GB Grid Generation") mapping_GB_grid_generator.pl
    if os.path.exists("INCHEM_GB"):
        os.system("mv INCHEM_GB INCHEM")
        os.system("perl mapping_GB_grid_generator.pl %s 1 %s" % (pdbname, centerfile))
        os.system("mv INCHEM INCHEM_GB")

    # SA Grid Generation;("SA Grid Generation") mapping_SA_grid_generator.pl
    if os.path.exists("INCHEM_SA"):
        os.system("mv INCHEM_SA INCHEM")
        os.system("perl mapping_SA_grid_generator.pl %s 1 %s" % (pdbname, centerfile))
        os.system("mv INCHEM INCHEM_SA")

def Mapping(pdbfile, bricklist, centerfile):
    # Docking and Mapping;("Docking and Mapping") Mapping.pl
    pdbname = pdbfile.split('.')[0]
    centername = centerfile.split('.')[0]
    os.system("python Mapping.py %s 1 %s %s" % (pdbname, bricklist, centername))

def Clustering(pdbfile, bricklist):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()
        
    # Cluster the binding modes of the chemical probes and output result.("Clustering") Cluster.pl
    pdbname = pdbfile.split('.')[0]
    with open(bricklist, 'r') as fr:
        lines = fr.readlines()
    if len(lines) <= 10:
        os.system("perl cluster.pl %s 1 %s  1 3 1 1 cluster_out" % (pdbname, bricklist))
    elif len(lines) <= 100:
        os.system("perl cluster.pl %s 1 %s  1 3 2 1 cluster_out" % (pdbname, bricklist))
    elif len(lines) <= 1000:
        os.system("perl cluster.pl %s 1 %s  2 3 2 1 cluster_out" % (pdbname, bricklist))
    elif len(lines) <= 10000:
        os.system("perl cluster.pl %s 1 %s  3 3 2 1 cluster_out" % (pdbname, bricklist))
    elif len(lines) <= 50000:
        os.system("perl cluster.pl %s 1 %s  3 3 10 1 cluster_out" % (pdbname, bricklist))
    elif len(lines) > 50000:
        os.system("perl cluster.pl %s 1 %s  3 2 1000 0.5 cluster_out" % (pdbname, bricklist))

    # # deal with result
    # for file_name in os.listdir('.'):
    #     if re.match(r'cluster_out_\d+',file_name):
    #         cmd.load(file_name)
    #         cmd.split_states(file_name.split('.')[0], prefix="result")
    #         cmd.select("result*")
    #         result_i = file_name.split('.')[0].split('_')[2]
    #         result_name = "result_"+ result_i +".pdb"
    #         cmd.save(result_name)
    #         cmd.delete("all")
    #         cmd.load(result_name)
    #         cmd.load(pdbfile)
    #         cmd.color("green", pdbname)
    #         cmd.select('result_res',result_name.split(".")[0]+' around 5')
    #         cmd.color("red", "result_res")
    #         cmd.save("result_res_"+ result_i + ".pdb","result_res")
    #         cmd.save("result_"+ result_i +".pse")
    #         cmd.delete("all")

    # final_filename = pdbname + "_result"
    # result_num = glob.glob("result*")
    # if len(result_num) >0:
    #     os.system("mkdir "+ final_filename)
    #     os.system("mv result* " + final_filename)
    #     os.system("zip -r "+final_filename+".zip "+final_filename)
    #     with open('result.txt','w') as f:
    #         f.write('success')
    # else:
    #     with open('result.txt','w') as f:
    #         f.write("error")