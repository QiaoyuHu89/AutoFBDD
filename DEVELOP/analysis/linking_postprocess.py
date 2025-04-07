import os
import re
from socket import if_indextoname
import sys
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")

import glob
import math
import time
import shutil
import argparse
from numpy import mean
from rdkit import Chem
from rdkit.Chem import Descriptors
from pymol import cmd
from multiprocessing import Pool
from collections import defaultdict
from data_pre_for_linking import delete_H


# Calculate the distance between generated mol and surrounding res
def cal_sur_dist(gen_mol, dummy_surround_mol):
    gen_conf = gen_mol.GetConformer()
    res_conf = dummy_surround_mol.GetConformer()
    mol_pos_list = []
    res_pos_list = []
    for atom in gen_mol.GetAtoms():
        mol_pos_list.append(list(gen_conf.GetAtomPosition(atom.GetIdx())))
    for atom in dummy_surround_mol.GetAtoms():
        if atom.GetSymbol() != '*':
            res_pos_list.append(list(res_conf.GetAtomPosition(atom.GetIdx())))

    dist_list = []
    for i in mol_pos_list:
        for j in res_pos_list:
            dist = 0
            for k in [0, 1, 2]:
                dist += (i[k] - j[k])**2
            dist = round(math.sqrt(dist), 2)
            dist_list.append(dist)
    dist_list.sort()
    return dist_list[0]

# Add hydrogens using pymol
def addH(file):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quite flag
        pymol2.SingletonPyMOL().start()

    cmd.load(file)
    cmd.h_add("all")
    cmd.save(file)
    cmd.delete("all")
    
# Add hydrogens and convert sdf to mol2 using chimera
def addH_convert(file):
    os.system('python chimera_addH.py -i ' + file)
    
# Add hydrogens, charge and convert sdf to mol2 using chimera
def addC_convert(file):
    os.system('python chimera_addC.py -i ' + file)

# Calculate energy difference between lowest conformation and current conformation of a molecule
def cal_energy_diff(sdf_file):
    # addH_convert(sdf_file)
    # mol2_file = sdf_file.replace('sdf', 'mol2')
    # pdb_file = sdf_file.replace('sdf', 'pdb')
    # os.system('babel ' + pdb_file + ' ' + mol2_file)
    # os.system('rm ' + pdb_file)
    mol2_file = sdf_file.replace('sdf', 'mol2')
    os.system('babel ' + sdf_file + ' ' + mol2_file + ' -h')
    mol2_result = mol2_file.split('.')[0] + '_result.mol2'
    energy_diff_list = []
    n, m = 0, 0
    while n < 10:
        os.system('./Cyndi -i ' + mol2_file + ' -o ' + mol2_result)
        energy_list = []
        with open(mol2_result, 'r') as fr:
            lines = fr.readlines()
        if lines != []:
            n += 1
        else:
            m += 1
            if m < 20:
                continue
            else:
                break
        for line in lines:
            if line.startswith('# Energy:'):
                energy_list.append(line.split()[2])
        energy_diff_list.append(float(energy_list[0].replace('kJ/mol', '')))
    if len(energy_diff_list) >= 5:
        ave_energy_diff = mean(energy_diff_list)
    else:
        ave_energy_diff = None
    return mol2_file, mol2_result, ave_energy_diff

def split_mol2(file, n):
    with open(file, 'r') as fr:
        lines = fr.readlines()
    idx = []
    mol2 = []
    dict = defaultdict(list)
    for k,va in [(v,i) for i,v in enumerate(lines)]:
        dict[k].append(va)
    for line in lines:
        if line.startswith('# Name:') and line not in mol2:
            mol2.append(line)
            if len(dict[line]) == 1:
                idx.append(dict[line][0])
            else:
                for j in range(len(dict[line])):
                    idx.append(dict[line][j])
    idx.sort()

    if len(idx) > n:
        for i in range(n):
            with open(file.split('.')[0] + '_' + str(i+1) + '.mol2', 'w') as fw:
                for line in lines[idx[i]:idx[i+1]]:
                    fw.write(line)
    else:
        for i in range(len(idx)-1):
            with open(file.split('.')[0] + '_' + str(i+1) + '.mol2', 'w') as fw:
                for line in lines[idx[i]:idx[i+1]]:
                    fw.write(line)
        with open(file.split('.')[0] + '_' + str(len(idx)) + '.mol2', 'w') as fw:
                for line in lines[idx[-1]:]:
                    fw.write(line)

def main(name, pdbfile, centerfile, eval_folder, num_pose, iter_idx):
    log_file = open(os.path.basename(name) + '_final_log.txt', 'a')
    sys.stdout = log_file
    print('Dealing with ' + os.path.basename(name))
    if not os.path.exists(name + '/final'):
        os.mkdir(name + '/final')

    dummy_surround_mol = Chem.MolFromMol2File('dummy_surround_mol_' + str(iter_idx) + '/' + os.path.basename(name) + 
                                              '_dummy_surround.mol2')

    if os.path.exists(name + '/RMSD/aligned_RMSD_mols.sdf') and os.path.getsize(name + '/RMSD/aligned_RMSD_mols.sdf') !=0 \
    and os.path.exists(name + '/pharm/aligned_pharm_mols.sdf') and os.path.getsize(name + '/pharm/aligned_pharm_mols.sdf') != 0:
        writer_1 = Chem.SDWriter(name + '/final/aligned_mols.sdf')
        aligned_RMSD_sdf = name + '/RMSD/aligned_RMSD_mols.sdf'
        aligned_RMSD_mols = Chem.SDMolSupplier(aligned_RMSD_sdf)
        aligned_pharm_sdf = name + '/pharm/aligned_pharm_mols.sdf'
        aligned_pharm_mols = Chem.SDMolSupplier(aligned_pharm_sdf)
        aligned_smi_list = []
        for i in range(len(aligned_RMSD_mols)):
            writer_1.write(aligned_RMSD_mols[i])
            aligned_smi_list.append(Chem.MolToSmiles(aligned_RMSD_mols[i]))
        for j in range(len(aligned_pharm_mols)):
            if Chem.MolToSmiles(aligned_pharm_mols[j]) not in aligned_smi_list:
                aligned_smi_list.append(Chem.MolToSmiles(aligned_pharm_mols[j]))
                writer_1.write(aligned_pharm_mols[j])
        writer_1.close()
    elif os.path.exists(name + '/RMSD/aligned_RMSD_mols.sdf') and os.path.getsize(name + '/RMSD/aligned_RMSD_mols.sdf') !=0 \
    and not os.path.exists(name + '/pharm/aligned_pharm_mols.sdf'):
        os.system('cp ' + name + '/RMSD/aligned_RMSD_mols.sdf ' + name + '/final/aligned_mols.sdf')
    elif os.path.exists(name + '/RMSD/aligned_RMSD_mols.sdf') and os.path.getsize(name + '/RMSD/aligned_RMSD_mols.sdf') !=0 \
    and os.path.exists(name + '/pharm/aligned_pharm_mols.sdf') and os.path.getsize(name + '/pharm/aligned_pharm_mols.sdf') == 0:
        os.system('cp ' + name + '/RMSD/aligned_RMSD_mols.sdf ' + name + '/final/aligned_mols.sdf')

    # Perform distance filtering
    if os.path.exists(name + '/final/aligned_mols.sdf') and os.path.getsize(name + '/final/aligned_mols.sdf') != 0:
        writer_2 = Chem.SDWriter(name + '/final/distance_filter_mols.sdf')
        aligned_sdf = name + '/final/aligned_mols.sdf'
        aligned_mols = Chem.SDMolSupplier(aligned_sdf)
        print('Number of generated mols original: ' + str(len(aligned_mols)))
        print('#'*100)
        for i in range(len(aligned_mols)):
            min_dist = cal_sur_dist(aligned_mols[i], dummy_surround_mol)
            if min_dist >= 1.5:
                writer_2.write(aligned_mols[i])
                print(Chem.MolToSmiles(aligned_mols[i]) + ' is remained, min dist: ' + str(min_dist))
            else:
                print(Chem.MolToSmiles(aligned_mols[i]) + ' is discarded, min dist: ' + str(min_dist))
        writer_2.close()
        print('#'*100)
    elif os.path.exists(name + '/final/aligned_mols.sdf') and os.path.getsize(name + '/final/aligned_mols.sdf') == 0:
        print('Number of generated mols original: 0')
        print('#'*100)

    # Perform energy filtering
    if os.path.exists(name + '/final/distance_filter_mols.sdf') and os.path.getsize(name + '/final/distance_filter_mols.sdf') != 0:
        writer_3 = Chem.SDWriter(name + '/final/energy_filter_mols.sdf')
        distance_sdf = name + '/final/distance_filter_mols.sdf'
        distance_mols = Chem.SDMolSupplier(distance_sdf)
        print('Number of generated mols original: ' + str(len(aligned_mols)))
        print('Number of generated mols after distance filtering: ' + str(len(distance_mols)))
        print('#'*100)
        for i in range(len(distance_mols)):
            sdf_file = name + '/final/' + str(i) + '.sdf'
            writer_4 = Chem.SDWriter(sdf_file)
            writer_4.write(distance_mols[i])
            mol2_file, mol2_result, energy_diff = cal_energy_diff(sdf_file)
            if energy_diff is not None:
                if energy_diff >= -30:
                    writer_3.write(distance_mols[i])
                    print(Chem.MolToSmiles(distance_mols[i]) + ' is remained, energy difference: ' + str(energy_diff))
                else:
                    print(Chem.MolToSmiles(distance_mols[i]) + ' is discarded, energy difference: ' + str(energy_diff))
            else:
                print(Chem.MolToSmiles(distance_mols[i]) + ' is discarded because Cyndi can not generate valid conformation!')
            os.system('rm ' + sdf_file)
            os.system('rm ' + mol2_file)
            os.system('rm ' + mol2_result)
        energy_sdf = name + '/final/energy_filter_mols.sdf'
        energy_mols = Chem.SDMolSupplier(energy_sdf)
        writer_3.close()
        writer_4.close()
        print('#'*100)
    elif os.path.exists(name + '/final/distance_filter_mols.sdf') and os.path.getsize(name + '/final/distance_filter_mols.sdf') == 0:
        print('Number of generated mols original: ' + str(len(aligned_mols)))
        print('Number of generated mols after distance filtering: 0')
        print('#'*100)

    # Perform ifitdock evaluation
    pwd = os.getcwd()
    if os.path.exists(name + '/final/energy_filter_mols.sdf') and os.path.getsize(name + '/final/energy_filter_mols.sdf') != 0:
        energy_sdf = name + '/final/energy_filter_mols.sdf'
        energy_mols = Chem.SDMolSupplier(energy_sdf)
        print('Number of generated mols original: ' + str(len(aligned_mols)))
        print('Number of generated mols after distance filtering: ' + str(len(distance_mols)))
        print('Number of generated mols after energy filtering: ' + str(len(energy_mols)))
        print('#'*100)

        if not os.path.exists(name + '/ifit_eval'):
            os.mkdir(name + '/ifit_eval')
        os.system('cp ./ifit_eval/* ' + name + '/ifit_eval')
        if not os.path.exists(name + '/ifit_eval/' + eval_folder):
            os.mkdir(name + '/ifit_eval/' + eval_folder)
        for i in range(len(energy_mols)):
            w = Chem.SDWriter(name + '/ifit_eval/' + eval_folder + '/' + str(i) + '.sdf')
            w.write(energy_mols[i])
            w.close()
        os.chdir(name + '/ifit_eval')
        os.system('python ifitdock_eval.py -i ' + pdbfile + ' -c ' + centerfile + ' -f ' + eval_folder)
        screening_sorted_file = pdbfile.split('.')[0] + '_screening_sorted.mol2'
        if os.path.exists(screening_sorted_file) and os.path.getsize(screening_sorted_file) != 0:
            split_mol2(screening_sorted_file, num_pose)
            os.system('cp ' + screening_sorted_file.split('.')[0] + '_* ../final')

        if not os.path.exists(pwd + '/brickfolder_' + str(iter_idx+1)):
            os.mkdir(pwd + '/brickfolder_' + str(iter_idx+1))
        screening_mol2 = glob.glob('*_screening_sorted_*.mol2')
        if screening_mol2 != []:
            for mol2 in screening_mol2:
                try:
                    with open(mol2, 'r') as fr:
                        lines = fr.readlines()
                    for line in lines:
                        if line.startswith('delta_G'):
                            delta_G = float(line.strip().split(' = ')[1])
                    sdf_file = delete_H(mol2)
                    mol = Chem.SDMolSupplier(sdf_file)[0]
                    molwt = Descriptors.MolWt(mol)
                    if delta_G <= -10.0 and molwt <= 500:
                        os.system('cp ' + mol2 + ' ' + pwd + '/brickfolder_' + str(iter_idx+1) + '/' + name.split('/')[-1] + 
                                '_screening' + mol2.split('.')[0].split('_')[-1] + '.mol2') 
                except:
                    pass
        print('ifitdock evaluation of ' + name + ' is done!')
        os.chdir(pwd)
    elif os.path.exists(name + '/final/energy_filter_mols.sdf') and os.path.getsize(name + '/final/energy_filter_mols.sdf') == 0:
        print('Number of generated mols original: ' + str(len(aligned_mols)))
        print('Number of generated mols after distance filtering: ' + str(len(distance_mols)))
        print('Number of generated mols after energy filtering: 0')
        print('#'*100)
    print('='*100 + '\n')
    log_file.close()
    os.system('mv ' + os.path.basename(name) + '_final_log.txt' + ' ' + name + '/final')
    
    # # Perform MD simulation and mmpbsa
    # pwd = os.getcwd()
    # if os.path.getsize(name + '/final/energy_filter_mols.sdf') != 0:
    #     binding_energy_dict = {}
    #     writer_5 = Chem.SDWriter(name + '/final/final_mols.sdf')
    #     energy_mols = Chem.SDMolSupplier(name + '/final/energy_filter_mols.sdf')
    #     if not os.path.exists(name + '/MD'):
    #         os.mkdir(name + '/MD')
    #     for i in range(len(energy_mols)):
    #         if not os.path.exists(name + '/MD/' + str(i)):
    #             os.mkdir(name + '/MD/' + str(i))
    #         energy_sdf_file = name + '/MD/' + str(i) + '/' + str(i) + '.sdf'
    #         w = Chem.SDWriter(energy_sdf_file)
    #         w.write(energy_mols[i])
    #         addC_convert(energy_sdf_file)
    #         energy_mol2_file = name + '/MD/' + str(i) + '/' + str(i) + '.mol2'
    #         energy_pdb_file = name + '/MD/' + str(i) + '/' + str(i) + '.pdb'
    #         os.system('babel ' + energy_mol2_file + ' ' + energy_pdb_file)
    #         with open(energy_mol2_file, 'r') as fr:
    #             lines = fr.readlines()
    #         charge = 0
    #         for line in lines[lines.index('@<TRIPOS>ATOM\n')+1:lines.index('@<TRIPOS>BOND\n')]:
    #             charge += float(line.strip().split()[-1])
    #         charge = round(charge)
    #         shutil.copy('Protein.pdb', name + '/MD/' + str(i))
    #         shutil.copy(energy_pdb_file, name + '/MD/' + str(i) + '/Ligand.pdb')
    #         shutil.copy('AutoMD_complex.sh', name + '/MD/' + str(i))
    #         shutil.copy('acpype.py', name + '/MD/' + str(i))
    #         os.chdir(name + '/MD/' + str(i))
    #         os.system('bash AutoMD_complex.sh ' + charge)
            
    #         if os.path.isfile('summary_energy.dat'):
    #             with open('summary_energy.dat','r') as fr:
    #                 lines = fr.readlines()
    #             for line in lines:
    #                 if line.startswith(' Binding energy '):
    #                     binding_energy = float(line.strip().split()[3])
    #             binding_energy_dict[i] = binding_energy
    #         os.chdir(pwd)

    #     binding_energy_ordered_dict = sorted(binding_energy_dict.items(), key=lambda x: x[1])
    #     if len(binding_energy_ordered_dict) >= 10:
    #         for i in range(10):
    #             writer_5.write(energy_mols[binding_energy_ordered_dict[i][0]])
    #     else:
    #         for i in range(len(binding_energy_ordered_dict)):
    #             writer_5.write(energy_mols[binding_energy_ordered_dict[i][0]])


if __name__=="__main__":
    # Usage: python linking_postprocess.py -i target.pdb -c center.txt -e eval_folder -n 20 -p 3 -t iter_idx
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-i', type=str, help='input pdb file')
    parser.add_argument('-c', type=str, help='pocket center file')
    parser.add_argument('-e', type=str, help='evaluation folder')
    parser.add_argument('-n', type=int, help='the number of cpu used in parallel run')
    parser.add_argument('-p', type=int, help='the number of poses saved')
    parser.add_argument('-t', type=int, help='iteration index')
    args = parser.parse_args()
    
    pdbfile = args.i
    centerfile = args.c
    eval_folder = args.e
    num_cpu = args.n
    num_pose = args.p
    iter_idx = args.t
    
    start = time.time()
    p = Pool(num_cpu)
    name_list = glob.glob('eval_results_' + str(iter_idx) + '/success/*')
    name_list.sort()
    for name in name_list:
        p.apply_async(main, args=(name, pdbfile, centerfile, eval_folder, num_pose, iter_idx, ))
    print('waiting for all processes')
    p.close()
    p.join()
    
    end = time.time()
    print('='*100)
    print('Running time: %s Seconds' % (end-start))
    print('Bricks linking postprocess done!')
