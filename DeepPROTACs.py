import os
import sys

AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")

import math
import glob
import json
import shutil
import argparse
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import Descriptors
from multiprocessing import Pool
from frag_utils import get_linker, remove_dummys


# obtain the surrounding environment of a ligand
def obtain_surround(ligand, protein):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()

    cmd.load(ligand)
    cmd.load(protein)
    cmd.create('whole', 'all')
    protein_name = os.path.splitext(protein)[0]
    cmd.delete(os.path.splitext(ligand)[0])
    cmd.delete(protein_name)
    cmd.select('ligand', 'organic')
    cmd.select('surrounding', 'byres ligand around 5')
    surround_mol = os.path.splitext(ligand)[0] + '_surround.mol2'
    cmd.save(surround_mol, 'surrounding')
    cmd.delete('all')

    return surround_mol

# sort rosetta docking results
def sort_docking(top_docking=30):
    with open('local.fasc') as fr:
        lines = fr.readlines()
    lines = lines[2:]
    docking_dict = {}
    for line in lines:
        docking_dict[line.split()[-1]] = float(line.split()[1])
    docking_list = sorted(docking_dict.items(), key = lambda x: x[1])
    results_list = []
    id_list = []
    count = 0
    for docking in docking_list:
        if count < top_docking:
            if docking[0].split('_')[0].split('.')[1] not in id_list:
                id_list.append(docking[0].split('_')[0].split('.')[1])
                results_list.append(docking)
                count += 1
            else:
                continue
        else:
            break
    return results_list

# select PROTACs based on DeepPROTACs model
def deep_PROTACs(results_list, subfolder):
    final_PROTACs = subfolder + '_PROTACs.sdf'
    writer = Chem.SDWriter(final_PROTACs)
    for result in results_list:
        folder = result[0]
        if os.path.isfile(folder + '/final/aligned_mols.sdf') and os.path.getsize(folder + '/final/aligned_mols.sdf') != 0:
            os.system('cp ' + AutoFBDD_FOL + '/DeepPROTACs/* ./')
            print(subfolder + ': ' + folder)
            generated_protacs = Chem.SDMolSupplier(folder + '/final/aligned_mols.sdf')
            clean_frag = Chem.SDMolSupplier(folder + '/' + folder + '_frag.sdf')[0]
            with open(folder + '/molecules_' + folder + '_out.json', 'r') as fr:
                data = json.load(fr)
            smi_in = data[0]['smiles_in']
            for i in range(len(generated_protacs)):
                gen_mol = generated_protacs[i]
                linker = get_linker(gen_mol, clean_frag, smi_in)
                linker = remove_dummys(linker)
                if linker != '':
                    with open('linker_' + str(i) + '.smi', 'w') as fw:
                        fw.write(linker)
                    if not os.path.exists(folder + '/final/' + str(i)):
                        os.mkdir(folder + '/final/' + str(i))
                    #os.chdir(folder)
                    #obtain_surround('warhead.mol2', 'target.pdb')
                    #obtain_surround('E3_ligand.mol2', 'E3_ligase.pdb')
                    #os.chdir('../')
                    os.system('cp ' + folder + '/E3_ligand.mol2 ' + folder + '/final/' + str(i) + '/ligase_ligand.mol2')
                    os.system('cp ' + folder + '/E3_ligand_surround.mol2 ' + folder + '/final/' + str(i) + '/ligase_pocket.mol2')
                    os.system('cp ' + folder + '/warhead.mol2 ' + folder + '/final/' + str(i) + '/target_ligand.mol2')
                    os.system('cp ' + folder + '/warhead_surround.mol2 ' + folder + '/final/' + str(i) + '/target_pocket.mol2')
                    os.system('mv linker_' + str(i) + '.smi ' + folder + '/final/' + str(i) + '/linker.smi')

                    os.system('python deepprotacs.py -p ' + folder + '/final/' + str(i))
                    with open('result.txt', 'r') as fr:
                        lines = fr.readlines()
                    percent = float(lines[0])
                    if percent > 50:
                        writer.write(gen_mol)
                        print(subfolder + ', ' + folder + ': ' + Chem.MolToSmiles(gen_mol) + ' has passed DeepPROTACs evaluation')
    return final_PROTACs

def main(folder, target, E3_ligase, warhead, top_docking):
    DeepPROTACs_log = 'DeepPROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0] + '_log.txt'
    log = open(folder + '/' + DeepPROTACs_log, 'a', buffering=1)
    os.chdir(folder)
    PROTACs_folder = 'PROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0]
    warhead = os.path.basename(warhead)
    subfolder = os.path.splitext(warhead)[0] + '_' + os.path.splitext(E3_ligase)[0]
    if os.path.exists(PROTACs_folder + '/' + subfolder + '/Patchdock_Results'):
        os.chdir(PROTACs_folder + '/' + subfolder + '/Patchdock_Results')
        results_list = sort_docking(top_docking)
        os.chdir('../')
        final_PROTACs = deep_PROTACs(results_list, subfolder)
        log.write('INFO: ' + final_PROTACs + ' has been generated.\n')
        os.chdir(AutoFBDD_FOL + '/' + folder)
    os.chdir(AutoFBDD_FOL)
    
def main2(folder, target, E3_ligase, top_docking):
    DeepPROTACs_log = 'DeepPROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0] + '_log.txt'
    log = open(folder + '/' + DeepPROTACs_log, 'w', buffering=1)
    log.write('INFO: DeepPROTACs has started.\n')
    os.chdir(folder)
    warheads = glob.glob(AutoFBDD_FOL + '/' + folder + '/final_results/*.mol2')
    PROTACs_folder = 'PROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0]
    for warhead in warheads:
        warhead = os.path.basename(warhead)
        subfolder = os.path.splitext(warhead)[0] + '_' + os.path.splitext(E3_ligase)[0]
        if os.path.exists(PROTACs_folder + '/' + subfolder + '/Patchdock_Results'):
            os.chdir(PROTACs_folder + '/' + subfolder + '/Patchdock_Results')
            results_list = sort_docking(top_docking)
            os.chdir('../')
            final_PROTACs = deep_PROTACs(results_list, subfolder)
            os.system('cp ' + final_PROTACs + ' ' + AutoFBDD_FOL + '/' + folder + '/final_results')
            log.write('INFO: ' + final_PROTACs + ' has been generated.\n')
            os.chdir(AutoFBDD_FOL + '/' + folder)
    os.chdir(AutoFBDD_FOL)
    log.write('INFO: DeepPROTACs has finished.\n')
    log.close()


if __name__ == '__main__':
    # Get the arguments
    # Usage: python DeepPROTACs.py --folder folder --target target.pdb --E3_ligase E3_ligase.pdb --top_docking 10 (--parallel --num_cpu 10)
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('--folder', type=str, help='main folder')
    parser.add_argument('--target', type=str, help='target protein')
    parser.add_argument('--E3_ligase', type=str, help='E3 ligase')
    parser.add_argument('--top_docking', type=int, default=30, help='the number of top rosetta docking results extracted')
    parser.add_argument('--parallel', default='False', action='store_true', help='parallel mode or not')
    parser.add_argument('--num_cpu', type=int, default=20, help='the number of cpu used in parallel run')
    args = parser.parse_args()

    folder = args.folder
    target = args.target
    E3_ligase = args.E3_ligase
    top_docking = args.top_docking
    parallel = args.parallel
    n = args.num_cpu
    
    if parallel:
        p = Pool(n)
        DeepPROTACs_log = 'DeepPROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0] + '_log.txt'
        log = open(folder + '/' + DeepPROTACs_log, 'a', buffering=1)
        log.write('INFO: DeepPROTACs has started.\n')
        warheads = glob.glob(AutoFBDD_FOL + '/' + folder + '/final_results/*.mol2')
        for warhead in warheads:
            p.apply_async(main, args=(folder, target, E3_ligase, warhead, top_docking, ))
        print('waiting for all processes')
        p.close()
        p.join()
        PROTACs_folder = 'PROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0]
        os.system('cp ' + AutoFBDD_FOL + '/' + folder +  '/' + PROTACs_folder + '/*/*PROTACs.sdf ' + AutoFBDD_FOL + '/' + folder + '/final_results')
        log.write('INFO: DeepPROTACs has finished.\n')
        log.close()
    else:
        main2(folder, target, E3_ligase, top_docking)

