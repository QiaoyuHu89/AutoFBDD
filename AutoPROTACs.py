import os
import sys
from tokenize import Name

AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
PatchDock = os.environ["PATCHDOCK"]
FOLDER = AutoFBDD_FOL + "/Rosetta/"
ROSETTA_FOL = os.environ["ROSETTA_FOL"]
SCRIPTS = ROSETTA_FOL + "/main/source/bin/rosetta_scripts.default.linuxgccrelease"

sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")

import math
import glob
import json
import shutil
import argparse
from copy import deepcopy
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import Descriptors
from multiprocessing import Pool
from rdkit.Chem import rdMolAlign
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
from data.prepare_data_linker_design import read_file_edit, preprocess
from data_pre_for_linking_pharms import obtain_dummy, cal_dummy_dist, cal_dis_ang, obtain_linker_pharma
from brick_linking_pharms import linking_bricks
from assess_linking_pharms import delete_du, assess_mols_RMSD, assess_mols_pharm, extract_sdf


# Delete hydrogens of brick mol2 file and save as sdf file
def delete_H(brick):
    mol_mol2 = Chem.MolFromMol2File(brick)
    brick_pdb = os.path.splitext(brick)[0] + '.pdb'
    brick_sdf = os.path.splitext(brick)[0] + '.sdf'
    brick_new_sdf = os.path.splitext(brick)[0] + '_new.sdf'
    os.system('babel ' + brick + ' ' + brick_pdb)
    os.system('babel ' + brick_pdb + ' ' + brick_sdf)
    os.system('babel ' + brick + ' ' + brick_new_sdf)
    
    mol_sdf = Chem.SDMolSupplier(brick_sdf)[0]
    mol_new_sdf = Chem.SDMolSupplier(brick_new_sdf)[0]
    brick_noH_sdf = os.path.splitext(brick)[0] + '_noH.sdf'
    if mol_mol2 is not None:
        mol_mol2 = Chem.RemoveHs(mol_mol2)
        w = Chem.SDWriter(brick_noH_sdf)
        w.write(mol_mol2)
        w.close()
        print('The hydrogens of ' + brick + ' is deleted.')
        os.system('rm ' + brick_pdb + ' ' + brick_sdf + ' ' + brick_new_sdf)
        return brick_noH_sdf
    elif mol_sdf is not None:
        mol_sdf = Chem.RemoveHs(mol_sdf)
        w = Chem.SDWriter(brick_noH_sdf)
        w.write(mol_sdf)
        w.close()
        print('The hydrogens of ' + brick + ' is deleted.')
        os.system('rm ' + brick_pdb + ' ' + brick_sdf + ' ' + brick_new_sdf)
        return brick_noH_sdf
    elif mol_new_sdf is not None:
        mol_new_sdf = Chem.RemoveHs(mol_new_sdf)
        w = Chem.SDWriter(brick_noH_sdf)
        w.write(mol_new_sdf)
        w.close()
        print('The hydrogens of ' + brick + ' is deleted.')
        os.system('rm ' + brick_pdb + ' ' + brick_sdf + ' ' + brick_new_sdf)
        return brick_noH_sdf
    else:
        print('The hydrogens of ' + brick + ' can not be deleted.')
        os.system('rm ' + brick_pdb)
        os.system('rm ' + brick_sdf)
        return None
    
# Add hydrogens of ligand mol2 file and save as sdf file
def add_H(ligand):
    #tmp_pdb = ligand.split('.')[0] + '.pdb'
    #os.system('babel ' + ligand + ' ' + tmp_pdb + ' -h')
    ligand_sdf = os.path.splitext(ligand)[0] + '_H.sdf'
    #os.system('babel ' + tmp_pdb + ' ' + ligand_sdf)
    #os.remove(tmp_pdb)
    os.system('babel ' + ligand + ' ' + ligand_sdf + ' -h')
    return ligand_sdf

# align generated molecules to original fragments
def align_mols(best_RMSD_mols, frag_sdf, output_sdf):
    gen_sdfs = Chem.SDMolSupplier(best_RMSD_mols)
    ref_mol = Chem.SDMolSupplier(frag_sdf)[0]
    writer = Chem.SDWriter(output_sdf)
    gen_smis = []
    for i in range(len(gen_sdfs)):
        gen_mol = gen_sdfs[int(i)]
        score = sascorer.calculateScore(gen_mol)
        gen_smi = Chem.MolToSmiles(gen_mol)
        if gen_mol is not None and score <= 6.0 and gen_smi not in gen_smis:
            gen_smis.append(gen_smi)
            pyO3A = rdMolAlign.GetO3A(gen_mol, ref_mol).Align()
            try:
                writer.write(gen_mol)
            except:
                continue
            print('write one aligned molecule')
        else:
            pass
    writer.close()
    print('aligned_RMSD_mols.sdf is created successfully!')

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

# obtain the anchor based on the sum of distance between the ligand and surrounding res
def obtain_anchor(ligand, surround_mol):
    ligand_sdf = delete_H(ligand)
    if ligand_sdf is not None:
        ligand_mol = Chem.SDMolSupplier(ligand_sdf)[0]
        # ligand_mol = Chem.MolFromMol2File(ligand, removeHs=False, sanitize=False)
        ligand_conf = ligand_mol.GetConformer()
        surround = Chem.MolFromMol2File(surround_mol)
        surround_conf = surround.GetConformer()
        ligand_pos_list = []
        ligand_id_list = []
        atom_pos_list = []
        for atom in ligand_mol.GetAtoms():
            ligand_pos_list.append(list(ligand_conf.GetAtomPosition(atom.GetIdx())))
            ligand_id_list.append(atom.GetIdx())
        for atom in surround.GetAtoms():
            atom_pos_list.append(list(surround_conf.GetAtomPosition(atom.GetIdx())))

        dist_dict = {}
        for i in ligand_pos_list:
            sum_dist = 0
            for j in atom_pos_list:
                dist = 0
                for k in [0, 1, 2]:
                    dist += (i[k] - j[k])**2
                dist = round(math.sqrt(dist), 2)
                sum_dist += dist
            index = ligand_pos_list.index(i)
            dist_dict[ligand_id_list[index]] = sum_dist
        dist_list = sorted(dist_dict.items(), key = lambda x: x[1], reverse=True)
        print(dist_list)
        
        for dist in dist_list:
            atom = ligand_mol.GetAtomWithIdx(dist[0])
            # if 'H' in [x.GetSymbol() for x in atom.GetNeighbors()]:
            if atom.GetSymbol() == 'C' and atom.GetExplicitValence() != 4:
                anchor_id = dist[0]
                break
        print(anchor_id)
        return anchor_id
    else:
        raise IOError('Hydrogens can not be deleted.')

def translate_anchors(old_mol2, new_sdf, old_anchor):
    OldMol = Chem.MolFromMol2File(old_mol2, sanitize=True)
    NewSdf = Chem.SDMolSupplier(new_sdf, sanitize=True)[0]
    print(Chem.MolToSmiles(OldMol))
    print(Chem.MolToSmiles(NewSdf))
    NewMatch = NewSdf.GetSubstructMatch(OldMol)
    if len(NewMatch) == 0:
        return -1
    print(NewMatch)
    return NewMatch[old_anchor]

# run PatchDock
def patchdock(structs, anchors, min_dist, max_dist, num_results=1000, threshold=4.0):
    [structA, structB] = structs
    [anchorA, anchorB] = anchors
    with open('Patchdock_cst', 'w') as f:
        f.write(' '.join([str(s) for s in [anchorA, anchorB, min_dist, max_dist]]) + '\n')
    os.system(PatchDock + '/buildParamsFine.pl ' + structA + ' ' + structB + ' ' + '2 ' + 'EI' + '; mv params.txt Patchdock_params.txt')
    os.system('sed -i \'s/#distanceConstraintsFile file_name/distanceConstraintsFile Patchdock_cst/g\' Patchdock_params.txt')
    os.system('sed -i \'s/clusterParams 0.1 4 2.0 4.0/clusterParams 0.1 4 2.0 ' + str(threshold) + '/g\' Patchdock_params.txt')
    os.system(PatchDock + '/patch_dock.Linux Patchdock_params.txt Patchdock_output')
    os.system(PatchDock + '/transOutput.pl Patchdock_output 1 ' + str(num_results))

    try:
        os.mkdir('Patchdock_Results')
    except FileExistsError:
        print("Patchdock_Results already exists, delete and recreate the directory.")
        shutil.rmtree("Patchdock_Results")
        os.mkdir('Patchdock_Results')

    results = [int(line.split('.')[1]) for line in glob.glob('Patchdock_output.*')]
    if len(results) == 0:
        return None
    Num_Results = max(results)
    for i in range(Num_Results):
        os.system('mv Patchdock_output.' + str(i + 1) + '.pdb Patchdock_Results/pd.' + str(i + 1) + '.pdb')
        os.system('sed -i \'s/E3l X/E3l Y/g\' Patchdock_Results/pd.' + str(i + 1) + '.pdb')
    return Num_Results

# run the mol_to_param Rosetta script
def mol_to_params(ligand, name, pdb, overwrite = True, conformers = False, nbr = -1):
    line = "python " + ROSETTA_FOL + "/main/source/scripts/python/public/molfile_to_params.py " + ligand + " -n " + name + " -p " + pdb
    if overwrite:
        line += " --clobber"
    if conformers:
        line += " --conformers-in-one-file"
    if not nbr == -1:
        line += " --nbr_atom=" + str(nbr)
    os.system(line)
    if conformers:
        return pdb + '.pdb', pdb + '.params'
    return pdb + '_0001.pdb', pdb + '.params'

# run a Rosetta local docking protocol
def local_docking(struct, chainsA, chainsB, ptA_params, ptB_params, nstruct = 50):
    return SCRIPTS + " -s " + struct + " -parser:protocol " + FOLDER + "docking.xml @" + FOLDER + "docking.flags @" + FOLDER + \
    "relax.flags -extra_res_fa " + ptA_params + " -extra_res_fa " + ptB_params + " -partners " + chainsA + "_" + chainsB + \
    " -scorefile local.fasc -nstruct " + str(nstruct) + " -overwrite"

# Rosetta local docking
def run_local_docking(i, chains, PT_params, Local):
    os.system(local_docking('pd.' + str(i + 1) + '.pdb', chains[0] + 'X', chains[1] + 'Y', PT_params[0], PT_params[1], Local))
    
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

# obtain two anchors id and distance from rosetta docking results
def docking_postprocess(ternary_complex, warhead, E3_ligand, warhead_anchor_id, ligand_anchor_id):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()
        
    cmd.load(ternary_complex)
    cmd.select('warhead', 'chain X')
    cmd.select('E3_ligand', 'chain Y')
    cmd.select('target', 'chain A')
    cmd.select('E3_ligase', 'chain B')
    cmd.select('ternary_protein', 'polymer')
    #warhead_pdb = 'warhead.pdb'
    #E3_ligand_pdb = 'E3_ligand.pdb'
    new_warhead = 'warhead.mol2'
    new_E3_ligand = 'E3_ligand.mol2'
    target = 'target.pdb'
    E3_ligase = 'E3_ligase.pdb'
    target_E3 = 'target_E3.pdb'
    #cmd.save(warhead_pdb, 'warhead')
    #cmd.save(E3_ligand_pdb, 'E3_ligand')
    cmd.save(new_warhead, 'warhead')
    cmd.save(new_E3_ligand, 'E3_ligand')
    cmd.save(target, 'target')
    cmd.save(E3_ligase, 'E3_ligase')
    cmd.save(target_E3, 'ternary_protein')
    cmd.delete('all')
    
    #new_warhead = 'warhead.mol2'
    #new_E3_ligand = 'E3_ligand.mol2'
    #os.system('babel ' + warhead_pdb + ' ' + new_warhead)
    #os.system('babel ' + E3_ligand_pdb + ' ' + new_E3_ligand)
    warhead_surround = obtain_surround(new_warhead, target)
    # warhead_anchor_id = obtain_anchor(warhead, warhead_surround)
    ligand_surround = obtain_surround(new_E3_ligand, E3_ligase)
    # ligand_anchor_id = obtain_anchor(E3_ligand, ligand_surround)
    
    os.system('babel ' + warhead + ' ' + warhead)
    os.system('babel ' + E3_ligand + ' ' + E3_ligand)
    os.system('babel ' + new_warhead + ' ' + new_warhead)
    os.system('babel ' + new_E3_ligand + ' ' + new_E3_ligand)
    mol1 = Chem.MolFromMol2File(warhead)
    mol2 = Chem.MolFromMol2File(E3_ligand)
    new_mol1 = Chem.MolFromMol2File(new_warhead)
    new_mol2 = Chem.MolFromMol2File(new_E3_ligand)
    if mol1 is not None and new_mol1 is not None:
        pyO3A_1 = rdMolAlign.GetO3A(mol1, new_mol1).Align()
        dist_list1 = []
        conf1 = mol1.GetConformer()
        new_conf1 = new_mol1.GetConformer()
        for atom in new_mol1.GetAtoms():
            warhead_anchor_coor = list(conf1.GetAtomPosition(warhead_anchor_id))
            atom_coor = list(new_conf1.GetAtomPosition(atom.GetIdx()))
            dist = 0
            for k in [0, 1, 2]:
                dist += (warhead_anchor_coor[k] - atom_coor[k])**2
            dist = round(math.sqrt(dist), 2)
            dist_list1.append(dist)
        new_warhead_anchor_id = dist_list1.index(min(dist_list1))
    if mol2 is not None and new_mol2 is not None:
        pyO3A_2 = rdMolAlign.GetO3A(mol2, new_mol2).Align()
        dist_list2 = []
        conf2 = mol2.GetConformer()
        new_conf2 = new_mol2.GetConformer()
        for atom in new_mol2.GetAtoms():
            ligand_anchor_coor = list(conf2.GetAtomPosition(ligand_anchor_id))
            atom_coor = list(new_conf2.GetAtomPosition(atom.GetIdx()))
            dist = 0
            for k in [0, 1, 2]:
                dist += (ligand_anchor_coor[k] - atom_coor[k])**2
            dist = round(math.sqrt(dist), 2)
            dist_list2.append(dist)
        new_ligand_anchor_id = dist_list2.index(min(dist_list2))
    if new_mol1 is not None and new_mol2 is not None:
        new_conf1 = new_mol1.GetConformer()
        new_conf2 = new_mol2.GetConformer()
        new_warhead_anchor_coor = list(new_conf1.GetAtomPosition(new_warhead_anchor_id))
        new_ligand_anchor_coor = list(new_conf2.GetAtomPosition(new_ligand_anchor_id))
        print(new_warhead_anchor_coor, new_ligand_anchor_coor)
        dist = 0
        for k in [0, 1, 2]:
            dist += (new_warhead_anchor_coor[k] - new_ligand_anchor_coor[k])**2
        dist = round(math.sqrt(dist), 2)
        delete_H(new_warhead)
        new_mol1 = Chem.SDMolSupplier(os.path.splitext(new_warhead)[0] + '_noH.sdf')[0]
        new_mol1_atom_num = new_mol1.GetNumAtoms()
        new_ligand_anchor_id = new_ligand_anchor_id + new_mol1_atom_num
    return new_warhead, new_E3_ligand, target_E3, new_warhead_anchor_id, new_ligand_anchor_id, dist

# linker generation
def generate_linker(results_list, warhead, E3_ligand, warhead_anchor_id, ligand_anchor_id):
    for result in results_list:
        folder = result[0]
        ternary_complex = folder + '.pdb'
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.system('cp Patchdock_Results/' + ternary_complex + ' ' + folder)
        os.system('cp ' + warhead + ' ' + E3_ligand + ' ' + folder)
        os.chdir(folder)
        folder_log = open(folder + '_log.txt', 'a')
        new_warhead, new_E3_ligand, target_E3, new_warhead_anchor_id, new_ligand_anchor_id, dist = docking_postprocess(ternary_complex, warhead, E3_ligand, warhead_anchor_id, ligand_anchor_id)
        du_sdf, du_surround_mol = obtain_dummy(new_warhead, new_E3_ligand, new_warhead_anchor_id, new_ligand_anchor_id, dist, target_E3, folder)
        min_dist = cal_dummy_dist(du_surround_mol)
        print(folder + ' min_dist: ' + str(min_dist))
        if min_dist >= 2.0:
            try:
                os.system('cp ' + AutoFBDD_FOL +'/DEVELOP/models/linker_design/pretrained_DEVELOP_model_pharms.pickle .')
                data_path, frag_sdf = cal_dis_ang(new_warhead, new_E3_ligand, new_warhead_anchor_id, new_ligand_anchor_id, folder)
                linker_pharma_sdf = folder + '_linker_pharma.sdf'
                pharm_count = obtain_linker_pharma(du_surround_mol, frag_sdf, linker_pharma_sdf)
                types = folder + '.types'
                with open(types, 'w') as fw:
                    fw.write('1 ' + linker_pharma_sdf + ' ' + frag_sdf)
                
                raw_data = read_file_edit(data_path, pharm_count, add_idx=True, calc_pharm_counts=True)
                preprocess(raw_data, "zinc", os.path.splitext(data_path)[0], './', False)
                json_name = 'molecules_' + os.path.splitext(data_path)[0] + '.json'
                
                os.system('rm ' + folder + '_sc*')
                os.system('rm ' + folder + '_bb*')
                folder_log.write(folder + ' is preprocessed.\n')
                
                with open(json_name, 'r') as fr:
                    lines = fr.readlines()
                if lines[0] != str([]):
                    if round(float(dist)) < 2:
                        min = str(0)
                        max = str(round(float(dist))+2)
                    else:
                        min = str(round(float(dist))-2)
                        max = str(round(float(dist))+2)
                    if not os.path.exists('generated_smi_1'):
                        os.mkdir('generated_smi_1')
                    linking_bricks(json_name, types, folder, min, max, str(1))
                folder_log.close()
            except Exception as e:
                os.system('rm ' + folder + '_sc*')
                os.system('rm ' + folder + '_bb*')
                folder_log.write(e + '\n')
                folder_log.write(folder + ' is discarded due to errors in data preparation process of linking.\n')
                print(folder + ' is discarded due to errors in data preparation process of linking.')
                folder_log.close()

        else:
            folder_log.write(folder + ' is discarded due to the block between two anchors.\n')
            print(folder + ' is discarded due to the block between two anchors.')
            folder_log.close()

        os.chdir('../')

def keep_chirality(generated_smi):
    frags = []
    ori_smis = []
    gen_smis = []
    gen_mols = []

    with open(generated_smi, 'r') as fr:
        smis_list = fr.readlines()
        for i in smis_list:
            frags.append(i.strip().split(' ')[0])
            ori_smis.append(i.strip().split(' ')[1])
            gen_smis.append(i.strip().split(' ')[2])

    for i in gen_smis:
        gen_mol = Chem.MolFromSmiles(i)
        gen_mols.append(gen_mol)

    def spam(n):
        out=[]
        for perm in getPerms(n):
            elem = [ int(i) for i in list(perm) ]
            out.append(elem)
        return out

    def getPerms(n):
        from itertools import permutations
        for i in getCandidates(n):
            for perm in set(permutations(i)):
                yield ''.join(perm)

    def getCandidates(n):
        for i in range(0, n+1):
            res = "1" * i + "0" * (n - i)
            yield res

    def GetStereoIsomers(mol):
        out = []
        smis = []
        outmol = deepcopy(mol)
        chiralCentres = Chem.FindMolChiralCenters(outmol, includeUnassigned=True)

        #return the molecule object when no chiral centres where identified
        if chiralCentres == []:
            return [mol]
        #All bit permutations with number of bits equals number of chiralCentres
        elif len(chiralCentres) >= 1 and len(chiralCentres) <= 10:
            elements = spam(len(chiralCentres))
            for isoId,element in enumerate(elements):       
                for centreId,i in enumerate(element):
                    atomId = chiralCentres[centreId][0]
                    if i == 0:
                        outmol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
                    elif i == 1:
                        outmol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
                out.append(outmol)
                smis.append(Chem.MolToSmiles(outmol,isomericSmiles=True))
            return out, smis
        else:
            return []

    ligand1_smi = frags[0].split('.')[0].replace('*', '')
    ligand2_smi = frags[0].split('.')[1].replace('*', '')
    ligand1 = Chem.MolFromSmiles(ligand1_smi)
    ligand2 = Chem.MolFromSmiles(ligand2_smi)
    gen_chiral_smis = []
    for mol in gen_mols:
        if len(GetStereoIsomers(mol)) > 1:
            for i in GetStereoIsomers(mol)[1]:
                gen_mol = Chem.MolFromSmiles(i)
                ligand1_list = gen_mol.GetSubstructMatch(ligand1,useChirality=True)
                ligand2_list = gen_mol.GetSubstructMatch(ligand2,useChirality=True)
                if len(ligand1_list) != 0 and len(ligand2_list) != 0:
                    gen_chiral_smis.append(Chem.MolToSmiles(gen_mol))
        elif len(GetStereoIsomers(mol)) == 1 and Chem.FindMolChiralCenters(ligand1, includeUnassigned=True) == [] \
        and Chem.FindMolChiralCenters(ligand2, includeUnassigned=True) == []:
            gen_chiral_smis.append(Chem.MolToSmiles(mol))

    with open(generated_smi, 'w') as fw:
        for gen_chiral_smi in gen_chiral_smis:
            line = frags[0] + ' ' + ori_smis[0] + ' ' + gen_chiral_smi + '\n'
            fw.write(line)
    return generated_smi

# evaluate generated PROTACs
def evaluate_PROTACs(results_list, n):
    for result in results_list:
        folder = result[0]
        if os.path.exists(folder + '/generated_smi_1'):
            os.chdir(folder)
            folder_log = open(folder + '_log.txt', 'a')
            os.system('cp ' + AutoFBDD_FOL + '/DEVELOP/analysis/evaluate_linking_mols*.py .')
            os.system('cp generated_smi_1/*.smi .')
            gen_smi = folder + '_generated_smiles.smi'
            gen_smi = keep_chirality(gen_smi)
            if os.path.getsize(gen_smi) == 0:
               folder_log.write('Can not generate mols with the same chirality with warhead and E3_ligand, skip this conformation!\n')
               os.chdir('../')
               continue
            else:
                folder_log.write('Generate mols with the same chirality with warhead and E3_ligand.\n')
                frag_sdf = folder + '_frag.sdf'
                delete_du(frag_sdf)
                try:
                    assess_mols_RMSD('ZINC', gen_smi, frag_sdf, AutoFBDD_FOL + '/DEVELOP/data/linker_design/data_zinc_train.txt', \
                        './', folder, str(n), 'True', 'None', AutoFBDD_FOL + '/DEVELOP/analysis/wehi_pains.csv', 'True')
                    os.system('mv ' + folder + ' ' + folder + '_RMSD')
                    log_file = folder + '_RMSD_log.txt'
                    extract_sdf(log_file, 'gen_RMSD_sdfs', 'gen_RMSD_smis', folder + '_RMSD')
                    best_RMSD_mols = 'gen_RMSD_sdfs/best_mols_RMSD_Frag.sdf'
                    output_sdf = 'aligned_RMSD_mols.sdf'
                    align_mols(best_RMSD_mols, frag_sdf, output_sdf)
                    folder_log.write(folder + ' has generated valid RMSD sdf.\n')
                    print(folder + ' has generated valid RMSD sdf.')
                except:
                    folder_log.write(folder + ' does not generate valid RMSD sdf.\n')
                    print(folder + ' does not generate valid RMSD sdf.')
                    folder_log.close()
                    os.chdir('../')
                    continue
                
                try:
                    assess_mols_pharm('ZINC', gen_smi, frag_sdf, AutoFBDD_FOL + '/DEVELOP/data/linker_design/data_zinc_train.txt', \
                        './', folder, str(n), 'True', 'None', AutoFBDD_FOL + '/DEVELOP/analysis/wehi_pains.csv', '.', 'True')
                    os.system('mv '  + folder + ' ' + folder + '_pharm')
                    log_file = folder + '_pharm_log.txt'
                    extract_sdf(log_file, 'gen_pharm_sdfs', 'gen_pharm_smis', folder + '_pharm')
                    best_RMSD_mols = 'gen_pharm_sdfs/best_mols_RMSD_Frag.sdf'
                    output_sdf = 'aligned_pharm_mols.sdf'
                    align_mols(best_RMSD_mols, frag_sdf, output_sdf)
                    folder_log.write(folder + ' has generated valid pharm sdf.\n')
                    print(folder + ' has generated valid pharm sdf.')
                except:
                    folder_log.write(folder + ' does not generate valid pharm sdf.\n')
                    print(folder + ' does not generate valid pharm sdf.')
                    folder_log.close()
                    os.chdir('../')
                    continue
                
                folder_log.write(folder + ' has been evaluated!\n')
                print(folder + ' has been evaluated!')
                folder_log.close()
                os.chdir('../')
        else:
            continue

# select PROTACs based on distance, energy filtering
def merge_sdf(results_list):
    for result in results_list:
        folder = result[0]
        if os.path.exists(folder + '/generated_smi_1'):
            os.chdir(folder)
            if not os.path.exists('final'):
                os.mkdir('final')
            if os.path.exists('aligned_RMSD_mols.sdf') and os.path.getsize('aligned_RMSD_mols.sdf') !=0 \
            and os.path.exists('aligned_pharm_mols.sdf') and os.path.getsize('aligned_pharm_mols.sdf') != 0:
                writer_1 = Chem.SDWriter('final/aligned_mols.sdf')
                aligned_RMSD_sdf = 'aligned_RMSD_mols.sdf'
                aligned_RMSD_mols = Chem.SDMolSupplier(aligned_RMSD_sdf)
                aligned_pharm_sdf = 'aligned_pharm_mols.sdf'
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
            elif os.path.exists('aligned_RMSD_mols.sdf') and os.path.getsize('aligned_RMSD_mols.sdf') !=0 \
            and not os.path.exists('aligned_pharm_mols.sdf'):
                os.system('cp aligned_RMSD_mols.sdf final/aligned_mols.sdf')
            elif os.path.exists('aligned_RMSD_mols.sdf') and os.path.getsize('aligned_RMSD_mols.sdf') !=0 \
            and os.path.exists('aligned_pharm_mols.sdf') and os.path.getsize('aligned_pharm_mols.sdf') == 0:
                os.system('cp aligned_RMSD_mols.sdf final/aligned_mols.sdf')
            os.chdir('../')

def main(folder, target, warhead_anchor_id, E3_ligase, E3_ligand, ligand_anchor_id, min_dist, max_dist, n, top_docking):
    AutoPROTACs_log = 'AutoPROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0] + '_log.txt'
    log = open(folder + '/' + AutoPROTACs_log, 'w', buffering=1)
    log.write('INFO: AutoPROTACs has started.\n')
    os.chdir(folder)
    PROTACs_folder = 'PROTACs_' + target.split('.')[0] + '_' + E3_ligase.split('.')[0]
    if not os.path.exists(PROTACs_folder):
        os.mkdir(PROTACs_folder)
    warheads = glob.glob(AutoFBDD_FOL + '/' + folder + '/final_results/*.mol2')
    for warhead in warheads[0:30]:
        try:
            subfolder = os.path.splitext(os.path.basename(warhead))[0] + '_' + os.path.splitext(E3_ligase)[0]
            if not os.path.exists(PROTACs_folder + '/' + subfolder):
                os.mkdir(PROTACs_folder + '/' + subfolder)
            os.system('cp ' + target + ' ' + E3_ligase + ' ' + E3_ligand + ' ' + PROTACs_folder + '/' + subfolder)
            os.system('cp ' + warhead + ' ' + PROTACs_folder + '/' + subfolder)
            warhead = os.path.basename(warhead)
            os.chdir(PROTACs_folder + '/' + subfolder)

            chains = ['A', 'B']
            warhead_surround = obtain_surround(warhead, target)
            ligand_surround = obtain_surround(E3_ligand, E3_ligase)
            if warhead_anchor_id == -1:
                warhead_anchor_id = obtain_anchor(warhead, warhead_surround)
            else:
                warhead_anchor_id = warhead_anchor_id - 1
            if ligand_anchor_id == -1:
                ligand_anchor_id = obtain_anchor(E3_ligand, ligand_surround)
            else:
                ligand_anchor_id = ligand_anchor_id - 1
            anchors = [warhead_anchor_id, ligand_anchor_id]
            warhead_H = add_H(warhead)
            E3_ligand_H = add_H(E3_ligand)
            # warhead_anchor_id = translate_anchors(warhead, warhead_H, warhead_anchor_id)
            # ligand_anchor_id = translate_anchors(E3_ligand, E3_ligand_H, ligand_anchor_id)
            # anchors = [warhead_anchor_id, ligand_anchor_id]
            warhead_pdb, warhead_param = mol_to_params(warhead_H, 'War', 'War')
            E3_ligand_pdb, E3_ligand_param = mol_to_params(E3_ligand_H, 'E3l', 'E3l')
            PT_params = [warhead_param, E3_ligand_param]
            os.system('cat ' + warhead_pdb + ' ' + target + ' > target_warhead.pdb')
            os.system('cat ' + E3_ligand_pdb + ' ' + E3_ligase + ' > E3_ligase_ligand.pdb')
            Structs = ['target_warhead.pdb', 'E3_ligase_ligand.pdb']

            # Patchdock
            log.write('INFO: Running PatchDock with the constrains for ' + subfolder + '.\n')
            Num_Results = patchdock(Structs, [a + 1 for a in anchors], min_dist, max_dist, 1000, 2.0)
            if Num_Results == None:
               log.write('INFO: PatchDock did not find any global docking solution within the geometrical constraints for ' + subfolder + '.\n')
               continue

            # Rosetta Local Docking
            log.write('INFO: Run Rosetta local docking on the top n PatchDock results for ' + subfolder + '.\n')
            os.system('cp ' + warhead_param + ' ' + E3_ligand_param + ' Patchdock_Results/')
            os.chdir('Patchdock_Results/')
            p = Pool(n)
            for i in range(Num_Results):
                p.apply_async(run_local_docking, args=(i, chains, PT_params, 10, ))
            print('waiting for all processes')
            p.close()
            p.join()

            # sort rosetta docking results
            log.write('INFO: Sort rosetta docking results.\n')
            results_list = sort_docking(top_docking)

            # linker generation
            log.write('INFO: Generating linker between warhead and E3_ligand.\n')
            os.chdir('../')
            generate_linker(results_list, warhead, E3_ligand, warhead_anchor_id, ligand_anchor_id)

            # evaluate PROTACs
            log.write('INFO: Evaluate generated PROTACs.\n')
            evaluate_PROTACs(results_list, n)

            # merge generated PROTACs
            log.write('INFO: Merge generated PROTACs.\n')
            merge_sdf(results_list)

            os.chdir(AutoFBDD_FOL + '/' + folder)
        except Exception as e:
            print(str(e))
            os.chdir(AutoFBDD_FOL + '/' + folder)
            continue
    os.chdir(AutoFBDD_FOL)
    log.write('INFO: AutoPROTACs has finished.\n')
    log.close()


if __name__=="__main__":
    # Get the arguments
    # Usage (select anchor automatically): python AutoPROTACs.py --folder folder --target target.pdb --E3_ligase E3_ligase.pdb --E3_ligand E3_ligand.mol2 \
    # --min_dist 5 --max_dist 15 --num_cpu 20 --top_docking 100
    # Usage (select anchor manually) python AutoPROTACs.py --folder folder --target target.pdb --warhead_anchor_id num --E3_ligase E3_ligase.pdb \
    # --E3_ligand E3_ligand.mol2 --ligand_anchor_id num --min_dist 5 --max_dist 15 --num_cpu 20 --top_docking 100
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('--folder', type=str, help='main folder')
    parser.add_argument('--target', type=str, help='target protein')
    parser.add_argument('--warhead_anchor_id', type=int, default=-1, help='warhead anchor ID')
    parser.add_argument('--E3_ligase', type=str, help='E3 ligase')
    parser.add_argument('--E3_ligand', type=str, help='E3 ligand')
    parser.add_argument('--ligand_anchor_id', type=int, default=-1, help='ligand anchor ID')
    parser.add_argument('--min_dist', type=int, default=5, help='the minimum length of linker')
    parser.add_argument('--max_dist', type=int, default=15, help='the maximum length of linker')
    parser.add_argument('--num_cpu', type=int, default=20, help='the number of cpu used in parallel run')
    parser.add_argument('--top_docking', type=int, default=30, help='the number of top rosetta docking results extracted')
    args = parser.parse_args()

    folder = args.folder
    target = args.target
    warhead_anchor_id = args.warhead_anchor_id
    E3_ligase = args.E3_ligase
    E3_ligand = args.E3_ligand
    ligand_anchor_id = args.ligand_anchor_id
    min_dist = args.min_dist
    max_dist = args.max_dist
    n = args.num_cpu
    top_docking = args.top_docking

    main(folder, target, warhead_anchor_id, E3_ligase, E3_ligand, ligand_anchor_id, min_dist, max_dist, n, top_docking)
