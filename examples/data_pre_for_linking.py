import os
import sys
from telnetlib import DO
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")

import math
import time
import copy
import shutil
import argparse
import frag_utils
import numpy as np
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from multiprocessing import Pool
from data.prepare_data_linker_design import read_file_edit, preprocess
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


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
        os.system('rm ' + brick_pdb + ' ' + brick_sdf + ' ' + brick_new_sdf)
        return None

# calculate the coordinates of a series of points on a line given the coordinates of two end points
def cal_points_coor(point1, point2):
    x1, y1, z1 = point1[0], point1[1], point1[2]
    x2, y2, z2 = point2[0], point2[1], point2[2]
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    num_point = int(round(dist, 2) // 1.0)
    increment = (x2 - x1) / (num_point + 1)
    dummy_atoms_coor = []
    for i in range(1, num_point+1):
        x = x1 + increment*i
        y = ((x - x1)/(x2 - x1)) * (y2 - y1) + y1
        z = ((x - x1)/(x2 - x1)) * (z2 - z1) + z1
        dummy_atoms_coor.append([round(x, 4), round(y, 4), round(z, 4)])
    return dummy_atoms_coor

# calculate the coordinates of a specific points on a line given the coordinates of two end points
def cal_point_coor(point1, point2):
    x1, y1, z1 = point1[0], point1[1], point1[2]
    x2, y2, z2 = point2[0], point2[1], point2[2]
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    x = 3.0/dist * (x2 - x1) + x1
    y = 3.0/dist * (y2 - y1) + y1
    z = 3.0/dist * (z2 - z1) + z1
    dummy_atom_coor = [round(x, 4), round(y, 4), round(z, 4)]
    return dummy_atom_coor

# obtain the dummy line between two bricks and its surrounding environment
def obtain_dummy(brick1, brick2, atom_idx_1, atom_idx_2, dist, target, name):
    brick1_sdf = delete_H(brick1)
    brick2_sdf = delete_H(brick2)
    if brick1_sdf is not None and brick2_sdf is not None:
        brick1_mol = Chem.SDMolSupplier(brick1_sdf)[0]
        brick2_mol = Chem.SDMolSupplier(brick2_sdf)[0]
        num_point = int(round(dist, 2) // 1.0)

        combo_no_exit = Chem.CombineMols(brick1_mol, brick2_mol)
        combo = Chem.CombineMols(brick1_mol, brick2_mol)
        for _ in range(num_point):
            combo = Chem.CombineMols(combo, Chem.MolFromSmiles("*"))
        combo_2d = Chem.Mol(combo)
        edcombo = Chem.EditableMol(combo_2d)
        mol_to_link = edcombo.GetMol()
        Chem.SanitizeMol(mol_to_link)
        num_heavy_atoms = mol_to_link.GetNumHeavyAtoms()
        mol_to_link = Chem.AddHs(mol_to_link)
        AllChem.ConstrainedEmbed(mol_to_link, combo_no_exit, randomseed=42)
        mol_to_link = Chem.RemoveHs(mol_to_link)

        conf = mol_to_link.GetConformer()
        ref_conf = combo_no_exit.GetConformer()
        ini_pos = []
        for i in [atom_idx_1, atom_idx_2]:
            pos = list(ref_conf.GetAtomPosition(i))
            ini_pos.append(pos)
        point1, point2 = ini_pos[0], ini_pos[1]
        dummy_pos_list = cal_points_coor(point1, point2)
        dummy_list = []
        for i in range(num_point):
            dummy_num = num_heavy_atoms + i
            dummy_list.append(dummy_num)
        for i in range(len(dummy_list)):
            atom = dummy_list[i]
            pos = dummy_pos_list[i]
            conf.SetAtomPosition(atom, pos)
        conf.SetId(0)
        _ = mol_to_link.AddConformer(conf)

        dummy_sdf = name + '_dummy.sdf'
        w = Chem.SDWriter(dummy_sdf)
        w.write(mol_to_link)
        w.close()

        if cmd._COb is None:
            import pymol2
            import pymol.invocation
            pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
            pymol2.SingletonPyMOL().start()

        cmd.load(dummy_sdf)
        cmd.load(target)
        cmd.create('whole', 'all')
        sdf_name = os.path.splitext(dummy_sdf)[0]
        target_name = os.path.splitext(target)[0]
        cmd.delete(sdf_name)
        cmd.delete(target_name)
        cmd.select('organic_atoms', 'organic and name H+B+C+N+O+S+P+F+Cl+Br+I')
        cmd.select('dummy', 'organic and not organic_atoms')
        cmd.remove('organic_atoms')
        cmd.select('surrounding', 'byres dummy expand 6')
        dummy_surround_mol = sdf_name + '_surround.mol2'
        cmd.save(dummy_surround_mol, 'surrounding')
        cmd.delete('all')

        return dummy_sdf, dummy_surround_mol
    else:
        raise IOError('Hydrogens can not be deleted.')

# calculate the distance list between points given two lists of coordinates
def cal_dis_coord(pos_list_1, pos_list_2):
    dist_list = []
    for i in pos_list_1:
        for j in pos_list_2:
            dist = 0
            for k in [0, 1, 2]:
                dist += (i[k] - j[k])**2
            dist = round(math.sqrt(dist), 2)
            dist_list.append(dist)
    dist_list.sort()
    return dist_list

# calculate the minimum distance between the dummy line and surrounding res given a mol2 file
def cal_dummy_dist(mol2):
    mol = Chem.MolFromMol2File(mol2)
    conf = mol.GetConformer()
    dummy_pos_list = []
    atom_pos_list = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == '*':
            dummy_pos_list.append(list(conf.GetAtomPosition(atom.GetIdx())))
        else:
            atom_pos_list.append(list(conf.GetAtomPosition(atom.GetIdx())))
    dist_list = cal_dis_coord(dummy_pos_list, atom_pos_list)
    return dist_list[0]

# check block between two bricks in a brick folder
def check_block_folder(app_bricks, brick_folder, target, name):
    if not os.path.exists('dummy_sdf'):
        os.mkdir('dummy_sdf')
        
    if not os.path.exists('dummy_surround_mol'):
        os.mkdir('dummy_surround_mol')

    new_app_bricks = os.path.splitext(app_bricks)[0] + '_new.list'
    with open(app_bricks, 'r') as fr:
        lines = fr.readlines()
    for line in lines[1:]:
        newline = line.strip().split(', ')
        brick1, brick2, atom_idx_1, atom_idx_2, dist = brick_folder + '/' + newline[0], brick_folder + '/' + newline[1], int(newline[2]), int(newline[3]), float(newline[4])

        du_sdf, du_surround_mol = obtain_dummy(brick1, brick2, atom_idx_1, atom_idx_2, dist, target, name)
        min_dist = cal_dummy_dist(du_surround_mol)

        with open(new_app_bricks, 'a') as fw:
            if min_dist >= 2.0:
                fw.write(line)
                shutil.move(du_sdf, 'dummy_sdf')
                shutil.move(du_surround_mol, 'dummy_surround_mol')
                print(brick1 + ' and ' + brick2 + ' are appropriate.')
            else:
                os.system('rm ' + du_sdf)
                os.system('rm ' + du_surround_mol)
                print(brick1 + ' and ' + brick2 + ' are not appropriate.')
    print('='*100)
    print('Check block done!')
    return new_app_bricks

# connect two bricks directly without any other atom
def connect_mol(brick1, brick2, atom_idx_1, atom_idx_2, name):
    brick1_sdf = delete_H(brick1)
    brick2_sdf = delete_H(brick2)
    try:
        brick1_mol = Chem.SDMolSupplier(brick1_sdf)[0]
        brick2_mol = Chem.SDMolSupplier(brick2_sdf)[0]
        combo = Chem.CombineMols(brick1_mol, brick2_mol)
        edcombo = Chem.EditableMol(combo)
        edcombo.AddBond(atom_idx_1, atom_idx_2, order=Chem.rdchem.BondType.SINGLE)
        comb_mol = edcombo.GetMol()
        Chem.SanitizeMol(comb_mol)
        AllChem.ConstrainedEmbed(comb_mol, combo, randomseed=42)

        connect_sdf = name + '_connect.sdf'
        connect_mol2 = name + '_connect.mol2'
        w = Chem.SDWriter(connect_sdf)
        w.write(comb_mol)
        w.close()
        os.system('obabel ' + connect_sdf + ' -O' + connect_mol2 + ' -h')
        return connect_mol2, connect_sdf
    except Exception as e:
        print(str(e))
        raise IOError('Can not connect two bricks directly!')

# calculate the distance and angle between two bricks given the anchors
def cal_dis_ang(brick1, brick2, atom_idx_1, atom_idx_2, name):
    brick1_sdf = delete_H(brick1)
    brick2_sdf = delete_H(brick2)
    try:
        brick1_mol = Chem.SDMolSupplier(brick1_sdf)[0]
        brick2_mol = Chem.SDMolSupplier(brick2_sdf)[0]

        combo_no_exit = Chem.CombineMols(brick1_mol, brick2_mol)
        combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))
        edcombo = Chem.EditableMol(combo)
        num_heavy_atoms = combo.GetNumHeavyAtoms()
        edcombo.AddBond(num_heavy_atoms, atom_idx_1, order=Chem.rdchem.BondType.SINGLE)
        edcombo.AddBond(num_heavy_atoms+1, atom_idx_2, order=Chem.rdchem.BondType.SINGLE)
        mol_to_link = edcombo.GetMol()
        Chem.SanitizeMol(mol_to_link)

        # Convert exit vectors to carbons for conformer generation
        du = Chem.MolFromSmiles('*')
        mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link, du, Chem.MolFromSmiles('C'), True)[0]
        Chem.SanitizeMol(mol_to_link_carbon)
        # Generate conformer
        mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
        AllChem.ConstrainedEmbed(mol_to_link_carbon, combo_no_exit, randomseed=42)
        mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)

        # Add this conformer to the two unlinked fragments
        conf = mol_to_link.GetConformer()
        ref_conf = mol_to_link_carbon.GetConformer()
        for i in range(mol_to_link_carbon.GetNumAtoms()):
            pos = list(ref_conf.GetAtomPosition(i))
            conf.SetAtomPosition(i, pos)
        conf.SetId(0)
        _ = mol_to_link.AddConformer(conf)

        frag_sdf = name + '_frag.sdf'

        w = Chem.SDWriter(frag_sdf)
        w.write(mol_to_link)
        w.close()

        # Get distance and angle between fragments
        dist, ang = frag_utils.compute_distance_and_angle(mol_to_link, "", Chem.MolToSmiles(mol_to_link))

        # Write data to file
        data_path = name + '_out.txt'
        with open(data_path, 'w') as f:
            f.write("%s %s %s" % (Chem.MolToSmiles(mol_to_link), dist, ang))
        return data_path, frag_sdf, dist
    except:
        raise IOError('Can not calculate dis and ang of two bricks!')

# calculate the distance and angle between two bricks in a brick folder
def cal_dis_ang_folder(app_bricks, brick_folder, name):
    if not os.path.exists('data_path'):
        os.mkdir('data_path')

    with open(app_bricks, 'r') as fr:
        lines = fr.readlines()
    for line in lines[1:]:
        newline = line.strip().split(', ')
        brick1, brick2, atom_idx_1, atom_idx_2 = brick_folder + '/' + newline[0], brick_folder + '/' + newline[1], int(newline[2]), int(newline[3])
        data_path, _, _ = cal_dis_ang(brick1, brick2, atom_idx_1, atom_idx_2, name)
        raw_data = read_file_edit(data_path, add_idx=True, calc_pharm_counts=True)
        preprocess(raw_data, "zinc", os.path.splitext(data_path)[0], './', False)
        json_name = 'molecules_' + os.path.splitext(data_path)[0] + '.json'
        shutil.move(data_path, 'data_path')
        shutil.move(json_name, 'data_path')
        print(brick1 + ' and ' + brick2 + ' are preprocessed.')
    print('='*100)
    print('Preprocess bricks done!')

# calculate angle among three points
def cal_angle(point1, point2, point3):
    vec12 = np.array([(point1[0]-point2[0]), (point1[1]-point2[1]), (point1[2]-point2[2])])
    vec13 = np.array([(point1[0]-point3[0]), (point1[1]-point3[1]), (point1[2]-point3[2])])
    cos_1 = vec12.dot(vec13)/(np.sqrt(vec12.dot(vec12))*np.sqrt(vec13.dot(vec13)))
    angle_1 = np.arccos(cos_1)*180/np.pi

    vec21 = np.array([(point2[0]-point1[0]), (point2[1]-point1[1]), (point2[2]-point1[2])])
    vec23 = np.array([(point2[0]-point3[0]), (point2[1]-point3[1]), (point2[2]-point3[2])])
    cos_2 = vec21.dot(vec23)/(np.sqrt(vec21.dot(vec21))*np.sqrt(vec23.dot(vec23)))
    angle_2 = np.arccos(cos_2)*180/np.pi
    return round(angle_1, 2), round(angle_2, 2)

# obtain the sidechains and backbone of a given molecule
def obtain_sc_bb(mol2_file):
    with open(mol2_file, 'r') as fr:
        lines = fr.readlines()
    sc_list = []
    for line in lines[(lines.index('@<TRIPOS>SUBSTRUCTURE\n') + 2):]:
        sidechain = line.strip().split()[1]
        sc_list.append(sidechain)
    
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()

    sc_bb_file = []
    cmd.load(mol2_file)
    for item in sc_list:
        sc_file = mol2_file.split('_dummy_')[0] + '_sc_' + item + '.mol2'
        cmd.select(item, 'sidechain and resi ' + item[3:6])
        cmd.save(sc_file, item)
        sc_bb_file.append(sc_file)
    bb_file = mol2_file.split('_dummy_')[0] + '_bb.mol2'
    cmd.select('bb', 'backbone')
    cmd.save(bb_file, 'bb')
    cmd.delete('all')
    sc_bb_file.append(bb_file)
    return sc_bb_file

# calculate the minimum distance and distance list between the dummy line and pharmacophore point
def cal_pharm_dummy_dist(pharm_coords, dummy_pos_list):
    dist_list = []
    for dummy_pos in dummy_pos_list:
        dist = 0
        for k in [0, 1, 2]:
            dist += (dummy_pos[k] - pharm_coords[k])**2
        dist = round(math.sqrt(dist), 2)
        dist_list.append(dist)
    return min(dist_list), dist_list

# obtain side chain pharmacophore coordinates and feature environment
def obtain_sc_pharm_ent(sc_new_feats, sc_conf, dummy_pos_list, anchor_1, anchor_2):
    pharm_coords_list = []
    dist_list = []
    for sc_ent in sc_new_feats:
        tmp_coords = []
        tmp_idxs = []
        for idx in sc_ent[1]:
            tmp_coords.append(list(sc_conf.GetAtomPosition(idx)))
            tmp_idxs.append(idx)
        if len(tmp_idxs) > 0:
            tmp_pharm_coords = list(np.average(tmp_coords, axis=0))
            angle_1, angle_2 = cal_angle(anchor_1, anchor_2, tmp_pharm_coords)
            if angle_1 <= 95 and angle_2 <= 95:
                pharm_coords_list.append(tmp_pharm_coords)
    if len(pharm_coords_list) > 1:
        for tmp_pharm_coords in pharm_coords_list:
            dist, _ = cal_pharm_dummy_dist(tmp_pharm_coords, dummy_pos_list)
            dist_list.append(dist)
        sc_pharm_coords = pharm_coords_list[dist_list.index(min(dist_list))]
        sc_ent = sc_new_feats[dist_list.index(min(dist_list))]
    elif len(pharm_coords_list) == 1:
        sc_pharm_coords = pharm_coords_list[0]
        sc_ent = sc_new_feats[0]
    else:
        sc_pharm_coords = None
        sc_ent = None
    return sc_pharm_coords, sc_ent

# create linker pharmacophore model based on side chain
def create_sc_linker_pharm(dummy_pos_list, sc_pharm_coords, atom_pos_list, antisymbol, pharms, sc_ent, fake_lig, coords, dummy_count):
    min_dist, dist_list = cal_pharm_dummy_dist(sc_pharm_coords, dummy_pos_list)
    dist_dict = {}
    for i,v in enumerate(dist_list):
        dist_dict[i] = v
    dist_dict_order = sorted(dist_dict.items(), key=lambda x:x[1], reverse=False)
    if min_dist <= 6.0:
        for item in dist_dict_order:
            if dummy_count[item[0]] == min(dummy_count):
                dummy_count[item[0]] += 1
                dummy_coords = dummy_pos_list[item[0]]
                linker_pharm_coords = cal_point_coor(sc_pharm_coords, dummy_coords)
                dummy_pos_list_new = cal_points_coor(linker_pharm_coords, dummy_coords)
                if dummy_pos_list_new != []:
                    dummy_min_dist = cal_dis_coord(dummy_pos_list_new, atom_pos_list)[0]
                    pharm_min_dist = cal_dis_coord([linker_pharm_coords], atom_pos_list)[0]
                    if dummy_min_dist > 1.0 and pharm_min_dist > 1.0:
                        fake_lig += '.' + antisymbol[pharms.index(sc_ent[0])]
                        coords.append(linker_pharm_coords)
                else:
                    fake_lig += '.' + antisymbol[pharms.index(sc_ent[0])]
                    coords.append(linker_pharm_coords)
                break
    return coords, fake_lig

# obtain backbone pharmacophore coordinates and feature environment
def obtain_bb_pharm_ent(bb_new_feats, bb_conf):
    pharm_coords_list = []
    bb_ents = []
    for bb_ent in bb_new_feats:
        tmp_coords = []
        tmp_idxs = []
        for idx in bb_ent[1]:
            tmp_coords.append(list(bb_conf.GetAtomPosition(idx)))
            tmp_idxs.append(idx)
        if len(tmp_idxs) > 0:
            tmp_pharm_coords = list(np.average(tmp_coords, axis=0))
            pharm_coords_list.append(tmp_pharm_coords)
            bb_ents.append(bb_ent)
    return pharm_coords_list, bb_ents

# create linker pharmacophore model based on backbone
def create_bb_linker_pharm(dummy_pos_list, bb_pharm_coords, atom_pos_list, antisymbol, pharms, bb_ent, fake_lig, coords, bb_mol, bb_conf, dummy_count):
    if bb_ent[0] == 'Donor':
        for neighbor in bb_mol.GetAtomWithIdx(bb_ent[1][0]).GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                H_coords = list(bb_conf.GetAtomPosition(neighbor.GetIdx()))
        angle_list = []
        for dummy_pos in dummy_pos_list:
            angle, _ = cal_angle(H_coords, bb_pharm_coords, dummy_pos)
            angle_list.append(angle)
        min_dist, dist_list = cal_pharm_dummy_dist(bb_pharm_coords, dummy_pos_list)
        angle_dict = {}
        for i,v in enumerate(angle_list):
            angle_dict[i] = v
        angle_dict_order = sorted(angle_dict.items(), key=lambda x:x[1], reverse=True)
        if max(angle_list) >= 120 and min_dist <= 6.0:
            for item in angle_dict_order:
                if dummy_count[item[0]] == min(dummy_count):
                    dummy_count[item[0]] += 1
                    dummy_coords = dummy_pos_list[item[0]]
                    linker_pharm_coords = cal_point_coor(bb_pharm_coords, dummy_coords)
                    dummy_pos_list_new = cal_points_coor(linker_pharm_coords, dummy_coords)
                    if dummy_pos_list_new != []:
                        dummy_min_dist = cal_dis_coord(dummy_pos_list_new, atom_pos_list)[0]
                        pharm_min_dist = cal_dis_coord([linker_pharm_coords], atom_pos_list)[0]
                        if dummy_min_dist > 1.0 and pharm_min_dist > 1.0:
                            fake_lig += '.' + antisymbol[pharms.index(bb_ent[0])]
                            coords.append(linker_pharm_coords)
                    else:
                        fake_lig += '.' + antisymbol[pharms.index(bb_ent[0])]
                        coords.append(linker_pharm_coords)
                    break
        return coords, fake_lig

    elif bb_ent[0] == 'Acceptor':
        for neighbor in bb_mol.GetAtomWithIdx(bb_ent[1][0]).GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                C_coords = list(bb_conf.GetAtomPosition(neighbor.GetIdx()))
        angle_list = []
        for dummy_pos in dummy_pos_list:
            _, angle = cal_angle(C_coords, bb_pharm_coords, dummy_pos)
            angle_list.append(angle)
        min_dist, dist_list = cal_pharm_dummy_dist(bb_pharm_coords, dummy_pos_list)
        dist_dict = {}
        for i,v in enumerate(dist_list):
            dist_dict[i] = v
        dist_dict_order = sorted(dist_dict.items(), key=lambda x:x[1], reverse=False)
        if max(angle_list) >= 90 and min_dist <= 6.0:
            for item in dist_dict_order:
                if dummy_count[item[0]] == min(dummy_count):
                    dummy_count[item[0]] += 1
                    dummy_coords = dummy_pos_list[item[0]]
                    linker_pharm_coords = cal_point_coor(bb_pharm_coords, dummy_coords)
                    dummy_pos_list_new = cal_points_coor(linker_pharm_coords, dummy_coords)
                    if dummy_pos_list_new != []:
                        dummy_min_dist = cal_dis_coord(dummy_pos_list_new, atom_pos_list)[0]
                        pharm_min_dist = cal_dis_coord([linker_pharm_coords], atom_pos_list)[0]
                        if dummy_min_dist > 1.0 and pharm_min_dist > 1.0:
                            fake_lig += '.' + antisymbol[pharms.index(bb_ent[0])]
                            coords.append(linker_pharm_coords)
                    else:
                        fake_lig += '.' + antisymbol[pharms.index(bb_ent[0])]
                        coords.append(linker_pharm_coords)
                    break
        return coords, fake_lig

# obtain the pharmacophore model of the linker between two bricks
def obtain_linker_pharma(mol2_file, frag_sdf, output_sdf):
    mol = Chem.MolFromMol2File(mol2_file)
    conf = mol.GetConformer()
    try:
        frag_mol = Chem.SDMolSupplier(frag_sdf)[0]
    except:
        os.system('babel ' + frag_sdf + ' ' + frag_sdf)
        frag_mol = Chem.SDMolSupplier(frag_sdf)[0]
    frag_conf = frag_mol.GetConformer()
    pharms = ['Donor', 'Acceptor', 'Aromatic']
    feats = factory.GetFeaturesForMol(mol)
    new_feats = []
    for feat in feats:
        if feat.GetFamily() in pharms:
                new_feats.append([feat.GetFamily(), feat.GetAtomIds()])

    coords = []
    anchor_coords = []
    for i, atom in enumerate(frag_mol.GetAtoms()):
        if atom.GetAtomicNum() == 0:
            coords.append(list(frag_conf.GetAtomPosition(i)))
            anchor = atom.GetNeighbors()[0]
            anchor_coords.append(list(frag_conf.GetAtomPosition(anchor.GetIdx())))
    anchor_1, anchor_2 = anchor_coords[0], anchor_coords[1]

    dummy_pos_list = []
    atom_pos_list = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == '*':
            dummy_pos_list.append(list(conf.GetAtomPosition(atom.GetIdx())))
        else:
            atom_pos_list.append(list(conf.GetAtomPosition(atom.GetIdx())))
    dummy_count = [0] * len(dummy_pos_list)

    # Create map - aromatic is average position
    symbol = ['O', 'N', 'C'] # pharmacophore symbol on receptor
    antisymbol = ['N', 'O', 'C'] # pharmacophore symbol on linker
    fake_lig = '*.*'
    
    # # Uncomment below to Create receptor pharmacophore model (comment create linker pharmacophore model)
    # for ent in new_feats:
    #     tmp_coords = []
    #     tmp_idxs = []
    #     for idx in ent[1]:
    #         tmp_coords.append(list(conf.GetAtomPosition(idx)))
    #         tmp_idxs.append(idx)
    #     if len(tmp_idxs) > 0 :
    #         pharm_coords = list(np.average(tmp_coords, axis=0))
    #         fake_lig += '.' + symbol[pharms.index(ent[0])]
    #         coords.append(pharm_coords)

    # Create linker pharmacophore model
    sc_bb_file = obtain_sc_bb(mol2_file)
    # Create linker pharmacophore model according to side chain
    for sc_file in sc_bb_file[:-1]:
        res = os.path.splitext(sc_file)[0].split('_sc_')[1][0:3]
        sc_mol = Chem.MolFromMol2File(sc_file)
        if sc_mol is None:
            continue
        sc_conf = sc_mol.GetConformer()
        sc_feats = factory.GetFeaturesForMol(sc_mol)
        sc_new_feats = []
        for sc_feat in sc_feats:
            if sc_feat.GetFamily() in pharms:
                sc_new_feats.append([sc_feat.GetFamily(), sc_feat.GetAtomIds()])
        sc_dict = {}
        for sc_ent in sc_new_feats:
            for idx in sc_ent[1]:
                if idx in sc_dict:
                    sc_dict[idx].append(sc_ent[0])
                else:
                    sc_dict[idx] = [sc_ent[0]]

        if res == 'TRP':
            aro_idxs = []
            sc_copy_feats = copy.deepcopy(sc_new_feats)
            for sc_ent in sc_copy_feats:
                if sc_ent[0] == 'Aromatic':
                    for i in list(sc_ent[1]):
                        aro_idxs.append(i)
                    sc_new_feats.remove(sc_ent)
            sc_new_feats.append(['Aromatic', tuple(set(aro_idxs))])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, dummy_pos_list, anchor_1, anchor_2)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_linker_pharm(dummy_pos_list, sc_pharm_coords, atom_pos_list, antisymbol, pharms, sc_ent, fake_lig, coords, dummy_count)
        elif res == 'SER' or res == 'THR' or res == 'TYR':
            for key, val in sc_dict.items():
                if len(val) == 2:
                    sc_new_feats.remove(['Acceptor', (key, )])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, dummy_pos_list, anchor_1, anchor_2)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_linker_pharm(dummy_pos_list, sc_pharm_coords, atom_pos_list, antisymbol, pharms, sc_ent, fake_lig, coords, dummy_count)
        elif res == 'ARG':
            donor_idxs = []
            sc_copy_feats = copy.deepcopy(sc_new_feats)
            for sc_ent in sc_copy_feats:
                if sc_ent[0] == 'Donor':
                    donor_idxs.append(list(sc_ent[1])[0])
                    sc_new_feats.remove(sc_ent)
            sc_new_feats.append(['Donor', tuple(set(donor_idxs))])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, dummy_pos_list, anchor_1, anchor_2)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_linker_pharm(dummy_pos_list, sc_pharm_coords, atom_pos_list, antisymbol, pharms, sc_ent, fake_lig, coords, dummy_count)
        elif res == 'HIS':
            for key, val in sc_dict.items():
                 if len(val) == 3:
                     sc_new_feats.remove(['Donor', (key,)])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, dummy_pos_list, anchor_1, anchor_2)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_linker_pharm(dummy_pos_list, sc_pharm_coords, atom_pos_list, antisymbol, pharms, sc_ent, fake_lig, coords, dummy_count)
        else:
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, dummy_pos_list, anchor_1, anchor_2)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_linker_pharm(dummy_pos_list, sc_pharm_coords, atom_pos_list, antisymbol, pharms, sc_ent, fake_lig, coords, dummy_count)

    # Create linker pharmacophore model according to backbone
    bb_file = sc_bb_file[-1]
    bb_mol = Chem.MolFromMol2File(bb_file, removeHs=False)
    bb_conf = bb_mol.GetConformer()
    bb_feats = factory.GetFeaturesForMol(bb_mol)
    bb_new_feats = []
    for bb_feat in bb_feats:
        if bb_feat.GetFamily() in pharms:
            bb_new_feats.append([bb_feat.GetFamily(), bb_feat.GetAtomIds()])
    pharm_coords_list, bb_ents = obtain_bb_pharm_ent(bb_new_feats, bb_conf)
    for i in range(len(pharm_coords_list)):
        angle_1, angle_2 = cal_angle(anchor_1, anchor_2, pharm_coords_list[i])
        if angle_1 <= 95 and angle_2 <= 95:
            coords, fake_lig = create_bb_linker_pharm(dummy_pos_list, pharm_coords_list[i], atom_pos_list, antisymbol, pharms, bb_ents[i], fake_lig, coords, bb_mol, bb_conf, dummy_count)

    # Create object
    linker_pharm = Chem.MolFromSmiles(fake_lig)
    linker_conf = Chem.Conformer(linker_pharm.GetNumAtoms())
    for i in range(len(coords)):
        linker_conf.SetAtomPosition(i, coords[i])
    linker_conf.SetId(0)
    cid = linker_pharm.AddConformer(linker_conf)
    sdwriter = Chem.SDWriter(output_sdf)
    sdwriter.write(linker_pharm)
    sdwriter.close()

    # Obtain linker pharmacophore dict
    linker_feats = factory.GetFeaturesForMol(linker_pharm)
    dic = {'O': 0, 'N': 0, 'C': 0}
    for linker_atom in linker_pharm.GetAtoms():
        if linker_atom.GetSymbol() in symbol:
            dic[linker_atom.GetSymbol()] += 1
    return [dic['O'], dic['N'], dic['C']]


def main(brick1, brick2, atom_idx_1, atom_idx_2, dist, target, name, iter_idx):
    if not os.path.exists('dummy_sdf_' + str(iter_idx)):
        os.mkdir('dummy_sdf_' + str(iter_idx))
        
    if not os.path.exists('dummy_surround_mol_' + str(iter_idx)):
        os.mkdir('dummy_surround_mol_' + str(iter_idx))
        
    if not os.path.exists('data_path_' + str(iter_idx)):
        os.mkdir('data_path_' + str(iter_idx))

    if not os.path.exists('brickfolder_' + str(iter_idx+1)):
        os.mkdir('brickfolder_' + str(iter_idx+1))

    try:
        du_sdf, du_surround_mol = obtain_dummy(brick1, brick2, atom_idx_1, atom_idx_2, dist, target, name)
        min_dist = cal_dummy_dist(du_surround_mol)
    except Exception as e:
        print(str(e))
        os.system('rm *dummy*')
        print(brick1 + ' and ' + brick2 + ' are discarded due to no dummy line between them.')

    if min_dist >= 1.5:
        try:
            data_path, frag_sdf, dist = cal_dis_ang(brick1, brick2, atom_idx_1, atom_idx_2, name)
            if dist >= 2.0:
                linker_pharma_sdf = name + '_linker_pharma.sdf'
                pharm_count = obtain_linker_pharma(du_surround_mol, frag_sdf, linker_pharma_sdf)
                types = name + '.types'
                with open(types, 'w') as fw:
                    fw.write('1 data_path_' + str(iter_idx) + '/' + linker_pharma_sdf + ' data_path_' + str(iter_idx) + '/' + frag_sdf)
                
                raw_data = read_file_edit(data_path, pharm_count, add_idx=True, calc_pharm_counts=True)
                preprocess(raw_data, "zinc", os.path.splitext(data_path)[0], './', False)
                json_name = 'molecules_' + os.path.splitext(data_path)[0] + '.json'
                    
                shutil.move(du_sdf, 'dummy_sdf_' + str(iter_idx))
                shutil.move(du_surround_mol, 'dummy_surround_mol_' + str(iter_idx))
                shutil.move(frag_sdf, 'data_path_' + str(iter_idx))
                shutil.move(linker_pharma_sdf, 'data_path_' + str(iter_idx))
                shutil.move(types, 'data_path_' + str(iter_idx))
                shutil.move(json_name, 'data_path_' + str(iter_idx))
                os.system('rm *out.txt *_sc_* *_bb* *dummy*')
                print(brick1 + ' and ' + brick2 + ' are appropriate and preprocessed.')
            else:
                con_mol, con_sdf = connect_mol(brick1, brick2, atom_idx_1, atom_idx_2, name)
                os.remove(data_path)
                os.remove(frag_sdf)
                mol = Chem.SDMolSupplier(con_sdf)[0]
                molwt = Descriptors.MolWt(mol)
                if molwt <= 500:
                    os.system('mv ' + con_mol + ' brickfolder_' + str(iter_idx+1) + '/' + con_mol.split('_connect')[0] + '_screeningC.mol2')
                    print(brick1 + ' and ' + brick2 + ' are directly connected.')
                    os.remove(con_sdf)
                else:
                    os.remove(con_mol)
                    os.remove(con_sdf)
        except Exception as e:
            print(str(e))
            os.system('rm *dummy* *_sc_* *_bb* *out.txt *frag.sdf *linker_pharma.sdf *.types *.json')
            print(brick1 + ' and ' + brick2 + ' are discarded due to errors in data preparation process of linking.')

    else:
        os.system('rm *dummy*')
        print(brick1 + ' and ' + brick2 + ' are not appropriate.')


if __name__=="__main__":
    # Get the arguments
    # Usage: python data_pre_for_linking.py -l app_bricks.list -f brickfolder -t WDR5.pdb -i iter_idx
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-l', type=str, help='appropriate bricks list')
    parser.add_argument('-f', type=str, help='brick folder')
    parser.add_argument('-t', type=str, help='target protein')
    # parser.add_argument('-n', type=int, help='number of core')
    parser.add_argument('-i', type=int, help='iteration index')
    args = parser.parse_args()

    app_bricks = args.l
    brickfolder = args.f
    target = args.t
    # n = args.n
    iter_idx = args.i

    # start = time.time()
    # p = Pool(n)
    with open(app_bricks, 'r') as fr:
        lines = fr.readlines()
    for line in lines[1:]:
        newline = line.strip().split(', ')
        brick1, brick2, atom_idx_1, atom_idx_2, dist = brickfolder + '/' + newline[0], brickfolder + '/' + newline[1], int(newline[2]), int(newline[3]), float(newline[4])
        if iter_idx == 1:
            name = os.path.splitext(newline[0])[0].split('_')[2] + '_' + \
                    os.path.splitext(newline[0])[0].split('_')[4] + '_' + \
                    os.path.splitext(newline[1])[0].split('_')[2] + '_' + \
                    os.path.splitext(newline[1])[0].split('_')[4]
        else:
            name = os.path.splitext(newline[0])[0] + '_' + \
                    os.path.splitext(newline[1])[0].split('_')[2] + '_' + \
                    os.path.splitext(newline[1])[0].split('_')[4]
        try:
            main(brick1, brick2, atom_idx_1, atom_idx_2, dist, target, name, iter_idx)
        except:
            continue
    #     p.apply_async(main, args=(brick1, brick2, atom_idx_1, atom_idx_2, dist, target, name, iter_idx, ))
    # print('waiting for all processes')
    # p.close()
    # p.join()
    # end = time.time()
    # print('='*100)
    # print('Running time: %s Seconds' % (end-start))
