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
from copy import deepcopy
from collections import Counter
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool
from data.prepare_data_scaffold_elaboration import read_file_edit, preprocess_edit
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

# calculate the distance between two points
def cal_points_dist(point1, point2):
    x1, y1, z1 = point1[0], point1[1], point1[2]
    x2, y2, z2 = point2[0], point2[1], point2[2]
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return dist

# calculate the coordinates of a series of points on a line given the coordinates of two end points
def cal_points_coor(point1, point2):
    x1, y1, z1 = point1[0], point1[1], point1[2]
    x2, y2, z2 = point2[0], point2[1], point2[2]
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    num_point = int(round(dist, 2) // 1.0 - 1)
    increment = (x2 - x1) / (num_point + 1)
    atoms_coor = []
    for i in range(1, num_point+1):
        x = x1 + increment*i
        y = ((x - x1)/(x2 - x1)) * (y2 - y1) + y1
        z = ((x - x1)/(x2 - x1)) * (z2 - z1) + z1
        atoms_coor.append([round(x, 4), round(y, 4), round(z, 4)])
    return atoms_coor

# calculate the coordinates of a middle points on a line given the coordinates of two end points
def cal_midpoint_coor(point1, point2):
    x1, y1, z1 = point1[0], point1[1], point1[2]
    x2, y2, z2 = point2[0], point2[1], point2[2]
    x = 0.5 * (x2 - x1) + x1
    y = 0.5 * (y2 - y1) + y1
    z = 0.5 * (z2 - z1) + z1
    atom_coor = [round(x, 4), round(y, 4), round(z, 4)]
    return atom_coor

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

# obtain the surrounding environment of a brick
def obtain_surround(brick, target, name):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()

    cmd.load(brick)
    cmd.load(target)
    cmd.create('whole', 'all')
    target_name = os.path.splitext(target)[0]
    cmd.delete(os.path.splitext(brick)[0])
    cmd.delete(target_name)
    cmd.select('ligand', 'organic')
    cmd.select('surrounding', 'byres ligand around 6')
    surround_mol = name + '_surround.mol2'
    cmd.save(surround_mol, 'surrounding')
    cmd.delete('all')

    return surround_mol

# calculate the minimium distance between the brick and surrounding res
def cal_brick_res_dist(brick, surround_mol):
    # brick_sdf = os.path.splitext(brick)[0] + '.sdf'
    # os.system('babel ' + brick + ' ' + brick_sdf)
    brick_sdf = delete_H(brick)
    if brick_sdf is not None:
        brick = Chem.SDMolSupplier(brick_sdf)[0]
        brick_conf = brick.GetConformer()
        surround = Chem.MolFromMol2File(surround_mol)
        surround_conf = surround.GetConformer()
        brick_pos_list = []
        atom_pos_list = []
        for atom in brick.GetAtoms():
            brick_pos_list.append(list(brick_conf.GetAtomPosition(atom.GetIdx())))
        for atom in surround.GetAtoms():
            atom_pos_list.append(list(surround_conf.GetAtomPosition(atom.GetIdx())))

        dist_list = []
        for i in brick_pos_list:
            for j in atom_pos_list:
                dist = 0
                for k in [0, 1, 2]:
                    dist += (i[k] - j[k])**2
                dist = round(math.sqrt(dist), 2)
                dist_list.append(dist)
        dist_list.sort()
        return dist_list[0]
    else:
        raise IOError('Hydrogens can not be deleted.')

# calculate the minimum distance between the dummy line and surrounding res
def cal_dummy_dist(dummy_pos_list, surround_mol):
    surround = Chem.MolFromMol2File(surround_mol)
    surround_conf = surround.GetConformer()
    atom_pos_list = []
    for atom in surround.GetAtoms():
        atom_pos_list.append(list(surround_conf.GetAtomPosition(atom.GetIdx())))

    dist_list = []
    for i in dummy_pos_list:
        for j in atom_pos_list:
            dist = 0
            for k in [0, 1, 2]:
                dist += (i[k] - j[k])**2
            dist = round(math.sqrt(dist), 2)
            dist_list.append(dist)
    dist_list.sort()
    return dist_list[0]

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
        sc_file = mol2_file.split('_surround')[0] + '_sc_' + item + '.mol2'
        cmd.select(item, 'sidechain and resi ' + item[3:6])
        cmd.save(sc_file, item)
        sc_bb_file.append(sc_file)
    bb_file = mol2_file.split('_surround')[0] + '_bb.mol2'
    cmd.select('bb', 'backbone')
    cmd.save(bb_file, 'bb')
    cmd.delete('all')
    sc_bb_file.append(bb_file)
    return sc_bb_file

# calculate the minimum distance and distance list between the brick and pharmacophore point
def cal_pharm_brick_dist(pharm_coords, brick_pos_list):
    dist_list = []
    for brick_pos in brick_pos_list:
        dist = 0
        for k in [0, 1, 2]:
            dist += (brick_pos[k] - pharm_coords[k])**2
        dist = round(math.sqrt(dist), 2)
        dist_list.append(dist)
    return min(dist_list), dist_list

# obtain side chain pharmacophore coordinates and feature environment
def obtain_sc_pharm_ent(sc_new_feats, sc_conf, brick_pos_list):
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
            pharm_coords_list.append(tmp_pharm_coords)
    if len(pharm_coords_list) > 1:
        for tmp_pharm_coords in pharm_coords_list:
            dist, _ = cal_pharm_brick_dist(tmp_pharm_coords, brick_pos_list)
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

# create pharmacophore model based on side chain
def create_sc_pharm(brick_pos_list, sc_pharm_coords, surround_mol, antisymbol, pharms, sc_ent, fake_lig, coords):
    min_dist, dist_list = cal_pharm_brick_dist(sc_pharm_coords, brick_pos_list)
    if min_dist >= 4.0 and min_dist <= 8.0:
        brick_coords = brick_pos_list[dist_list.index(min_dist)]
        pharm_coords = cal_points_coor(sc_pharm_coords, brick_coords)[2]
        dummy_pos_list = cal_points_coor(pharm_coords, brick_coords)
        if dummy_pos_list != []:
            dummy_min_dist = cal_dummy_dist(dummy_pos_list, surround_mol)
            pharm_min_dist = cal_dummy_dist([pharm_coords], surround_mol)
            if dummy_min_dist >= 1.5 and pharm_min_dist >= 1.5:
                fake_lig += '.' + antisymbol[pharms.index(sc_ent[0])]
                coords.append(pharm_coords)
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

# create pharmacophore model based on backbone
def create_bb_pharm(brick_pos_list, bb_pharm_coords, surround_mol, antisymbol, pharms, bb_ent, fake_lig, coords, bb_mol, bb_conf):
    if bb_ent[0] == 'Donor':
        for neighbor in bb_mol.GetAtomWithIdx(bb_ent[1][0]).GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                H_coords = list(bb_conf.GetAtomPosition(neighbor.GetIdx()))
        angle_list = []
        for brick_pos in brick_pos_list:
            angle, _ = cal_angle(H_coords, bb_pharm_coords, brick_pos)
            angle_list.append(angle)
        min_dist, dist_list = cal_pharm_brick_dist(bb_pharm_coords, brick_pos_list)
        if max(angle_list) >= 120 and min_dist >= 4.0 and min_dist <= 8.0:
            brick_coords = brick_pos_list[dist_list.index(min_dist)]
            pharm_coords = cal_points_coor(bb_pharm_coords, brick_coords)[2]
            dummy_pos_list = cal_points_coor(pharm_coords, brick_coords)
            if dummy_pos_list != []:
                dummy_min_dist = cal_dummy_dist(dummy_pos_list, surround_mol)
                pharm_min_dist = cal_dummy_dist([pharm_coords], surround_mol)
                if dummy_min_dist >= 1.5 and pharm_min_dist >= 1.5:
                    fake_lig += '.' + antisymbol[pharms.index(bb_ent[0])]
                    coords.append(pharm_coords)
        return coords, fake_lig

    elif bb_ent[0] == 'Acceptor':
        for neighbor in bb_mol.GetAtomWithIdx(bb_ent[1][0]).GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                C_coords = list(bb_conf.GetAtomPosition(neighbor.GetIdx()))
        angle_list = []
        for brick_pos in brick_pos_list:
            _, angle = cal_angle(C_coords, bb_pharm_coords, brick_pos)
            angle_list.append(angle)
        min_dist, dist_list = cal_pharm_brick_dist(bb_pharm_coords, brick_pos_list)
        if max(angle_list) >= 90 and min_dist >= 4.0 and min_dist <= 8.0:
            brick_coords = brick_pos_list[dist_list.index(min_dist)]
            pharm_coords = cal_points_coor(bb_pharm_coords, brick_coords)[2]
            dummy_pos_list = cal_points_coor(pharm_coords, brick_coords)
            if dummy_pos_list != []:
                dummy_min_dist = cal_dummy_dist(dummy_pos_list, surround_mol)
                pharm_min_dist = cal_dummy_dist([pharm_coords], surround_mol)
                if dummy_min_dist >= 1.5 and pharm_min_dist >= 1.5:
                    fake_lig += '.' + antisymbol[pharms.index(bb_ent[0])]
                    coords.append(pharm_coords)
        return coords, fake_lig

# obtain the pharmacophore model of the receptor around the brick
def obtain_pharma(brick, surround_mol, output_sdf, iter_idx):
    brick_sdf = os.path.splitext(brick)[0] + '.sdf'
    os.system('babel ' + brick + ' ' + brick_sdf)
    brick_mol = Chem.SDMolSupplier(brick_sdf, removeHs=False, sanitize=False)[0]
    brick_conf = brick_mol.GetConformer()
    # brick_sdf = delete_H(brick)
    # brick_mol = Chem.SDMolSupplier(brick_sdf)[0]
    # brick_conf = brick_mol.GetConformer()
    surround = Chem.MolFromMol2File(surround_mol)
    surround_conf = surround.GetConformer()
    pharms = ['Donor', 'Acceptor', 'Aromatic']
    feats = factory.GetFeaturesForMol(surround)
    new_feats = []
    for feat in feats:
        if feat.GetFamily() in pharms:
                new_feats.append([feat.GetFamily(), feat.GetAtomIds()])

    coords = []
    brick_pos_list = []
    for atom in brick_mol.GetAtoms():
        brick_pos_list.append(list(brick_conf.GetAtomPosition(atom.GetIdx())))

    # Create map - aromatic is average position
    symbol = ['O', 'N', 'C'] # pharmacophore symbol on receptor
    antisymbol = ['N', 'O', 'C'] # pharmacophore symbol on linker
    fake_lig = ''

    # # Uncomment below to create receptor pharmacophore model (comment create pharmacophore model)
    # for ent in new_feats:
    #     tmp_coords = []
    #     tmp_idxs = []
    #     for idx in ent[1]:
    #         tmp_coords.append(list(surround_conf.GetAtomPosition(idx)))
    #         tmp_idxs.append(idx)
    #     if len(tmp_idxs) > 0 :
    #         pharm_coords = list(np.average(tmp_coords, axis=0))
    #         fake_lig += '.' + symbol[pharms.index(ent[0])]
    #         coords.append(pharm_coords)

    # Create pharmacophore model
    sc_bb_file = obtain_sc_bb(surround_mol)
    # Create pharmacophore model according to side chain
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
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, brick_pos_list)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_pharm(brick_pos_list, sc_pharm_coords, surround_mol, antisymbol, pharms, sc_ent, fake_lig, coords)
        elif res == 'SER' or res == 'THR' or res == 'TYR':
            for key, val in sc_dict.items():
                if len(val) == 2:
                    sc_new_feats.remove(['Acceptor', (key, )])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, brick_pos_list)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_pharm(brick_pos_list, sc_pharm_coords, surround_mol, antisymbol, pharms, sc_ent, fake_lig, coords)
        elif res == 'ARG':
            donor_idxs = []
            sc_copy_feats = copy.deepcopy(sc_new_feats)
            for sc_ent in sc_copy_feats:
                if sc_ent[0] == 'Donor':
                    donor_idxs.append(list(sc_ent[1])[0])
                    sc_new_feats.remove(sc_ent)
            sc_new_feats.append(['Donor', tuple(set(donor_idxs))])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, brick_pos_list)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_pharm(brick_pos_list, sc_pharm_coords, surround_mol, antisymbol, pharms, sc_ent, fake_lig, coords)
        elif res == 'HIS':
            for key, val in sc_dict.items():
                if len(val) == 3:
                    sc_new_feats.remove(['Donor', (key,)])
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, brick_pos_list)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_pharm(brick_pos_list, sc_pharm_coords, surround_mol, antisymbol, pharms, sc_ent, fake_lig, coords)
        else:
            sc_pharm_coords, sc_ent = obtain_sc_pharm_ent(sc_new_feats, sc_conf, brick_pos_list)
            if sc_pharm_coords is not None and sc_ent is not None:
                coords, fake_lig = create_sc_pharm(brick_pos_list, sc_pharm_coords, surround_mol, antisymbol, pharms, sc_ent, fake_lig, coords)

    # Create pharmacophore model according to backbone
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
        coords, fake_lig = create_bb_pharm(brick_pos_list, pharm_coords_list[i], surround_mol, antisymbol, pharms, bb_ents[i], fake_lig, coords, bb_mol, bb_conf)

    # Update pharmacophore model by clustering
    fake_lig = fake_lig[1:]
    coords_group = []
    for coord in coords:
        coords_group.append([coord])
    tmp_coords_list = deepcopy(coords_group)
    while tmp_coords_list != []:
        tmp_coords_list_copy = deepcopy(tmp_coords_list)
        tmp_coords_list = []
        for coord in coords:
            for tmp_coords in tmp_coords_list_copy:
                tmp_num = 0
                if coord not in tmp_coords:
                    for tmp_coord in tmp_coords:
                        dist = cal_points_dist(coord, tmp_coord)
                        if dist <= 5.0:
                            tmp_num += 1
                    if tmp_num == len(tmp_coords):
                        tmp_coords_copy = deepcopy(tmp_coords)
                        tmp_coords_copy.append(coord)
                        tmp_coords_list.append(tmp_coords_copy)
        coords_group.extend(tmp_coords_list)
    coords_group = [sorted(i) for i in coords_group]
    new_coords_group = []

    if coords_group != []:
        for coord_group in coords_group:
            if coord_group not in new_coords_group:
                new_coords_group.append(coord_group)
        new_coords_group = sorted(new_coords_group, key = lambda i: len(i), reverse=True)
        new_coords_len = len(new_coords_group[0])
        large_coords_group = [new_coord_group for new_coord_group in new_coords_group if len(new_coord_group) == new_coords_len]
        large_coords_dict = {}
        centroid_brick = np.mean(np.array(brick_pos_list), axis=0).tolist()
        for i in range(len(large_coords_group)):
            centroid = np.mean(np.array(large_coords_group[i]), axis=0).tolist()
            dist = cal_points_dist(centroid, centroid_brick)
            # dist, _ = cal_pharm_brick_dist(centroid, brick_pos_list)
            large_coords_dict[i] = dist
        sorted_large_coords_dict = sorted(large_coords_dict.items(), key=lambda i: i[1])
        final_coords = large_coords_group[sorted_large_coords_dict[0][0]]
        final_centroid = np.mean(np.array(final_coords), axis=0).tolist()
        _, dist_list = cal_pharm_brick_dist(final_centroid, brick_pos_list)
        sorted_dist_list = sorted(dist_list)
        for dist in sorted_dist_list:
            atom = brick_mol.GetAtomWithIdx(dist_list.index(dist))
            # if atom.GetSymbol() == 'C' and atom.GetExplicitValence() != 4:
            if 'H' in [x.GetSymbol() for x in atom.GetNeighbors()]:
                anchor_index = dist_list.index(dist)
                anchor_point = brick_pos_list[anchor_index]
                final_coords.append(anchor_point)
                with open('app_bricks_' + str(iter_idx) + '.list', 'a') as fw:
                    fw.write(brick.split('/')[1] + ', ' + str(dist) + '\n')
                break
        final_fake_lig = ''
        fake_lig_list = fake_lig.split('.')
        for final_coord in final_coords[:-1]:
            final_fake_lig += '.' + fake_lig_list[coords.index(final_coord)]
        final_fake_lig += '.*'
        final_fake_lig = final_fake_lig[1:]

        # Create object
        pharm = Chem.MolFromSmiles(final_fake_lig)
        conf = Chem.Conformer(pharm.GetNumAtoms())
        for i in range(len(final_coords)):
            conf.SetAtomPosition(i, final_coords[i])
        conf.SetId(0)
        cid = pharm.AddConformer(conf)
        sdwriter = Chem.SDWriter(output_sdf)
        sdwriter.write(pharm)
        sdwriter.close()

        # Obtain pharmacophore dict
        dic = {'O': 0, 'N': 0, 'C': 0}
        for atom in pharm.GetAtoms():
            if atom.GetSymbol() in symbol:
                dic[atom.GetSymbol()] += 1
        return [dic['O'], dic['N'], dic['C']], anchor_index
    else:
        raise IOError('Pharmacophore coordinates group list is empty!.')

# obtain brick sdf file with anchor point
def obtain_frag(brick, anchor_index, frag_sdf):
    brick_sdf = delete_H(brick)
    try:
        brick = Chem.SDMolSupplier(brick_sdf)[0]
        brick_conf = brick.GetConformer()
        brick_pos_list = []
        for atom in brick.GetAtoms():
            brick_pos_list.append(list(brick_conf.GetAtomPosition(atom.GetIdx())))
        
        combo = Chem.CombineMols(brick, Chem.MolFromSmiles("*"))
        edcombo = Chem.EditableMol(combo)
        num_heavy_atoms = combo.GetNumHeavyAtoms()
        edcombo.AddBond(num_heavy_atoms, anchor_index, order=Chem.rdchem.BondType.SINGLE)
        mol_to_link = edcombo.GetMol()
        Chem.SanitizeMol(mol_to_link)

        # Convert exit vectors to carbons for conformer generation
        du = Chem.MolFromSmiles('*')
        mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link, du, Chem.MolFromSmiles('C'), True)[0]
        Chem.SanitizeMol(mol_to_link_carbon)
        # Generate conformer
        mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
        AllChem.ConstrainedEmbed(mol_to_link_carbon, brick, randomseed=42)
        mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)

        # Add this conformer to the two unlinked fragments
        conf = mol_to_link.GetConformer()
        ref_conf = mol_to_link_carbon.GetConformer()
        for i in range(mol_to_link_carbon.GetNumAtoms()):
            pos = list(ref_conf.GetAtomPosition(i))
            conf.SetAtomPosition(i, pos)
        conf.SetId(0)
        _ = mol_to_link.AddConformer(conf)

        w = Chem.SDWriter(frag_sdf)
        w.write(mol_to_link)
        w.close()
    except:
        raise IOError('Frags with anchor can not be generated!')

def main(brick, target, name, iter_idx):
    if not os.path.exists('surround_mol_' + str(iter_idx)):
        os.mkdir('surround_mol_' + str(iter_idx))
        
    if not os.path.exists('data_path_' + str(iter_idx)):
        os.mkdir('data_path_' + str(iter_idx))

    brick_sdf = name + '_frag.sdf'
    try:
        surround_mol = obtain_surround(brick, target, name)
        min_dist = cal_brick_res_dist(brick, surround_mol)
    except:
        os.system('rm *_surround*')
        print(brick + ' is discarded due to no surround.')

    if min_dist >= 1.5:
        try:
            pharma_sdf = name + '_pharma.sdf'
            pharm_count, anchor_index = obtain_pharma(brick, surround_mol, pharma_sdf, iter_idx)
            obtain_frag(brick, anchor_index, brick_sdf)
            data_path = name + '_out.txt'
            with open(data_path, 'w') as f:
                f.write("%s" % (Chem.MolToSmiles(Chem.SDMolSupplier(brick_sdf)[0])))
            types = name + '.types'
            with open(types, 'w') as fw:
                fw.write('1 data_path_' + str(iter_idx) + '/' + brick_sdf + ' data_path_' + str(iter_idx) + '/' + pharma_sdf)

            raw_data = read_file_edit(data_path, pharm_count, add_idx=True, calc_pharm_counts=True)
            preprocess_edit(raw_data, "zinc", os.path.splitext(data_path)[0], './', False)
            json_name = 'molecules_' + os.path.splitext(data_path)[0] + '.json'

            shutil.move(surround_mol, 'surround_mol_' + str(iter_idx))
            os.system('rm ' + data_path)
            shutil.move(brick_sdf, 'data_path_' + str(iter_idx))
            shutil.move(pharma_sdf, 'data_path_' + str(iter_idx))
            shutil.move(types, 'data_path_' + str(iter_idx))
            shutil.move(json_name, 'data_path_' + str(iter_idx))
            os.system('rm *out.txt *_sc_* *_bb* *_surround*')
            print(brick + ' is appropriate and preprocessed.')
        except Exception as e:
            print(str(e))
            os.system('rm *_surround* *_sc_* *_bb* *out.txt *frag.sdf *pharma.sdf *.types *.json')
            print(brick + ' is discarded due to errors in data preparation process of growing.')

    else:
        os.system('rm *_surround*')
        print(brick + ' is not appropriate.')


if __name__=="__main__":
    # Get the arguments
    # Usage: python data_pre_for_growing.py -f brickfolder -t WDR5.pdb -i iter_idx
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-f', type=str, help='brick folder')
    parser.add_argument('-t', type=str, help='target protein')
    # parser.add_argument('-n', type=int, help='number of core')
    parser.add_argument('-i', type=int, help='iteration index')
    args = parser.parse_args()

    brickfolder = args.f
    target = args.t
    # n = args.n
    iter_idx = args.i

    # start = time.time()
    with open('app_bricks_' + str(iter_idx) + '.list', 'a') as fw:
        fw.write('brick, dis: \n')
    # p = Pool(n)
    brick_list = os.listdir(brickfolder)
    for brick in brick_list:
        brick = brickfolder + '/' + brick
        if iter_idx == 1:
            name = os.path.splitext(brick)[0].split('/')[1].split('_')[2] + '_' + os.path.splitext(brick)[0].split('/')[1].split('_')[4]
        else:
            name = os.path.splitext(brick)[0].split('/')[1]
        try:
            main(brick, target, name, iter_idx)
        except:
            continue
        # p.apply_async(main, args=(brick, target, name, iter_idx, ))
    # print('waiting for all processes')
    # p.close()
    # p.join()
    # end = time.time()
    # print('='*100)
    # print('Running time: %s Seconds' % (end-start))