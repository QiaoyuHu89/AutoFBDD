import sys
import os
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")
import time
import math
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Geometry import rdGeometry
from collections import defaultdict, Counter
from utils import dataset_info
from data_pre_for_linking import delete_H

# sort the brick cluster based on the cluster size
def cluster_sort(file):
    with open(file, 'r') as fr:
        lines = fr.readlines()
    newlines = lines.copy()
    cluster_dict = {}
    for i in range(3, len(lines)-3):
        cluster_num = lines[i].strip().split()[0] + ' ' + lines[i].strip().split()[1]
        cluster_size = lines[i].strip().split()[3]
        cluster_dict[cluster_num] = int(cluster_size)
    cluster_dict_order = sorted(cluster_dict.items(), key=lambda x: x[1], reverse=True)
    for i in range(len(cluster_dict_order)):
        print(cluster_dict_order[i][0] + ': ' + str(cluster_dict_order[i][1]))
        for j in range(3, len(lines)-3):
            cluster_num = lines[j].strip().split()[0] + ' ' + lines[j].strip().split()[1]
            if cluster_num == cluster_dict_order[i][0]:
                newlines[i+3] = lines[j]
    newfile = file.split('.')[0] + '_order.' + file.split('.')[1]
    with open(newfile, 'w') as fw:
        fw.writelines(newlines)
    print(cluster_dict_order)
    return cluster_dict_order

# obtain the top k bricks in specified cluster
def obtain_top_bricks(cluster, k):
    with open(cluster, 'r') as fr:
        mol2_lines = fr.readlines()
        mol2_dict = defaultdict(list)
        for key,val in [(v,i) for i,v in enumerate(mol2_lines)]:
            mol2_dict[key].append(val)
        mol2_index = []
        mol2_name = []
        for mol2_line in mol2_lines:
            if mol2_line.startswith('# Name') and mol2_line not in mol2_name:
                mol2_name.append(mol2_line)
                if len(mol2_dict[mol2_line]) == 1:
                    mol2_index.append(mol2_dict[mol2_line][0])
                else:
                    for j in range(len(mol2_dict[mol2_line])):
                        mol2_index.append(mol2_dict[mol2_line][j])

    mol2_index.sort()
    bricks_list = []
    if k < len(mol2_index):
        for i in range(len(mol2_index)):
            brick = os.path.splitext(cluster)[0] + '_brick' + str(i+1) + '.mol2'
            with open(brick, 'w') as fw:
                for mol2_line in mol2_lines[mol2_index[i]: mol2_index[i+1]]:
                        fw.write(mol2_line)
            brick_noH_sdf = delete_H(brick)
            if brick_noH_sdf is not None:
                mol = Chem.SDMolSupplier(brick_noH_sdf)[0]
                for atom in mol.GetAtoms():
                    if "%s%i(%i)" % (atom.GetSymbol(), atom.GetTotalValence(), atom.GetFormalCharge()) not in dataset_info('zinc')['atom_types']:
                        os.system('rm ' + brick + ' ' + brick_noH_sdf)
                        print(brick + ' is deleted.')
                        break
                    elif atom.GetSymbol() == 'S':
                        neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
                        counter = Counter(neighbor_list)
                        if counter['O'] >= 3:
                            os.system('rm ' + brick + ' ' + brick_noH_sdf)
                            print(brick + ' is deleted.')
                            break
                else:
                    bricks_list.append(brick)
                    os.system('rm ' + brick_noH_sdf)
                if len(bricks_list) == k:
                    break
            else:
                os.system('rm ' + brick)
                print(brick + ' is deleted.')
    elif k >= len(mol2_index):
        print('Extract all the bricks in ' + cluster)
        for i in range(len(mol2_index)):
            if i < len(mol2_index) - 1:
                brick = os.path.splitext(cluster)[0] + '_brick' + str(i+1) + '.mol2'
                with open(brick, 'w') as fw:
                    for mol2_line in mol2_lines[mol2_index[i]: mol2_index[i+1]]:
                        fw.write(mol2_line)
                brick_noH_sdf = delete_H(brick)
                if brick_noH_sdf is not None:
                    mol = Chem.SDMolSupplier(brick_noH_sdf)[0]
                    for atom in mol.GetAtoms():
                        if "%s%i(%i)" % (atom.GetSymbol(), atom.GetTotalValence(), atom.GetFormalCharge()) not in dataset_info('zinc')['atom_types']:
                            os.system('rm ' + brick + ' ' + brick_noH_sdf)
                            print(brick + ' is deleted.')
                            break
                        elif atom.GetSymbol() == 'S':
                            neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
                            counter = Counter(neighbor_list)
                            if counter['O'] >= 3:
                                os.system('rm ' + brick + ' ' + brick_noH_sdf)
                                print(brick + ' is deleted.')
                                break
                    else:
                        bricks_list.append(brick)
                        os.system('rm ' + brick_noH_sdf)
                else:
                    os.system('rm ' + brick)
                    print(brick + ' is deleted.')
            elif i == len(mol2_index) - 1:
                brick = os.path.splitext(cluster)[0] + '_brick' + str(i+1) + '.mol2'
                with open(brick, 'w') as fw:
                    for mol2_line in mol2_lines[mol2_index[i]:]:
                        fw.write(mol2_line)
                brick_noH_sdf = delete_H(brick)
                if brick_noH_sdf is not None:
                    mol = Chem.SDMolSupplier(brick_noH_sdf)[0]
                    for atom in mol.GetAtoms():
                        if "%s%i(%i)" % (atom.GetSymbol(), atom.GetTotalValence(), atom.GetFormalCharge()) not in dataset_info('zinc')['atom_types']:
                            os.system('rm ' + brick + ' ' + brick_noH_sdf)
                            print(brick + ' is deleted.')
                            break
                        elif atom.GetSymbol() == 'S':
                            neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
                            counter = Counter(neighbor_list)
                            if counter['O'] >= 3:
                                os.system('rm ' + brick + ' ' + brick_noH_sdf)
                                print(brick + ' is deleted.')
                                break
                    else:
                        bricks_list.append(brick)
                        os.system('rm ' + brick_noH_sdf)
                else:
                    os.system('rm ' + brick)
                    print(brick + ' is deleted.')
    print(bricks_list)
    print('Top ' + str(k) + ' bricks in ' + cluster + ' is extracted and saved!' + '\n')
    return bricks_list

# calculate the distance of center of mass between two bricks
def cal_com(brick1, brick2):
    mol1 = Chem.MolFromMol2File(brick1, sanitize=False)
    mol2 = Chem.MolFromMol2File(brick2, sanitize=False)
    if mol1 is not None and mol2 is not None:
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        com1 = rdMolTransforms.ComputeCentroid(conf1)
        com2 = rdMolTransforms.ComputeCentroid(conf2)
        com_dis = rdGeometry.Point3D.Distance(com1, com2)
        com_dis = round(com_dis, 2)
        return com_dis
    else:
        return None

# calculate the nearest distance between two bricks
def cal_atom(brick1, brick2):
    mol1 = Chem.MolFromMol2File(brick1, sanitize=False)
    mol2 = Chem.MolFromMol2File(brick2, sanitize=False)
    if mol1 is not None and mol2 is not None:
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        mol1_coordinate = []
        for i in range(mol1.GetNumAtoms()):
            if mol1.GetAtomWithIdx(i).GetSymbol() != 'H':
                mol1_coordinate.append((conf1.GetAtomPosition(i)[0], conf1.GetAtomPosition(i)[1], conf1.GetAtomPosition(i)[2]))
        mol2_coordinate = []
        for j in range(mol2.GetNumAtoms()):
            if mol2.GetAtomWithIdx(j).GetSymbol() != 'H':
                mol2_coordinate.append((conf2.GetAtomPosition(j)[0], conf2.GetAtomPosition(j)[1], conf2.GetAtomPosition(j)[2]))
        
        atom_dis_list = []
        for coor1 in mol1_coordinate:
            for coor2 in mol2_coordinate:
                dis = 0
                for k in [0, 1, 2]:
                    dis += (coor1[k] - coor2[k])**2
                dis = round(math.sqrt(dis), 2)
                atom_dis_list.append(dis)
        atom_dis_list.sort()
        return atom_dis_list[0]
    else:
        return None

# calculate the nearest distance between atoms that can be connected within two bricks
def cal_atom_wH(brick1, brick2):
    mol1 = Chem.MolFromMol2File(brick1, removeHs=False, sanitize=False)
    mol2 = Chem.MolFromMol2File(brick2, removeHs=False, sanitize=False)
    if mol1 is not None and mol2 is not None:
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        atom1_idx = []
        for atom1 in mol1.GetAtoms():
            if 'H' in [x.GetSymbol() for x in atom1.GetNeighbors()]:
                atom1_idx.append(atom1.GetIdx())
                
        atom2_idx = []
        for atom2 in mol2.GetAtoms():
            if 'H' in [x.GetSymbol() for x in atom2.GetNeighbors()]:
                atom2_idx.append(atom2.GetIdx())
        
        mol1_coordinate = []
        for i in atom1_idx:
            mol1_coordinate.append((conf1.GetAtomPosition(i)[0], conf1.GetAtomPosition(i)[1], conf1.GetAtomPosition(i)[2]))
        mol2_coordinate = []
        for j in atom2_idx:
            mol2_coordinate.append((conf2.GetAtomPosition(j)[0], conf2.GetAtomPosition(j)[1], conf2.GetAtomPosition(j)[2]))
            
        atom_dis_list = []
        for coor1 in mol1_coordinate:
            for coor2 in mol2_coordinate:
                dis = 0
                for k in [0, 1, 2]:
                    dis += (coor1[k] - coor2[k])**2
                dis = round(math.sqrt(dis), 2)
                atom_dis_list.append(dis)
        atom_dis_list.sort()
        if atom_dis_list != []:
            return atom_dis_list[0]
        else:
            return None
    else:
        return None

# obtain the anchor points of two bricks
def obtain_anchors(brick1, brick2):
    mol1 = Chem.MolFromMol2File(brick1, removeHs=False, sanitize=False)
    mol2 = Chem.MolFromMol2File(brick2, removeHs=False, sanitize=False)
    if mol1 is not None and mol2 is not None:
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        atom1_idx = []
        mol1_atom = []
        for atom1 in mol1.GetAtoms():
            if 'H' in [x.GetSymbol() for x in atom1.GetNeighbors()]:
                atom1_idx.append(atom1.GetIdx())
            if atom1.GetSymbol() != 'H':
                mol1_atom.append(atom1.GetSymbol())
        mol1_atom_num = len(mol1_atom)
                
        atom2_idx = []
        for atom2 in mol2.GetAtoms():
            if 'H' in [x.GetSymbol() for x in atom2.GetNeighbors()]:
                atom2_idx.append(atom2.GetIdx())
        
        dis_val = 1000
        for i in atom1_idx:
            for j in atom2_idx:
                dis = 0
                for k in [0, 1, 2]:
                    dis += (conf1.GetAtomPosition(i)[k] - conf2.GetAtomPosition(j)[k])**2
                dis = round(math.sqrt(dis), 2)
                if dis < dis_val:
                    dis_val = dis
                    (anchor1, anchor2) = (i, j + mol1_atom_num)
        return (anchor1, anchor2)
    else:
        return (None, None)

# obtain the list of two appropriate bricks within threshold distance among every cluster
def obtain_app_bricks(bricks_list, val):
    app_bricks = []
    for brick1 in bricks_list:
        for brick2 in bricks_list:
            if brick1.split('_')[2] != brick2.split('_')[2]:
                #com_dis = cal_com(brick1, brick2)
                atom_dis = cal_atom(brick1, brick2)
                atom_wH_dis = cal_atom_wH(brick1, brick2)
                if atom_dis is not None and atom_wH_dis is not None:
                    if  atom_dis <= val and atom_wH_dis <= val and 1.0 <= atom_dis and 1.0 <= atom_wH_dis:
                        if (brick1, brick2) not in app_bricks and (brick2, brick1) not in app_bricks:
                            print(brick1 + ' and ' + brick2 + ' are within appropriate distance.')
                            app_bricks.append((brick1, brick2))

    if app_bricks != []:
        with open('app_bricks_1.list', 'w') as fw:
            fw.write('brick1, brick2, anchor1, anchor2, dis: \n')
            for bricks in app_bricks:
                (anchor1, anchor2) = obtain_anchors(bricks[0], bricks[1])
                atom_wH_dis = cal_atom_wH(bricks[0], bricks[1])
                if anchor1 is not None and anchor2 is not None and atom_wH_dis is not None:
                    fw.write(str(bricks[0]) + ', ' + str(bricks[1]) + ', ' + str(anchor1) + ', ' + str(anchor2) + ', ' + 
                            str(atom_wH_dis) + '\n')
    return app_bricks


if __name__=="__main__":
    # Usage: python ifitdock_postprocess.py -f cluster_out_cluster_info.list -t 10 -k 3 -v 10
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-f', type=str, help='cluster info list')
    parser.add_argument('-t', type=int, help='obtain top t clusters')
    parser.add_argument('-k', type=int, help='obtain top k bricks in each cluster')
    parser.add_argument('-v', type=float, help='threshold distance value of two bricks for brick linking')
    args = parser.parse_args()
    
    file = args.f
    t = args.t
    k = args.k
    val = args.v
    
    cluster_dict_order = cluster_sort(file)
    
    all_bricks_list = []
    for cluster in cluster_dict_order[:t]:
        cluster_mol2 = 'cluster_out_' + cluster[0].split()[1] + '_sorted.mol2'
        bricks_list = obtain_top_bricks(cluster_mol2, k)
        all_bricks_list.extend(bricks_list)
    # for brick in all_bricks_list:
    #     brick_pdb = brick.split('.')[0] + '.pdb'
    #     os.system('babel ' + brick + ' ' + brick_pdb)
    #     os.system('babel ' + brick_pdb + ' ' + brick)
    #     os.system('rm ' + brick_pdb)
    
    app_bricks = obtain_app_bricks(all_bricks_list, val)
    print(len(app_bricks))
