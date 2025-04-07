import os
import re
import sys
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")
import glob
import argparse
from ifitdock_postprocess import *


def main(iter_idx, val):
    app_bricks = []
    if os.path.exists('brickfolder_' + str(iter_idx+1)):
        os.system('cp brickfolder_ori/* brickfolder_' + str(iter_idx+1))
        screening_mols = glob.glob('brickfolder_' + str(iter_idx+1) + '/*screening*')
        screeningC_mols = glob.glob('brickfolder_' + str(iter_idx+1) + '/*screeningC*')
        new_screening_mols = [mol for mol in screening_mols if mol not in screeningC_mols]
        screening_dict = {}
        for screening_mol in new_screening_mols:
            with open(screening_mol, 'r') as fr:
                lines = fr.readlines()
            for line in lines:
                if line.startswith('delta_G'):
                    delta_G = float(line.split(' ')[2])
            screening_dict[screening_mol] = delta_G
        screening_list = sorted(screening_dict.items(), key=lambda x: x[1], reverse=False)
        if len(screening_list) > 100:
            screening_list = screening_list[0:100]
        else:
            pass
        new_screening_list = []
        for screening_mol in screening_list:
            new_screening_list.append(screening_mol[0])

        bricks_list_1 = list(set(new_screening_list + screeningC_mols))
        bricks_list_2 = glob.glob('brickfolder_' + str(iter_idx+1) + '/cluster*')
        if bricks_list_1 != [] and bricks_list_2 != []:
            for brick1 in bricks_list_1:
                for brick2 in bricks_list_2:
                    cluster_list_1 = re.findall(r"\b\d+\b", ', '.join(brick1.split('/')[1].split('_')))
                    cluster_list_2 = re.findall(r"\b\d+\b", ', '.join(brick2.split('/')[1].split('_')))
                    if cluster_list_2[0] not in cluster_list_1:
                        # com_dis = cal_com(brick1, brick2)
                        atom_dis = cal_atom(brick1, brick2)
                        atom_wH_dis = cal_atom_wH(brick1, brick2)
                        if atom_dis is not None and atom_wH_dis is not None:
                            if atom_dis <= val and atom_wH_dis <= val and 1.0 <= atom_dis and 1.0 <= atom_wH_dis:
                                if (brick1, brick2) not in app_bricks and (brick2, brick1) not in app_bricks:
                                    print(brick1 + ' and ' + brick2 + ' are within appropriate distance.')
                                    app_bricks.append((brick1, brick2))

    if app_bricks != []:
        with open('app_bricks_' + str(iter_idx+1) + '.list', 'w') as fw:
            fw.write('brick1, brick2, anchor1, anchor2, dis: \n')
            for bricks in app_bricks:
                (anchor1, anchor2) = obtain_anchors(bricks[0], bricks[1])
                atom_wH_dis = cal_atom_wH(bricks[0], bricks[1])
                if anchor1 is not None and anchor2 is not None and atom_wH_dis is not None:
                    fw.write(str(bricks[0].split('/')[1]) + ', ' + str(bricks[1].split('/')[1]) + ', ' + str(anchor1) + ', ' + 
                            str(anchor2) + ', ' + str(atom_wH_dis) + '\n')
        print(len(app_bricks))
    else:
        print('There are no appropriate bricks in iteration ' + str(iter_idx) + '. AutoFBDD exits!')


if __name__=="__main__":    
    # Usage: python select_linked_bricks.py -i iter_idx -v 10
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-i', type=int, help='iteration index')
    parser.add_argument('-v', type=float, help='threshold distance value of two bricks for brick linking')
    args = parser.parse_args()
    
    iter_idx = args.i
    val = args.v
    
    main(iter_idx, val)
    print('Selection linked bricks done!')
