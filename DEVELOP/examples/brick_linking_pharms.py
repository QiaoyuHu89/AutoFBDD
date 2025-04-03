import os
import sys
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")

import math
import time
import glob
import shutil
import argparse
from multiprocessing import Pool
from collections import defaultdict
import DEVELOP
from DEVELOP import DenseGGNNChemModel


def get_data_dict(app_bricks, data_path, iter_idx):
    with open(app_bricks, 'r') as fr:
        lines = fr.readlines()

    types_list = glob.glob(data_path + '/*.types')
    types_list.sort()
    data_dict = {}
    for types in types_list:
        types_name = os.path.splitext(types)[0].split('/')[-1]
        for line in lines[1:]:
            newline = line.strip().split(', ')
            if iter_idx == 1:
                bricks_name = os.path.splitext(newline[0])[0].split('_')[2] + '_' + \
                        os.path.splitext(newline[0])[0].split('_')[4] + '_' + \
                        os.path.splitext(newline[1])[0].split('_')[2] + '_' + \
                        os.path.splitext(newline[1])[0].split('_')[4]
            else:
                bricks_name = os.path.splitext(newline[0])[0] + '_' + \
                        os.path.splitext(newline[1])[0].split('_')[2] + '_' + \
                        os.path.splitext(newline[1])[0].split('_')[4]

            if types_name == bricks_name:
                dist = newline[-1]
                data_dict[types_name] = dist
    return data_dict

def linking_bricks(json_file, types_file, name, min, max, iter_idx):
    args = defaultdict(None)
    args['--dataset'] = 'zinc'
    args['--config'] = '{"generation": true, \
                        "batch_size": 1, \
                        "number_of_generation_per_valid": 50, \
                        "train_file": "' + json_file + '", \
                        "valid_file": "' + json_file + '", \
                        "train_struct_file": "' + types_file + '", \
                        "valid_struct_file": "' + types_file + '", \
                        "struct_data_root": "./", \
                        "output_name": "' + name + '", \
                        "struct_data_len": 9, \
                        "min_atoms": ' + min + ', \
                        "max_atoms": ' + max + '}'
    args['--freeze-graph-model'] = False
    args['--restore'] = './pretrained_DEVELOP_model_pharms.pickle'
    # Setup model and generate molecules
    model = DenseGGNNChemModel(args)
    model.train()
    # Free up some memory
    model = ''
    os.system('rm *generated_smiles_zinc')
    os.system('rm *params_zinc.json')
    print(name + '_generated_smis.smi is generated!')
    os.system('mv ' + name + '*.smi generated_smi_' + str(iter_idx))

def main(app_bricks, data_path, iter_idx):
    data_dict = get_data_dict(app_bricks, data_path, iter_idx)
    for key, val in data_dict.items():
        json_file = data_path + '/molecules_' + key + '_out.json'
        types_file = data_path + '/' + key + '.types'
        with open(json_file, 'r') as fr:
            lines = fr.readlines()
        if lines[0] != str([]):
            if round(float(val)) < 2:
                min = str(0)
                max = str(round(float(val))+2)
            else:
                min = str(round(float(val))-2)
                max = str(round(float(val))+2)
            linking_bricks(json_file, types_file, key, min, max, iter_idx)

if __name__=="__main__":
    # Get the arguments
    # Usage: python brick_linking.py -l app_bricks.list -d data_path -i iter_idx
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-l', type=str, help='appropriate bricks list')
    parser.add_argument('-d', type=str, help='data path folder')
    parser.add_argument('-i', type=int, help='iteration index')
    args = parser.parse_args()

    app_bricks = args.l
    data_path = args.d
    iter_idx = args.i
    
    if not os.path.exists('generated_smi_' + str(iter_idx)):
        os.mkdir('generated_smi_' + str(iter_idx))

    start = time.time()
    main(app_bricks, data_path, iter_idx)
    end = time.time()
    print('='*100)
    print('Running time: %s Seconds' % (end-start))
    print('brick linking done!')
