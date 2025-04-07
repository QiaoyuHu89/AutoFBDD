import os
import shutil
from rdkit import Chem

with open('obj-KRASG12D-inhi.smi', 'r') as fr:
    line_list = fr.readlines()

if not os.path.exists('KRAS_G12D_smi'):
    os.mkdir('KRAS_G12D_smi')

for line in line_list:
    smi = line.split('\t')[0]
    name = line.split('\t')[1].split('\n')[0]
    # name_sdf = name + '.sdf'
    # mol = Chem.MolFromSmiles(smi)
    # writer = Chem.SDWriter(name_sdf)
    # writer.write(mol)
    name_smi = name + '.smi'
    with open(name_smi, 'w') as fw:
        fw.write(smi)

    # os.system('babel ' + name + '.sdf ' + name + '.mol2')
    # os.remove(name + '.sdf')
    shutil.move(name_smi, 'KRAS_G12D_smi')
    print(name + ' is extracted!')