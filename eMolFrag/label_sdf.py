import os
import shutil

name_list = os.listdir('KRAS_G12D_sdf')
os.chdir('KRAS_G12D_sdf')
for name in name_list:
    # print(name)
    # with open(name, 'a+') as fa:
    #     fa.write(name)
    name_mol2 = name.split('.')[0] + '.mol2'
    os.system('babel ' + name + ' ' + name_mol2)
    shutil.move(name_mol2, '../KRAS_G12D_mol2')
    print('Convert ' + name_mol2 + ' successfully!')
