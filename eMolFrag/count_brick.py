import os
import glob

folders = os.listdir('zinc_drug_brick')

num_1 = 0
num_2 = 0
for folder in folders:
        sub_num_1 = 0
        sub_folders = os.listdir('zinc_drug_out/' + folder)
        for sub_folder in sub_folders:
                brick_list = glob.glob('zinc_drug_out/' + folder + '/' + sub_folder + '/output-brick/*.sdf')
                if len(brick_list) != 0:
                        sub_num_1 += len(brick_list)
        num_1 += sub_num_1
        print('number of brick in ' + folder + ' before removing redundancy is ' + str(sub_num_1))
        sub_num_2 = 0
        brick_list = glob.glob('zinc_drug_brick/' + folder  + '/output-brick/*.sdf')
        if len(brick_list) != 0:
                sub_num_2 += len(brick_list)
        num_2 += sub_num_2
        print('number of brick in ' + folder + ' after removing redundancy is ' + str(sub_num_2))
        print('-' * 50)
print('=' * 100)
print('counting finish')
print('total number of brick before removing redundancy is: ', num_1)
print('total number of brick after first round of removing redundancy is: ', num_2)
