import os
import glob

folders = sorted(os.listdir('zinc_drug_like'))

num_1 = 0
for folder in folders:
	if os.path.isdir('zinc_drug_like/' + folder):
		num_2 = 0
		sub_folders = os.listdir('zinc_drug_like/' + folder)
		for sub_folder in sub_folders:
			mol2_list = glob.glob('zinc_drug_like/' + folder + '/' + sub_folder + '/*new*.mol2')
			num_2 += len(mol2_list)
			num_1 += len(mol2_list)
		print('number of new mol2 in ' + folder + ' is ' + str(num_2))
print('=' * 100)
print('counting finish')
print('total number of new mol2 is: ', num_1)
