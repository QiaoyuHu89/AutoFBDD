import os
import time
import glob

folders = sorted(os.listdir('zinc_drug_like_out'))

num_1 = 0
for folder in folders:
	num_2 = 0
	sub_folders = os.listdir('zinc_drug_like_out/' + folder)
	for sub_folder in sub_folders:
		sdf_list = glob.glob('zinc_drug_like_out/' + folder + '/' + sub_folder + '/output-chop-comb/b*.sdf')
		num_2 += len(sdf_list)
		num_1 += len(sdf_list)
	print('number of brick sdf in ' + folder + ' is ' + str(num_2))
print('=' * 100)
print('counting finish')
print('time:', time.asctime(time.localtime(time.time())))
print('total number of brick sdf is: ', num_1)
