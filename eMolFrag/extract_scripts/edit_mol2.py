import os
import glob
from multiprocessing import Pool
import time

def edit_mol2(folder):
	mol2_files = glob.glob('zinc_drug_like/' + folder + '/*/*new*.mol2')

	for mol2_file in mol2_files:
		with open(mol2_file, 'r') as fr:
			mol2_lines = fr.readlines()
		os.system('rm ' + mol2_file)
		with open(mol2_file, 'a') as fw:
			fw.writelines(mol2_lines[:-2])
		print('Edit ' + mol2_file)

start = time.time()
p = Pool(20)
folders = os.listdir('zinc_drug_like')
for folder in folders:
	p.apply_async(edit_mol2, args=(folder, ))
print('waiting for all processes')
p.close()
p.join()
end = time.time()
print("Total time {} s".format((end - start)))
print('=' * 100)
print('edit done!')
