import os
import glob
from multiprocessing import Pool, cpu_count
import time
from collections import defaultdict

def extract(mol2_file):
	# gz_file = line.strip().split(' ')[-1]
	# os.system('gunzip ' + gz_file)
	# mol2_file = line.strip().split(' ')[-1].replace('.gz', '')
	# new_mol2_file = mol2_file.split('.')[0] + '.' + mol2_file.split('.')[1] + '_new.mol2'
	# os.system('babel ' + mol2_file + ' ' + new_mol2_file + ' -m')
	with open(mol2_file, 'r') as fr:
		mol2_lines = fr.readlines()
		mol2_dict = defaultdict(list)
		for k,va in [(v,i) for i,v in enumerate(mol2_lines)]:
			mol2_dict[k].append(va)
		mol2_index = []
		mol2_zinc = []
		for mol2_line in mol2_lines:
			if mol2_line.startswith('ZINC') and mol2_line not in mol2_zinc:
				mol2_zinc.append(mol2_line)
				if len(mol2_dict[mol2_line]) == 1:
					mol2_index.append(mol2_dict[mol2_line][0])
				else:
					for j in range(len(mol2_dict[mol2_line])):
						mol2_index.append(mol2_dict[mol2_line][j])
		# print(mol2_index)
		# print(len(mol2_index))
		mol2_num = len(mol2_index)
	for i in range(mol2_num):
		if i < mol2_num - 1:
			with open(mol2_file.split('.')[0] + '.' + mol2_file.split('.')[1] + '_new' + str(i+1) + '.mol2', 'w') as fw:
				for mol2_line in mol2_lines[mol2_index[i]-1: mol2_index[i+1]-1]:
					fw.write(mol2_line)
		elif i == mol2_num - 1:
			with open(mol2_file.split('.')[0] + '.' + mol2_file.split('.')[1] + '_new' + str(i+1) + '.mol2', 'w') as fw:
				for mol2_line in mol2_lines[mol2_index[i]-1:]:
					fw.write(mol2_line)
	print(mol2_file + ' is extracted and separated!')


start = time.time()
p = Pool(10)
# with open('ZINC.wget', 'r') as fr:
# 	lines = fr.readlines()
mol2_files = glob.glob('JA/*/*.mol2')

# for line in lines:
for mol2_file in mol2_files:
	p.apply_async(extract, args=(mol2_file, ))
print('waiting for all processes')
p.close()
p.join()
end = time.time()
print("Total time {} s".format((end - start)))
print('=' * 100)
print('extraction done!')
