import os
from multiprocessing import Pool, cpu_count
import time

def delete(line):
	mol2_file = line.strip().split(' ')[-1].replace('.gz', '')
	new_mol2_file = mol2_file.split('.')[0] + '.' + mol2_file.split('.')[1] + '_new*.mol2'
	os.system('rm ' + new_mol2_file)
	print(new_mol2_file + ' is deleted!')


start = time.time()
p = Pool(28)
with open('ZINC-downloader-3D-mol2.gz.wget', 'r') as fr:
	lines = fr.readlines()

for line in lines:
	p.apply_async(delete, args=(line,))
print('waiting for all processes')
p.close()
p.join()
end = time.time()
print("Total time {} s".format((end - start)))	
print('=' * 100)
print('deletion done!')
