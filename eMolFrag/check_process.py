import os
import glob

root = 'FH'
folders = os.listdir('zinc_drug_like_out/' + root)
folders.sort()

for folder in folders:
	process_file = 'zinc_drug_like_out/' + root + '/' + folder + '/output-log/Process.log'
	with open(process_file, 'r') as fr:
		lines = fr.readlines()
	print(root + '/' + folder + ': ' + lines[-1].strip())
print('=' * 100)
print('checking finish')
