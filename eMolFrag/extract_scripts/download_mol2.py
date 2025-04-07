import os

with open('ZINC.wget', 'r') as fr:
	lines = fr.readlines()

for line in lines:
	os.system(line)
	print(line + ' is executed!')

print('=' * 100)
print('download done!')
