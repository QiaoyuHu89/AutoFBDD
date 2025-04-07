import os
import glob
import sys
import shutil

args = sys.argv
folder = args[1]
outputfolder = args[2]
outputlog = args[3]

if not os.path.exists('zinc_drug_brick/' + folder):
    os.mkdir('zinc_drug_brick/' + folder)

if not os.path.exists(outputfolder):
    os.mkdir(outputfolder)
    
if not os.path.exists(outputlog):
    os.mkdir(outputlog)

brick_list = glob.glob('zinc_drug_out/' + folder + '*/output-brick/*.sdf')
for brick in brick_list:
    brick_name = os.path.basename(brick)
    new_brick = brick.replace('output-brick', 'output-chop-comb')
    shutil.copy(new_brick, outputfolder)
    BrickListAll = brick.split('/')[0] + '/' + brick.split('/')[1] + '/' + brick.split('/')[2] + '/output-log/BrickListAll.txt'
    with open(BrickListAll, 'r') as fr:
        lines = fr.readlines()
    for line in lines:
        if line.strip().split(' ')[0].split('/')[-1] == brick_name:
            with open(outputlog + 'BrickListAll.txt', 'a') as fw:
                new_line = line.strip().split(' ')[0].split('/')[0] + '/' + line.strip().split(' ')[0].split('/')[1] + '/' + \
                        line.strip().split(' ')[0].split('/')[2] + '/' + line.strip().split(' ')[0].split('/')[3] + '/' + \
                        line.strip().split(' ')[0].split('/')[4] + '/zinc_drug_brick/' + folder + 'whole-brick/' + brick_name + ' ' + \
                        line.strip().split(' ')[1] + ' ' + line.strip().split(' ')[2] + ' ' + line.strip().split(' ')[3] + ' ' + \
                        line.strip().split(' ')[4] + ' ' + line.strip().split(' ')[5] + ' ' + line.strip().split(' ')[6] + ' ' + \
                        line.strip().split(' ')[7] + ' ' + line.strip().split(' ')[8] + '\n'
                fw.write(new_line)
    print(brick + ' is extracted.')
print('Extraction done!')
