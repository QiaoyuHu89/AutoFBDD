import os
import glob
import argparse

def main(brickfolder):
    files = glob.glob(brickfolder + '/*')
    files.sort()
    if not os.path.exists(brickfolder + '_processed/'):
        os.mkdir(brickfolder + '_processed/')
    for file in files:
        with open(file, 'r') as fr:
            lines = fr.readlines()
        if '@<TRIPOS>SUBSTRUCTURE\n' in lines:
            idx = lines.index('@<TRIPOS>SUBSTRUCTURE\n')
            with open(brickfolder + '_processed/' + file.split('/')[1], 'w') as fw:
                for line in lines[0:idx]:
                   fw.write(line)
        else:
            os.system('cp ' + file + ' ' + brickfolder + '_processed/' + file.split('/')[1])
        print(file + ' is processed')
    print('=' * 100)
    print('All files are processed!')
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-f', type=str, help='brick folder')
    args = parser.parse_args()

    brickfolder = args.f
    main(brickfolder)
