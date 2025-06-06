import os
import glob
import time
import shutil
import argparse
from multiprocessing import Pool

def move(brick, outputfolder):
    shutil.move(brick, outputfolder)
    print(brick + ' is extracted.')

if __name__=="__main__":
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-f', type=str, help='brick folder')
    args = parser.parse_args()
    
    brickfolder = args.f
    inputfolder = 'zinc_drug_like_out/' + brickfolder + '/'
    outputfolder = '/home/bailab/other/hqy/zinc_drug_brick/whole_brick/' + brickfolder
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    start = time.time()
    p = Pool(20)
    brick_list = glob.glob(inputfolder + '*/output-chop-comb/b*.sdf')
    for brick in brick_list:
        p.apply_async(move, args=(brick, outputfolder, ))
    print('waiting for all processes')
    p.close()
    p.join()
    end = time.time()
    print("Total time {} s".format((end - start)))
    print('='*100)
    print('The total number of brick before extraction in ' + inputfolder + ' is: ' + str(len(brick_list)))
    brick_list = glob.glob(inputfolder + '*/output-chop-comb/b*.sdf')
    print('The total number of brick after extraction in ' + inputfolder + ' is: ' + str(len(brick_list)))
    print('Extraction done!')
