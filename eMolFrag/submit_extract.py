import os
from multiprocessing import Pool, cpu_count
import time

def submit_extract(folder, outputfolder, outputlog):
    os.system('python /public/home/hqy/eMolFrag/extract_brick.py ' + 
            folder + ' ' + outputfolder + ' ' + outputlog)

start = time.time()
p = Pool(40)
folders = os.listdir('zinc_drug_like_new_out')
for folder in folders:
    folder = folder + '/'
    outputfolder = 'zinc_drug_bricks/' + folder + 'whole-brick/'
    outputlog = 'zinc_drug_bricks/' + folder + 'output-log/'
    p.apply_async(submit_extract, args=(folder, outputfolder, outputlog, ))
print('waiting for all processes')
p.close()
p.join()
end = time.time()
print("Total time {} s".format((end - start)))
print('=' * 100)
print('extraction done!')
