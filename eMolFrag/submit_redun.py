import os
from multiprocessing import Pool, cpu_count
import time

def submit_redun(outputDir, outputFolderPath_log, tcBorder, pool):
    os.system('python /public/home/hqy/eMolFrag/eMolFrag/eMolFrag_2017_06_19_01/rm_block_redun.py ' + 
            outputDir + ' ' + outputFolderPath_log + ' ' + tcBorder + ' ' + pool)

start = time.time()
p = Pool(400)
folders = os.listdir('zinc_drug_like_new_out')
for folder in folders:
    outputDir = 'zinc_drug_bricks/' + folder + '/'
    outputFolderPath_log = 'zinc_drug_bricks/' + folder + 'output-log/'
    tcBorder = '1.0'
    pool = '10'
    p.apply_async(submit_redun, args=(outputDir, outputFolderPath_log, tcBorder, pool, ))
print('waiting for all processes')
p.close()
p.join()
end = time.time()
print("Total time {} s".format((end - start)))
print('=' * 100)
print('remove bricks redundancy done!')
