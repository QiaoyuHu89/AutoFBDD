import os
from multiprocessing import Pool, cpu_count
import time

def submit_eMolFrag(inputfolder, outputfolder):
    os.system('python /public/home/hqy/eMolFrag/eMolFrag/eMolFrag_2017_06_19_01/eMolFrag_edit.py -i ' + 
            inputfolder + ' -o ' + outputfolder + ' -p 10 -m 0 -c 0')

start = time.time()
p = Pool(400)
folders = os.listdir('zinc_drug_like')
for folder in folders:
    sub_folders = os.listdir('zinc_drug_like/' + folder)
    for sub_folder in sub_folders:
        input_folder = 'zinc_drug_like/' + folder + '/' + sub_folder
        if not os.path.exists('zinc_drug_like_new_out/' + folder):
            os.mkdir('zinc_drug_like_new_out/' + folder)
        output_folder = 'zinc_drug_like_new_out/' + folder + '/' + sub_folder + '_out'
        p.apply_async(submit_eMolFrag, args=(input_folder, output_folder, ))
print('waiting for all processes')
p.close()
p.join()
end = time.time()
print("Total time {} s".format((end - start)))
print('=' * 100)
print('eMolFrag done!')
