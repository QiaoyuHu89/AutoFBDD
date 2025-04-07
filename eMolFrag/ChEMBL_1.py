import os
import time

def submit_eMolFrag(inputfolder, outputfolder):
    os.system('python /home/bailab/data/hqy/PROTACs/eMolFrag/eMolFrag/eMolFrag_2017_06_19_01/eMolFrag_chembl.py -i ' + 
            inputfolder + ' -o ' + outputfolder + ' -p 8 -m 0 -c 0')

start = time.time()
folder = 'ChEMBL_part_1'
input_folder = 'ChEMBL/' + folder
output_folder = 'ChEMBL_out/' + folder + '_out'
submit_eMolFrag(input_folder, output_folder)
end = time.time()
print("Total time {} s".format((end - start)))
print('=' * 100)
print('eMolFrag done!')
