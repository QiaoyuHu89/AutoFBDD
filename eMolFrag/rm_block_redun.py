import subprocess
import os
#import path
import glob
import os.path
import shutil
import sys
import time
from multiprocessing import Pool
from functools import partial
from rdkit import Chem
from rmRedBrick01 import RmBrickRed

def PrintLog(path, msg):
    # write log
    with open(path, 'at') as outLog:
        outLog.write(time.asctime( time.localtime(time.time()) ))
        outLog.write(msg)
        outLog.write('\n')

def RmBrickRedundancy(outputDir, outputFolderPath_log, tcBorder, pool):
    # Brick Part
    #Step 3: Form and group lists by atom numbers
    fileNameAndAtomNumList_R=[]
    with open(outputFolderPath_log+'BrickListAll.txt','r') as inList:
        fileNameAndAtomNumList_R=inList.readlines()

    FNAANLList_R=[] #file name and atom number list list
    for FNAAN in fileNameAndAtomNumList_R: #FNAAN: file name and atom number
        FNAANList=FNAAN.split() #FNAANList: file name and atom number list
        FNAANLList_R.append([FNAANList[0],FNAANList[1:]])

    atomNumPro_R=[]
    for tempValue in FNAANLList_R: #tempValue: [[filename],['T','','C','','N','','O','']]
        if tempValue[1] not in atomNumPro_R: #tempValue[1]: ['T','','C','','N','','O',''],Group Property
            atomNumPro_R.append(tempValue[1])

    fileNameGroup_R=[[y[0] for y in FNAANLList_R if y[1]==x] for x in atomNumPro_R]

    with open(outputFolderPath_log+'BrickGroupList.txt','w') as groupOut:
        for i in range(len(atomNumPro_R)):
            groupOut.write(' '.join(atomNumPro_R[i])+' - ')
            groupOut.write('File Num: ')
            groupOut.write(str(len(fileNameGroup_R[i])))
            groupOut.write('\n')

    # Log
    path = outputFolderPath_log+'Process.log'
    msg = ' Start Remove Brick Redundancy '
    PrintLog(path, msg)

    #Step 4: Generate similarity data and etc.
    fileNameGroup_Rs=sorted(fileNameGroup_R,key=lambda x:len(x),reverse=True) #process long list first

    partial_RmBrick=partial(RmBrickRed, outputDir, tcBorder)
    pool.map(partial_RmBrick,fileNameGroup_Rs)

    # Log
    path = outputFolderPath_log+'Process.log'
    msg = ' End Remove Brick Redundancy '
    PrintLog(path, msg)
    
if __name__ == '__main__':
    args = sys.argv
    outputDir = args[1]
    outputFolderPath_log = args[2]
    tcBorder = float(args[3])
    processNum = int(args[4])
    pool=Pool(processes=processNum)
    RmBrickRedundancy(outputDir, outputFolderPath_log, tcBorder, pool)
    shutil.rmtree(outputDir + 'whole-brick/')