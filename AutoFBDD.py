import os
import sys
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")
import time
import glob
import shutil
import argparse
from pymol import cmd
from multiprocessing import Pool
from linking_postprocess import split_mol2

def main_linking(brickmode, folder, pdbfile, centerfile, brickfolder, bricklist, n, sep_num, t, k, val, num_pose):
    log = open(folder + '/AutoFBDD_log.txt', 'w', buffering=1)
    log.write('INFO: AutoFBDD(linking) has started.\n')
    cwd = os.getcwd()
    if brickmode == 'mol':
        # eMolFrag
        log.write('INFO: brick generation by eMolFrag has started.\n')
        os.chdir(folder)
        if not os.path.exists('eMolFrag'):
            os.mkdir('eMolFrag')
        os.system('cp -r ' + brickfolder + ' eMolFrag/')
        os.chdir('eMolFrag')
        os.system('python ' + AutoFBDD_FOL + '/eMolFrag/process.py -f ' + brickfolder)
        os.system('python ' + AutoFBDD_FOL + '/eMolFrag/eMolFrag/eMolFrag_2017_06_19_01/eMolFrag_edit.py -i ' + \
            brickfolder + '_processed' + ' -o ' + brickfolder + '_out -p ' + str(n) + ' -m 0 -c 0')
        # ifitdock
        log.write('INFO: brick ifitdock has started.\n')
        os.chdir('../')
        if not os.path.exists('ifit_dock'):
            os.mkdir('ifit_dock')
        os.system('cp ../ifitdock/ifit_dock/* ifit_dock/')
        os.system('cp -r eMolFrag/' + brickfolder + '_out/output-brick ifit_dock/')
        os.system('cp ' + pdbfile + ' ' + centerfile + ' ifit_dock/')
        os.chdir('ifit_dock')
        os.system('python ifitdock.py -i ' + pdbfile + ' -c ' + centerfile + ' -f output-brick -l ' + \
            bricklist + ' -n ' + str(n) + ' -s ' + str(sep_num) + ' --parallel')

    elif brickmode == 'brick':
        # ifitdock
        log.write('INFO: brick ifitdock has started.\n')
        os.chdir(folder)
        if not os.path.exists('ifit_dock'):
            os.mkdir('ifit_dock')
        os.system('cp ../ifitdock/ifit_dock/* ifit_dock/')
        if brickfolder not in ['zinc', 'chembl', 'pdbbind', 'all', 'all_0.4', 'all_adme']:
            os.system('cp -r ' + pdbfile + ' ' + centerfile + ' ' + brickfolder + ' ifit_dock/')
        else:
            os.system('cp ' + pdbfile + ' ' + centerfile + ' ifit_dock/')
            #os.system('cp -r ' + AutoFBDD_FOL + '/brick_library/' + brickfolder + ' ifit_dock/')
        os.chdir('ifit_dock')
        os.system('python ifitdock.py -i ' + pdbfile + ' -c ' + centerfile + ' -f ' + brickfolder + ' -l ' + \
                bricklist + ' -n ' + str(n) + ' -s ' + str(sep_num) + ' --parallel')

    else:
        log.write('INFO: Unrecognized brickmode! Please select: mol or brick.\n')
        print('Unrecognized brickmode! Please select: mol or brick.')
        log.close()
        sys.exit()
    
    # ifitdock postprocess
    log.write('INFO: ifitdock postprocess has started.\n')
    os.system('python ifitdock_postprocess.py -f cluster_out_cluster_info.list -t ' + str(t) + ' -k ' + str(k) + ' -v ' + str(val))
    if not os.path.exists('brickfolder_ori'):
        os.mkdir('brickfolder_ori')
    bricks = glob.glob('cluster_out*sorted*brick*')
    if bricks != []:
        for brick in bricks:
            shutil.copy(brick, 'brickfolder_ori')
    else:
        log.write('INFO: There are no bricks generated through ifitdock clustering.\n')
        log.close()
        sys.exit()
    
    os.chdir('../')
    if not os.path.exists('brick_linking'):
        os.mkdir('brick_linking')
    shutil.copy(pdbfile, 'brick_linking/')
    shutil.copy('ifit_dock/app_bricks_1.list', 'brick_linking/')
    shutil.copytree('ifit_dock/brickfolder_ori', 'brick_linking/brickfolder_ori')
    shutil.copytree('ifit_dock/brickfolder_ori', 'brick_linking/brickfolder_1')
    
    os.chdir('brick_linking')
    iter_idx = 1
    while os.path.isfile('app_bricks_' + str(iter_idx) + '.list'):
        # prepare data for brick linking
        log.write('INFO: Data prepare for brick linking in iteration ' + str(iter_idx) + '.\n')
        shutil.copy(AutoFBDD_FOL + '/DEVELOP/examples/data_pre_for_linking_pharms.py', './')
        os.system('python data_pre_for_linking_pharms.py -l app_bricks_' + str(iter_idx) + '.list -f brickfolder_' + str(iter_idx) \
                    + ' -t ' + pdbfile + ' -i ' + str(iter_idx))
        
        # perform brick linking
        log.write('INFO: brick linking in iteration ' + str(iter_idx) + '.\n')
        shutil.copy(AutoFBDD_FOL + '/DEVELOP/examples/brick_linking_pharms.py', './')
        shutil.copy(AutoFBDD_FOL + '/DEVELOP/models/linker_design/pretrained_DEVELOP_model_pharms.pickle', './')
        os.system('python brick_linking_pharms.py -l app_bricks_' + str(iter_idx) + '.list -d data_path_' + str(iter_idx) \
                    + ' -i ' + str(iter_idx))
        
        # evaluate generated molecules
        log.write('INFO: Generated mols evaluation in iteration ' + str(iter_idx) + '.\n')
        os.system('cp ' + AutoFBDD_FOL + '/DEVELOP/analysis/* ./')
        os.system('python assess_linking_pharms.py -d data_path_' + str(iter_idx) + ' -a ZINC -f generated_smi_' + str(iter_idx) + \
                    ' -t ' + AutoFBDD_FOL + '/DEVELOP/data/linker_design/data_zinc_train.txt -s ./ -n ' + str(n) + \
                    ' -v True -r None -p ./wehi_pains.csv -i ' + str(iter_idx))
        
        # brick linking postprocess
        log.write('INFO: brick linking postprocess in iteration ' + str(iter_idx) + '.\n')
        os.system('cp -r ' + AutoFBDD_FOL + '/ifitdock/ifit_eval ./')
        os.system('cp ' + AutoFBDD_FOL + '/' + folder + '/' + pdbfile + ' ' + AutoFBDD_FOL + '/' + folder + '/' + centerfile + ' ./ifit_eval')
        os.system('python linking_postprocess.py -i ' + pdbfile + ' -c ' + centerfile + ' -e eval_folder -n ' + str(n) + \
                    ' -p ' + str(num_pose) + ' -t ' + str(iter_idx))

        # select linked bricks
        log.write('INFO: select linked bricks in iteration ' + str(iter_idx) + '.\n')
        os.system('python select_linked_bricks.py -i ' + str(iter_idx) + ' -v ' + str(val))

        # updata iter_idx
        iter_idx += 1

    # generate delta_G for directly connected molecules
    log.write('INFO: generate delta_G of directly connected molecules.\n')
    screeningC_list = glob.glob('brickfolder_*/*screeningC.mol2')
    if screeningC_list != []:
        if not os.path.exists('screeningC_eval'):
            os.mkdir('screeningC_eval')
        if not os.path.exists('screeningC_eval/eval_folder'):
            os.mkdir('screeningC_eval/eval_folder')
        if not os.path.exists('brickfolder_C'):
            os.mkdir('brickfolder_C')
        os.system('cp ifit_eval/* screeningC_eval')
        
        for mol2 in screeningC_list:
            os.system('obabel ' + mol2 + ' -OscreeningC_eval/eval_folder/' + os.path.basename(mol2).split('.')[0] + '.sdf')
        os.chdir('screeningC_eval')
        os.system('python ifitdock_eval.py -i ' + pdbfile + ' -c ' + centerfile + ' -f eval_folder')
        screening_sorted_file = pdbfile.split('.')[0] + '_screening_sorted.mol2'
        if os.path.exists(screening_sorted_file) and os.path.getsize(screening_sorted_file) != 0:
            split_mol2(screening_sorted_file, num_pose)
        os.system('cp ' + screening_sorted_file.split('.')[0] + '_*mol2 ../brickfolder_C')

    # select final generated molecules according to delta_G
    log.write('INFO: select final generated molecules according to delta_G.\n')
    if screeningC_list != []:
        os.chdir('../../')
    else:
        os.chdir('../')
    if not os.path.exists('final_results'):
        os.mkdir('final_results')
    shutil.copy(pdbfile, 'final_results/')
    screening_mols = glob.glob('brick_linking/brickfolder_*/*screening*mol2')
    screeningC_mols = glob.glob('brick_linking/brickfolder_*/*screeningC.mol2')
    new_screening_mols = [mol for mol in screening_mols if mol not in screeningC_mols]
    screening_dict = {}
    for screening_mol in new_screening_mols:
        with open(screening_mol, 'r') as fr:
            lines = fr.readlines()
        for line in lines:
            if line.startswith('delta_G'):
                delta_G = float(line.split(' ')[2])
        screening_dict[screening_mol] = delta_G
    screening_list = sorted(screening_dict.items(), key=lambda x: x[1], reverse=False)
    
    if len(screening_list) > 100:
        screening_list = screening_list[0:100]
    else:
        pass
    for screening_mol in screening_list:
        final_mol = 'final_results/' + os.path.splitext(pdbfile)[0] + '_' + str(screening_list.index(screening_mol)) + '.mol2'
        shutil.copy(screening_mol[0], final_mol)
        with open('screening_list.txt', 'a') as fw:
            fw.write(screening_mol[0] + ' ' + final_mol + '\n')
        # os.system('babel ' + final_mol + ' ' + os.path.splitext(final_mol)[0] + '.pdb')
        # os.system('babel ' + os.path.splitext(final_mol)[0] + '.pdb ' + final_mol)
        # os.system('rm ' + os.path.splitext(final_mol)[0] + '.pdb')
        cmd.load(final_mol)
        cmd.remove('hydrogens')
        cmd.save(final_mol)
        cmd.delete('all')
    for screening_mol in screening_list:
        final_mol = 'final_results/' + os.path.splitext(pdbfile)[0] + '_' + str(screening_list.index(screening_mol)) + '.mol2'
        cmd.load(final_mol)
    cmd.load(pdbfile)
    cmd.remove('hydrogens')
    cmd.save(os.path.splitext(pdbfile)[0] + '_results.pse')
    cmd.delete('all')
    shutil.move(os.path.splitext(pdbfile)[0] + '_results.pse', 'final_results/')
    os.chdir(AutoFBDD_FOL)
    print('AutoFBDD done!')
    log.write('INFO: AutoFBDD has finished.\n')
    log.close()

def main_growing(brickmode, folder, pdbfile, centerfile, brickfolder, bricklist, n, sep_num, t, k, val, num_pose):
    log = open(folder + '/AutoFBDD_log.txt', 'w', buffering=1)
    log.write('INFO: AutoFBDD(growing) has started.\n')
    cwd = os.getcwd()
    if brickmode == 'mol':
        # eMolFrag
        log.write('INFO: brick generation by eMolFrag has started.\n')
        os.chdir(folder)
        if not os.path.exists('eMolFrag'):
            os.mkdir('eMolFrag')
        os.system('cp -r ' + brickfolder + ' eMolFrag/')
        os.chdir('eMolFrag')
        os.system('python ' + AutoFBDD_FOL + '/eMolFrag/process.py -f ' + brickfolder)
        os.system('python ' + AutoFBDD_FOL + '/eMolFrag/eMolFrag/eMolFrag_2017_06_19_01/eMolFrag_edit.py -i ' + \
            brickfolder + '_processed' + ' -o ' + brickfolder + '_out -p ' + str(n) + ' -m 0 -c 0')
        # ifitdock
        log.write('INFO: brick ifitdock has started.\n')
        os.chdir('../')
        if not os.path.exists('ifit_dock'):
            os.mkdir('ifit_dock')
        os.system('cp ../ifitdock/ifit_dock/* ifit_dock/')
        os.system('cp -r eMolFrag/' + brickfolder + '_out/output-brick ifit_dock/')
        os.system('cp ' + pdbfile + ' ' + centerfile + ' ifit_dock/')
        os.chdir('ifit_dock')
        os.system('python ifitdock.py -i ' + pdbfile + ' -c ' + centerfile + ' -f output-brick -l ' + \
            bricklist + ' -n ' + str(n) + ' -s ' + str(sep_num) + ' --parallel')

    elif brickmode == 'brick':
        # ifitdock
        log.write('INFO: brick ifitdock has started.\n')
        os.chdir(folder)
        if not os.path.exists('ifit_dock'):
            os.mkdir('ifit_dock')
        os.system('cp ../ifitdock/ifit_dock/* ifit_dock/')
        if brickfolder not in ['zinc', 'chembl', 'pdbbind', 'all', 'all_0.4', 'all_adme']:
            os.system('cp -r ' + pdbfile + ' ' + centerfile + ' ' + brickfolder + ' ifit_dock/')
        else:
            os.system('cp ' + pdbfile + ' ' + centerfile + ' ifit_dock/')
            #os.system('cp -r ' + AutoFBDD_FOL + '/brick_library/' + brickfolder + ' ifit_dock/')
        os.chdir('ifit_dock')
        os.system('python ifitdock.py -i ' + pdbfile + ' -c ' + centerfile + ' -f ' + brickfolder + ' -l ' + \
                bricklist + ' -n ' + str(n) + ' -s ' + str(sep_num) + ' --parallel')

    else:
        log.write('INFO: Unrecognized brickmode! Please select: mol or brick.\n')
        print('Unrecognized brickmode! Please select: mol or brick')
        log.close()
        sys.exit()

    # ifitdock postprocess
    log.write('INFO: ifitdock postprocess has started.\n')
    os.system('python ifitdock_postprocess.py -f cluster_out_cluster_info.list -t ' + str(t) + ' -k ' + str(k) + ' -v ' + str(val))
    if not os.path.exists('brickfolder_ori'):
        os.mkdir('brickfolder_ori')
    bricks = glob.glob('cluster_out*sorted*brick*')
    if bricks != []:
        for brick in bricks:
            shutil.copy(brick, 'brickfolder_ori')
    else:
        log.write('INFO: There are no bricks generated through ifitdock clustering.\n')
        log.close()
        sys.exit()
    
    os.chdir('../')
    if not os.path.exists('brick_growing'):
        os.mkdir('brick_growing')
    shutil.copy(pdbfile, 'brick_growing/')
    shutil.copy('ifit_dock/app_bricks_1.list', 'brick_growing/')
    shutil.copytree('ifit_dock/brickfolder_ori', 'brick_growing/brickfolder_ori')
    shutil.copytree('ifit_dock/brickfolder_ori', 'brick_growing/brickfolder_1')
    
    os.chdir('brick_growing')
    iter_idx = 1
    while os.path.isfile('app_bricks_' + str(iter_idx) + '.list'):
        os.system('rm app_bricks_' + str(iter_idx) + '.list')
        # prepare data for brick growing
        log.write('INFO: Data prepare for brick growing in iteration ' + str(iter_idx) + '.\n')
        shutil.copy(AutoFBDD_FOL + '/DEVELOP/examples/data_pre_for_growing.py', './')
        os.system('python data_pre_for_growing.py -f brickfolder_' + str(iter_idx) + ' -t ' + pdbfile \
                    + ' -i ' + str(iter_idx))
        
        # perform brick growing
        log.write('INFO: brick growing in iteration ' + str(iter_idx) + '.\n')
        shutil.copy(AutoFBDD_FOL + '/DEVELOP/examples/brick_growing.py', './')
        shutil.copy(AutoFBDD_FOL + '/DEVELOP/models/scaffold_elaboration/pretrained_DEVELOP_model.pickle', './')
        os.system('python brick_growing.py -l app_bricks_' + str(iter_idx) + '.list -d data_path_' + str(iter_idx) \
                    + ' -i ' + str(iter_idx))
        
        # evaluate generated molecules
        log.write('INFO: Generated mols evaluation in iteration ' + str(iter_idx) + '.\n')
        os.system('cp ' + AutoFBDD_FOL + '/DEVELOP/analysis/* ./')
        os.system('python assess_growing_mols.py -d data_path_' + str(iter_idx) + ' -a ZINC -f generated_smi_' + str(iter_idx) + \
                    ' -t ' + AutoFBDD_FOL + '/DEVELOP/data/scaffold_elaboration/data_zinc_train.txt -s ./ -n ' + str(n) + \
                    ' -v True -r None -p ./wehi_pains.csv -i ' + str(iter_idx))
        
        # brick growing postprocess
        log.write('INFO: brick growing postprocess in iteration ' + str(iter_idx) + '.\n')
        os.system('cp -r ' + AutoFBDD_FOL + '/ifitdock/ifit_eval ./')
        os.system('cp ' + AutoFBDD_FOL + '/' + folder + '/' + pdbfile + ' ' + AutoFBDD_FOL + '/' + folder + '/' + centerfile + ' ./ifit_eval')
        os.system('python growing_postprocess.py -i ' + pdbfile + ' -c ' + centerfile + ' -e eval_folder -n ' + str(n) + \
                    ' -p ' + str(num_pose) + ' -t ' + str(iter_idx))

        # select grown bricks
        log.write('INFO: select grown bricks in iteration ' + str(iter_idx) + '.\n')
        if os.path.exists('brickfolder_' + str(iter_idx+1)):
            bricks_list = glob.glob('brickfolder_' + str(iter_idx+1) + '/*screening*')
            if bricks_list != []:
                with open('app_bricks_' + str(iter_idx+1) + '.list', 'w') as fw:
                    for bricks in bricks_list:
                        fw.write(str(bricks.split('/')[1]) + '\n')
            print(len(bricks_list))

        # updata iter_idx
        iter_idx += 1

    # select final generated molecules according to delta_G
    log.write('INFO: select final generated molecules according to delta_G.\n')
    os.chdir('../')
    if not os.path.exists('final_results'):
        os.mkdir('final_results')
    shutil.copy(pdbfile, 'final_results/')
    screening_mols = glob.glob('brick_growing/brickfolder_*/*screening*mol2')
    screening_dict = {}
    for screening_mol in screening_mols:
        with open(screening_mol, 'r') as fr:
            lines = fr.readlines()
        for line in lines:
            if line.startswith('delta_G'):
                delta_G = float(line.split(' ')[2])
        screening_dict[screening_mol] = delta_G
    screening_list = sorted(screening_dict.items(), key=lambda x: x[1], reverse=False)
    
    if len(screening_list) > 100:
        screening_list = screening_list[0:100]
    else:
        pass
    for screening_mol in screening_list:
        final_mol = 'final_results/' + os.path.splitext(pdbfile)[0] + '_' + str(screening_list.index(screening_mol)) + '.mol2'
        shutil.copy(screening_mol[0], final_mol)
        # os.system('babel ' + final_mol + ' ' + os.path.splitext(final_mol)[0] + '.pdb')
        # os.system('babel ' + os.path.splitext(final_mol)[0] + '.pdb ' + final_mol)
        # os.system('rm ' + os.path.splitext(final_mol)[0] + '.pdb')
        cmd.load(final_mol)
        cmd.remove('hydrogens')
        cmd.save(final_mol)
        cmd.delete('all')
    for screening_mol in screening_list:
        final_mol = 'final_results/' + os.path.splitext(pdbfile)[0] + '_' + str(screening_list.index(screening_mol)) + '.mol2'
        cmd.load(final_mol)
    cmd.load(pdbfile)
    cmd.remove('hydrogens')
    cmd.save(os.path.splitext(pdbfile)[0] + '_results.pse')
    cmd.delete('all')
    shutil.move(os.path.splitext(pdbfile)[0] + '_results.pse', 'final_results/')
    os.chdir(AutoFBDD_FOL)
    print('AutoFBDD done!')
    log.write('INFO: AutoFBDD has finished.\n')
    log.close()


if __name__=="__main__":
    # Usage: python AutoFBDD.py --mode mode --brick_mode brickmode --folder folder --input_pdb target.pdb --center center.txt \
    # --brickfolder brickfolder --brickfile_list bricks_file.txt --num_cpu 20 --sep_bricks 5 --top_clusters 10 --top_bricks 3 --dis_val 10 --poses 3
    # Get the arguments for ifitdock
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('--mode', type=str, help='FBDD mode: linking or growing')
    parser.add_argument('--brick_mode', type=str, help='brick mode: brick or mol. If you want to use brick_library or provide your own brickfolder, '
    'please specify brick in brick_mode. If you have a list of small molecule drugs and want to cut them into bricks, please specify mol in brick_mode.')
    parser.add_argument('--folder', type=str, help='main folder')
    parser.add_argument('--input_pdb', type=str, help='input pdb file')
    parser.add_argument('--center', type=str, help='pocket center file')
    parser.add_argument('--brickfolder', type=str, help='brick folder. If you want to use brick_library, '
    'please specify zinc or chembl or pdbbind or all. Otherwise, provide your own brick or molecule folder and put it in the main folder.')
    parser.add_argument('--brickfile_list', type=str, help='brick file list, just give an empty txt file that can be used to store bricks')
    parser.add_argument('--num_cpu', type=int, default=20, help='the number of cpu used in parallel run')
    parser.add_argument('--sep_bricks', type=int, default=5, help='the number of bricks in each brick list file')
    # Get the arguments for ifitdock postprocess
    parser.add_argument('--top_clusters', type=int, default=10, help='the number of top clusters saved')
    parser.add_argument('--top_bricks', type=int, default=3, help='the number of top bricks saved in each top cluster')
    parser.add_argument('--dis_val', type=float, default=10.0, help='threshold distance value of two bricks for brick linking')
    # Get the arguments for appropriate molecules selection
    parser.add_argument('--poses', type=int, default=3, help='the number of poses saved for each pair of linked or grown bricks')
    args = parser.parse_args()
    
    mode = args.mode
    brickmode = args.brick_mode
    folder = args.folder
    pdbfile = args.input_pdb
    centerfile = args.center
    brickfolder = args.brickfolder
    bricklist = args.brickfile_list
    n = args.num_cpu
    sep_num = args.sep_bricks
    
    t = args.top_clusters
    k = args.top_bricks
    val = args.dis_val
    
    num_pose = args.poses
    
    if mode == 'linking':
        main_linking(brickmode, folder, pdbfile, centerfile, brickfolder, bricklist, n, sep_num, t, k, val, num_pose)
    elif mode == 'growing':
        main_growing(brickmode, folder, pdbfile, centerfile, brickfolder, bricklist, n, sep_num, t, k, val, num_pose)
    else:
        print('Unregonized molecules generation mode! Please select: linking or growing')
