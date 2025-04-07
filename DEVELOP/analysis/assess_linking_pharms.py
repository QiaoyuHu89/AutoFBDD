import os
import re
import sys
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")

import glob
import time
import shutil
import argparse
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

# delete the dummy atoms of bricks
def delete_du(brick):
    if cmd._COb is None:
        import pymol2
        import pymol.invocation
        pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
        pymol2.SingletonPyMOL().start()
        
    cmd.load(brick)
    cmd.select('organic_atoms', 'name H+B+C+N+O+S+P+F+Cl+Br+I')
    cmd.select('du', 'not organic_atoms')
    cmd.remove('du')
    cmd.save(brick)
    cmd.delete('all')
    mol = Chem.SDMolSupplier(brick)[0]
    if mol is not None:
        print(brick + ' du is deleted!')

# evaluate the generated molecules using 2D and 3D analysis according to RMSD
def assess_mols_RMSD(dataset, gen_smi, frag_sdf, train_set_path, save_path, output_name, \
                num_cores, verbose, restrict, pains_smarts_loc, PROTACs):
    os.system('python evaluate_linking_mols_edit.py ' + dataset + ' ' + gen_smi + ' ' + frag_sdf + ' ' + \
            train_set_path + ' ' + save_path + ' ' + output_name + ' ' + num_cores + ' ' + verbose + ' ' + \
            restrict + ' ' + pains_smarts_loc + ' ' + PROTACs + ' >> ' + output_name + '_RMSD_log.txt')
    
# evaluate the generated molecules using 2D and 3D analysis according to linker pharmacophore
def assess_mols_pharm(dataset, gen_smi, frag_sdf, train_set_path, save_path, output_name, \
                num_cores, verbose, restrict, pains_smarts_loc, data_path, PROTACs):
    os.system('python evaluate_linking_mols_pharms.py ' + dataset + ' ' + gen_smi + ' ' + frag_sdf + ' ' + \
            train_set_path + ' ' + save_path + ' ' + output_name + ' ' + num_cores + ' ' + verbose + ' ' + \
            restrict + ' ' + pains_smarts_loc + ' ' + data_path + ' ' + PROTACs + ' >> ' + output_name + '_pharm_log.txt')

# extract molecules that align well with original fragments
def extract_sdf(log_file, top_sdf_folder, top_smi_folder, output_name):
    with open(log_file, 'r') as fr:
        lines = fr.readlines()
    for line in lines:
        if line.startswith('Number of valid SMILES: '):
            num_valid_smi = line.split()[4]
        if line.startswith('best_idx_SC_RDKit_Full: '):
            best_idx_SC_RDKit_Full = line.split('  ')[1]
        elif  line.startswith('best_idx_SC_RDKit_Frag: '):
            best_idx_SC_RDKit_Frag = line.split('  ')[1]
        elif  line.startswith('best_idx_RMSD_Frag: '):
            best_idx_RMSD_Frag = line.split('  ')[1]
        elif line.startswith('best_mols_RMSD_Frag: '):
            best_mols_RMSD_Frag = line.split('  ')[1]

    if num_valid_smi != '0':
        idx_list_1 = best_idx_SC_RDKit_Full.strip().replace('[','').replace(']','').split(', ')
        print(idx_list_1)
        idx_list_2 = best_idx_SC_RDKit_Frag.strip().replace('[','').replace(']','').split(', ')
        print(idx_list_2)
        idx_list_3 = best_idx_RMSD_Frag.strip().replace('[','').replace(']','').split(', ')
        print(idx_list_3)

        if not os.path.exists(top_sdf_folder):
            os.mkdir(top_sdf_folder)
        gen_sdfs = Chem.SDMolSupplier(output_name + '/generated_mols.sdf')
        writer_1 = Chem.SDWriter(top_sdf_folder + '/best_mols_SC_RDKit_Full.sdf')
        for i in idx_list_1:
            gen_mol = gen_sdfs[int(i)]
            writer_1.write(gen_mol)
        writer_1.close()
        print('best_mols_SC_RDKit_Full.sdf is created successfully')

        writer_2 = Chem.SDWriter(top_sdf_folder + '/best_mols_SC_RDKit_Frag.sdf')
        for i in idx_list_2:
            gen_mol = gen_sdfs[int(i)]
            writer_2.write(gen_mol)
        writer_2.close()
        print('best_mols_SC_RDKit_Frag.sdf is created successfully')

        writer_3 = Chem.SDWriter(top_sdf_folder + '/best_mols_RMSD_Frag.sdf')
        for i in idx_list_3:
            gen_mol = gen_sdfs[int(i)]
            writer_3.write(gen_mol)
        writer_3.close()
        print('best_mols_RMSD_Frag.sdf is created successfully')
        
        smis_list = []
        if not os.path.exists(top_smi_folder):
            os.mkdir(top_smi_folder)

        best_smis_nums = []
        for i in range(len(re.findall('\(.+?\)\),', best_mols_RMSD_Frag))):
            best_smis_nums.append(re.findall('\(.+?\)\),', best_mols_RMSD_Frag)[i][0:-1])
        for j in best_smis_nums:
            sdf = top_sdf_folder + '/' + str(j.split(',')[2][2:]) + '.sdf'
            smi = top_smi_folder + '/' + str(j.split(',')[2][2:]) + '.smi'
            gen_mol = gen_sdfs[int(j.split(',')[2][2:])]
            writer = Chem.SDWriter(sdf)
            writer.write(gen_mol)
            writer.close()
            with open(smi, 'w') as fw:
                if j.split(',')[0][2:-1] not in smis_list:
                    fw.write(j.split(',')[0][2:-1])
                    smis_list.append(j.split(',')[0][2:-1])
        print(smis_list)
    else:
        print('No molecules pass 2D filters.')

# align generated molecules to original fragments
def align_mols(best_RMSD_mols, frag_sdf, output_sdf):
    gen_sdfs = Chem.SDMolSupplier(best_RMSD_mols)
    ref_mol = Chem.SDMolSupplier(frag_sdf)[0]
    writer = Chem.SDWriter(output_sdf)
    gen_smis = []
    for i in range(len(gen_sdfs)):
        gen_mol = gen_sdfs[int(i)]
        score = sascorer.calculateScore(gen_mol)
        gen_smi = Chem.MolToSmiles(gen_mol)
        if gen_mol is not None and score <= 5.0 and gen_smi not in gen_smis:
            gen_smis.append(gen_smi)
            pyO3A = rdMolAlign.GetO3A(gen_mol, ref_mol).Align()
            try:
                writer.write(gen_mol)
            except:
                continue
            print('write one aligned molecule')
        else:
            pass
    writer.close()
    print('aligned_RMSD_mols.sdf is created successfully!')

def main(data_path, dataset, gen_smi_folder, train_set_path, save_path, num_cores, verbose, restrict, pains_smarts_loc, PROTACs, iter_idx):
    if not os.path.exists('eval_results_' + str(iter_idx)):
        os.mkdir('eval_results_' + str(iter_idx))
    if not os.path.exists('eval_results_' + str(iter_idx) + '/success'):
        os.mkdir('eval_results_' + str(iter_idx) + '/success')
    if not os.path.exists('eval_results_' + str(iter_idx) + '/fail'):
        os.mkdir('eval_results_' + str(iter_idx) + '/fail')
    
    frag_list = [frag for frag in os.listdir(data_path) if frag.endswith('frag.sdf')]
    smi_list = os.listdir(gen_smi_folder)
    smi_list.sort()
    for smi in smi_list:
        name = os.path.splitext(smi)[0].split('_generated_')[0]
        if not os.path.exists(name):
            os.mkdir(name)
        shutil.copy(gen_smi_folder + '/' + smi, name)
        for frag in frag_list:
            if os.path.splitext(frag)[0].split('_frag')[0] == name:
                shutil.copy(data_path + '/' + frag, name)
        gen_smi = name + '/' + name + '_generated_smiles.smi'
        frag_sdf = name + '/' + name + '_frag.sdf'
        delete_du(frag_sdf)
        try:
            assess_mols_RMSD('ZINC', gen_smi, frag_sdf, train_set_path, save_path, name, \
                    str(num_cores), verbose, restrict, pains_smarts_loc, PROTACs)
            
            if not os.path.exists(name + '/RMSD'):
                os.mkdir(name + '/RMSD')
            shutil.move(name + '_RMSD_log.txt', name + '/RMSD')
            log_file = name + '/RMSD/' + name + '_RMSD_log.txt'

            extract_sdf(log_file, 'gen_sdfs', 'gen_smis', name)
            if os.path.isdir('gen_sdfs') and os.path.isdir('gen_smis'):
                shutil.move('gen_sdfs', name + '/RMSD')
                shutil.move('gen_smis', name + '/RMSD')
            
            shutil.move(name + '/generated_mols.sdf', name + '/RMSD')
            best_RMSD_mols = name + '/RMSD/gen_sdfs/best_mols_RMSD_Frag.sdf'
            output_sdf = name + '/RMSD/aligned_RMSD_mols.sdf'
            align_mols(best_RMSD_mols, frag_sdf, output_sdf)
        except:
            print(name + ' did not generate valid RMSD sdf.')
            shutil.move(name, 'eval_results_' + str(iter_idx) + '/fail')
            continue
        
        try:
            assess_mols_pharm(dataset, gen_smi, frag_sdf, train_set_path, save_path, name, \
                    str(num_cores), verbose, restrict, pains_smarts_loc, data_path, PROTACs)
            
            if not os.path.exists(name + '/pharm'):
                os.mkdir(name + '/pharm')
            shutil.move(name + '_pharm_log.txt', name + '/pharm')
            log_file = name + '/pharm/' + name + '_pharm_log.txt'

            extract_sdf(log_file, 'gen_sdfs', 'gen_smis', name)
            if os.path.isdir('gen_sdfs') and os.path.isdir('gen_smis'):
                shutil.move('gen_sdfs', name + '/pharm')
                shutil.move('gen_smis', name + '/pharm')
            
            shutil.move(name + '/generated_mols.sdf', name + '/pharm')
            best_RMSD_mols = name + '/pharm/gen_sdfs/best_mols_RMSD_Frag.sdf'
            output_sdf = name + '/pharm/aligned_pharm_mols.sdf'
            align_mols(best_RMSD_mols, frag_sdf, output_sdf)
        except:
            print(name + ' did not generate valid pharm sdf.')
        
        shutil.move(name, 'eval_results_' + str(iter_idx) + '/success')
        print(name + ' has been evaluated!')

if __name__=="__main__":
    # usage: python assess_linking_mols.py -d data_path -a ZINC -f generated_smi -t ../data/linker_design/data_zinc_train.txt 
    # -s ./ -n 10 -v True -r None -p ./wehi_pains.csv -i iter_idx
    # Get the arguments
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-d', type=str, help='data path folder')
    parser.add_argument('-a', type=str, help='data set: ZINC or CASF')
    parser.add_argument('-f', type=str, help='generated smiles folder')
    parser.add_argument('-t', type=str, help='train set path')
    parser.add_argument('-s', type=str, help='save path')
    parser.add_argument('-n', type=int, help='number of cores')
    parser.add_argument('-v', type=str, help='verbose')
    parser.add_argument('-r', type=str, help='restrict')
    parser.add_argument('-p', type=str, help='pains smarts file location')
    parser.add_argument('--protacs', default='False', action='store_true', help='assess PROTACs or not')
    parser.add_argument('-i', type=int, help='iteration index')
    
    args = parser.parse_args()

    data_path = args.d
    dataset = args.a
    gen_smi_folder = args.f
    train_set_path = args.t
    save_path = args.s
    num_cores = args.n
    verbose = args.v
    restrict = args.r
    pains_smarts_loc = args.p
    PROTACs = args.protacs
    iter_idx = args.i
    
    start = time.time()
    main(data_path, dataset, gen_smi_folder, train_set_path, save_path, num_cores, verbose, restrict, pains_smarts_loc, PROTACs, iter_idx)
    end = time.time()
    print('='*100)
    print('Running time: %s Seconds' % (end-start))
    print('Evaluating generated molecules done!')
