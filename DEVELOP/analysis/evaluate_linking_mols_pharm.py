#!/usr/bin/env python

# # Analysis script
# 
# Basic flow:
# - Load data
# - Check validity
# - Check 2D properties
# - Generate conformers
# - Check 3D similarity

# Imports

from unicodedata import name
import numpy as np
import re
import copy
import os, sys
AutoFBDD_FOL = os.environ['AutoFBDD_FOL']
sys.path.append(AutoFBDD_FOL + "/DEVELOP/")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/examples")
sys.path.append(AutoFBDD_FOL + "/DEVELOP/analysis/")
from collections import Counter

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem import rdMolAlign
from rdkit.Chem import Descriptors
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
import json

from joblib import Parallel, delayed

import rdkit_conf_parallel
import calc_SC_RDKit
import frag_utils


# Setup

if len(sys.argv) != 13:
    print("Not provided all arguments")
    quit()

data_set = sys.argv[1] # Options: ZINC, CASF
gen_smi_file = sys.argv[2] # Path to generated molecules
reference_sdf_loc = sys.argv[3] # Path to SDF containing reference molecules
train_set_path = sys.argv[4] # Path to training set
save_path = sys.argv[5] # Save path
version = sys.argv[6] # Name
n_cores = int(sys.argv[7]) # Number of cores to use
verbose = bool(sys.argv[8]) # Output results
if sys.argv[9] == "None":
    restrict = None
else:
    restrict = int(sys.argv[9]) # Set to None if don't want to restrict
pains_smarts_loc = sys.argv[10] # Path to PAINS SMARTS
data_path = sys.argv[11] # Path to data path
PROTACs = bool(sys.argv[12]) # Whether generate PROTACs or not

if verbose:
    print("##### Start Settings #####")
    print("Data set:", data_set)
    print("Generated smiles file:", gen_smi_file)
    print("Reference SDF:", reference_sdf_loc)
    print("Training set:", train_set_path)
    print("Save path:", save_path)
    print("Version:", version)
    print("Number of cores:", n_cores)
    print("Verbose:", verbose)
    print("Restrict data:", restrict)
    print("PAINS SMARTS location:", pains_smarts_loc)
    print("Data path:", data_path)
    print("PROTACs:", PROTACs)
    print("#####  End Settings  #####")


# Make output file
if not os.path.exists(save_path+version):
    os.makedirs(save_path+version)

# Prepare data

# Load molecules
# FORMAT: (Starting fragments (SMILES), Original molecule (SMILES), Generated molecule (SMILES))
generated_smiles = frag_utils.read_triples_file(gen_smi_file)

if restrict is not None and int(restrict) > 0:
        generated_smiles = generated_smiles[:restrict]

if verbose:
    print("Number of generated SMILES: %d" % len(generated_smiles))


in_mols = [smi[1] for smi in generated_smiles]
frag_mols = [smi[0] for smi in generated_smiles]
gen_mols = [smi[2] for smi in generated_smiles]

# Remove dummy atoms from starting points
clean_frags = Parallel(n_jobs=n_cores)(delayed(frag_utils.remove_dummys)(smi) for smi in frag_mols)


# Check valid
results = []
for in_mol, frag_mol, gen_mol, clean_frag in zip(in_mols, frag_mols, gen_mols, clean_frags):
    if len(Chem.MolFromSmiles(gen_mol).GetSubstructMatch(Chem.MolFromSmiles(clean_frag)))>0:
        results.append([in_mol, frag_mol, gen_mol, clean_frag])

# if verbose:
#     print("Number of generated SMILES: \t%d" % len(generated_smiles))
#     print("Number of valid SMILES: \t%d" % len(results))
#     print("%% Valid: \t\t\t%.2f%%" % (len(results)/len(generated_smiles)*100))


# Determine linkers of generated molecules
linkers = Parallel(n_jobs=n_cores)(delayed(frag_utils.get_linker)(Chem.MolFromSmiles(m[2]), Chem.MolFromSmiles(m[3]), m[1])                                    for m in results)

# Standardise linkers
for i, linker in enumerate(linkers):
    if linker == "":
        continue
    try:
        linker_canon = Chem.MolFromSmiles(re.sub('[0-9]+\*', '*', linker))
        Chem.rdmolops.RemoveStereochemistry(linker_canon)
        linkers[i] = MolStandardize.canonicalize_tautomer_smiles(Chem.MolToSmiles(linker_canon))
    except:
        continue
    
# Update results
for i in range(len(results)):
    results[i].append(linkers[i])
    
# Update results according to no oxygen linking to one oxygen or two nitrogens 
# or two sulfurs or one nitrogen one sulfur directly
results_copy = copy.deepcopy(results)
for result in results_copy:
    gen_mol = Chem.MolFromSmiles(result[2])
    for atom in gen_mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
            counter = Counter(neighbor_list)
            if 'O' in neighbor_list:
                results = [res for res in results if res != result]
                break
            elif counter['N'] == 2:
                results = [res for res in results if res != result]
                break
            elif counter['S'] == 2:
                results = [res for res in results if res != result]
                break
            elif 'N' in neighbor_list and 'S' in neighbor_list:
                results = [res for res in results if res != result]
                break
    
# Update results according to no nitrogen linking to two oxygens or one oxygen one nitrogen 
# or two sulfurs or one oxygen one sulfur or one nitrogen one sulfur or nitrogens without ring directly
results_copy = copy.deepcopy(results)
for result in results_copy:
    gen_mol = Chem.MolFromSmiles(result[2])
    for atom in gen_mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
            counter = Counter(neighbor_list)
            if counter['O'] == 2:
                results = [res for res in results if res != result]
                break
            elif 'N' in neighbor_list and 'O' in neighbor_list:
                results = [res for res in results if res != result]
                break
            elif counter['S'] == 2:
                results = [res for res in results if res != result]
                break
            elif 'O' in neighbor_list and 'S' in neighbor_list:
                results = [res for res in results if res != result]
                break
            elif 'N' in neighbor_list and 'S' in neighbor_list:
                results = [res for res in results if res != result]
                break
            elif counter['N'] == 2 and not atom.IsInRing():
                results = [res for res in results if res != result]
                break

# Update results according to no sulfur linking to one sulfur or two nitrogens or
# one oxygen but no double bonds
results_copy = copy.deepcopy(results)
for result in results_copy:
    gen_mol = Chem.MolFromSmiles(result[2])
    for atom in gen_mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            neighbor_list = [x.GetSymbol() for x in atom.GetNeighbors()]
            counter = Counter(neighbor_list)
            if 'S' in neighbor_list:
                results = [res for res in results if res != result]
                break
            elif counter['N'] == 2:
                results = [res for res in results if res != result]
                break
            elif 'O' in neighbor_list and rdkit.Chem.rdchem.BondType.DOUBLE not in [bond.GetBondType() for bond in atom.GetBonds()]:
                results = [res for res in results if res != result]
                break

# Update results according to molecular weight (MW <= 500 )
if not PROTACs:
    results_copy = copy.deepcopy(results)
    for result in results_copy:
        gen_mol = Chem.MolFromSmiles(result[2])
        molwt = Descriptors.MolWt(gen_mol)
        if molwt > 500:
            results = [res for res in results if res != result]
print('='*50)

# Update results according to the pharmacophore of linkers
with open(data_path + '/molecules_' + version + '_out.json', 'r') as json_file:
    json_dict = json.load(json_file)[0]
linker_pharm = json_dict['abs_dist'][-3:]
print('linker pharmacophore: ', linker_pharm)

results_copy = copy.deepcopy(results)
for result in results_copy:
    linker_mol = Chem.MolFromSmiles(result[-1])
    pharms = ['Donor', 'Acceptor', 'Aromatic']
    dic = {'Donor': 0, 'Acceptor': 0, 'Aromatic': 0}
    dic_list = []
    if linker_mol is not None:
        feats = factory.GetFeaturesForMol(linker_mol)
        for feat in feats:
            if feat.GetFamily() in pharms:
                dic[feat.GetFamily()] += 1
        dic_list.append(dic['Donor'])
        dic_list.append(dic['Acceptor'])
        dic_list.append(dic['Aromatic'])
        if dic_list != linker_pharm:
            results = [res for res in results if res != result]
            print(result[-1] + ' ' + str(dic_list) + ' is removed.')
        else:
            print(result[-1] + ' ' + str(dic_list) + ' is kept.')

if verbose:
    print("Number of generated SMILES: \t%d" % len(generated_smiles))
    print("Number of valid SMILES: \t%d" % len(results))
    print("%% Valid: \t\t\t%.2f%%" % (len(results)/len(generated_smiles)*100))

# Prepare training set database

# Load ZINC training set
linkers_train = []

with open(train_set_path, 'r') as f:
    for line in f:
        toks = line.strip().split(' ')
        linkers_train.append(toks[1])
        
if verbose:
    print("Number of training examples: %d" % len(linkers_train))


# Prepare unique set of linkers

# Remove stereochemistry
linkers_train_nostereo = []
for smi in list(set(linkers_train)):
    mol = Chem.MolFromSmiles(smi)
    Chem.RemoveStereochemistry(mol)
    linkers_train_nostereo.append(Chem.MolToSmiles(Chem.RemoveHs(mol)))
    
# Standardise / canonicalise training set linkers
linkers_train_nostereo = {smi.replace(':1', '').replace(':2', '') for smi in set(linkers_train_nostereo)}
linkers_train_canon = []
for smi in list(linkers_train_nostereo):
    linkers_train_canon.append(MolStandardize.canonicalize_tautomer_smiles(smi))

# Remove duplicates
linkers_train_canon_unique = list(set(linkers_train_canon))

if verbose:
    print("Number of unique linkers: %d" % len(linkers_train_canon_unique))

# 2D analysis
if len(results) > 0:
    # Create dictionary of results
    results_dict = {}
    for res in results:
        if res[0]+'.'+res[1] in results_dict: # Unique identifier - starting fragments and original molecule
            results_dict[res[0]+'.'+res[1]].append(tuple(res))
        else:
            results_dict[res[0]+'.'+res[1]] = [tuple(res)]

    # Check number of unique molecules
    if verbose:
        print("Unique molecules: %.2f%%" % (frag_utils.unique(results_dict.values())*100))

    # Check novelty of generated molecules
    count_novel = 0
    for res in results:
        if res[4] in linkers_train_canon_unique:
            continue
        else:
            count_novel +=1
            
    if verbose:
        print("Novel linkers: %.2f%%" % (count_novel/len(results)*100))

    # Check proportion recovered
    recovered = frag_utils.check_recovered_original_mol(list(results_dict.values()))
    if verbose:
        print("Recovered: %.2f%%" % (sum(recovered)/len(results_dict.values())*100))

    # Check if molecules pass 2D filters 
    filters_2d = frag_utils.calc_filters_2d_dataset(results, pains_smarts_loc=pains_smarts_loc, n_cores=n_cores)

    results_filt = []
    for res, filt in zip(results, filters_2d):
        if filt[0] and filt[1] and filt[2]:
            results_filt.append(res)
            
    if verbose:
        print("Pass all 2D filters: \t\t%.2f%%" % (len(results_filt)/len(results)*100))
        print("Valid and pass all 2D filters: \t%.2f%%" % (len(results_filt)/len(generated_smiles)*100))
        print("Pass synthetic accessibility (SA) filter: \t%.2f%%" % (len([f for f in filters_2d if f[0]])/len(filters_2d)*100))
        print("Pass ring aromaticity filter: \t\t\t%.2f%%" % (len([f for f in filters_2d if f[1]])/len(filters_2d)*100))
        print("Pass SA and ring filters: \t\t\t%.2f%%" % (len([f for f in filters_2d if f[0] and f[1]])/len(filters_2d)*100))
        print("Pass PAINS filters: \t\t\t\t%.2f%%" % (len([f for f in filters_2d if f[2]])/len(filters_2d)*100))

    # with open('unlabeled_smiles.smi', 'w') as fw:
    #     for i in results_filt:
    #         fw.write(i[2] + '\n')

    # Generate conformers
    _ = rdkit_conf_parallel.gen_confs([res[2] for res in results_filt], save_path+version+"/generated_mols.sdf", 
                                     smi_frags=[res[1] for res in results_filt], numcores=n_cores, jpsettings=True)

    # generated_sdfs = save_path+version+"/RMSD/generated_mols.sdf"
    # writer = Chem.SDWriter(save_path+version+"/generated_mols.sdf")
    # gen_filter_mols = [res[2] for res in results_filt]
    # for i in range(len(Chem.SDMolSupplier(generated_sdfs))):
    #     if Chem.MolToSmiles(Chem.SDMolSupplier(generated_sdfs)[i]) in gen_filter_mols:
    #         writer.write(Chem.SDMolSupplier(generated_sdfs)[i])

    # 3D analysis

    # Load conformers
    gen_sdfs = Chem.SDMolSupplier(save_path+version+"/generated_mols.sdf")
    ref_sdfs = Chem.SDMolSupplier(reference_sdf_loc) # Pre-computed reference conformers
    ref_mol = Chem.SDMolSupplier(reference_sdf_loc)[0]

    # Get indices and compounds
    if verbose:
        print("Loading reference conformers.")
    # ZINC
    if data_set == "ZINC":
        ref_names = []
        ref_idxs = []
        for j, ref_mol in enumerate(ref_sdfs):
            if ref_mol.GetProp("_Name") not in ref_names:
                ref_names.append(ref_mol.GetProp("_Name"))
                ref_idxs.append(j)
        # For ZINC dataset, test conformers are towards the end
        ref_names = ref_names[-1000:]
        ref_idxs = ref_idxs[-1000:]
        if verbose:
            print("Loaded ZINC reference conformers.")
    # CASF
    elif data_set == "CASF":
        ref_names = []
        ref_idxs = []
        errors = 0
        for num, sdf in enumerate(ref_sdfs):
            try:
                ref_names.append(Chem.MolToSmiles(sdf))
                ref_idxs.append(num)
            except:
                errors +=1
        if verbose:
            print("Loaded CASF refernce conformers. %d errors" % errors)

    # Get list of starting fragments and original molecules
    used = set([])
    ref_identifiers = [(res[1], res[0]) for res in results_filt if res[1]+'.'+res[0] not in used and (used.add(res[1]+'.'+res[0]) or True)]

    ref_smiles = [res[1] for res in ref_identifiers]
    ref_starting_points = [res[0] for res in ref_identifiers]

    # Get indices of compounds in SD file
    start_stop_idxs = []
    start = 0
    errors = 0

    curr_st_pt = ""
    for count, gen_mol in enumerate(gen_sdfs):
        try:
            # Check if seen this ligand before
            if gen_mol.GetProp("_Model") == str(0):
                stop = count
                if curr_st_pt in ref_starting_points:
                    start_stop_idxs.append((start, stop))
                start = int(stop) # deep copy
                curr_st_pt = gen_mol.GetProp("_StartingPoint")
        except:
            errors += 1
            continue

    # Add last
    start_stop_idxs.append((start, len(gen_sdfs)))
        
    if verbose:
        print("Number of molecules passed 2D filters: \t%d" % len(results_filt))
        print("Number of molecules with conformers: \t%d" % len(start_stop_idxs))
        print("Number of errors: \t\t\t%d" % errors)

    # Calculate SC_RDKit full scores
    names_full = []
    best_scores_full = []
    idx_best_poses_full = []
    names_start_pts_full = []

    with Parallel(n_jobs=n_cores, backend='multiprocessing') as parallel:
        for i in range(-(-len(start_stop_idxs)//n_cores)):
            jobs = []
            for core in range(n_cores):
                if i*n_cores+core < len(start_stop_idxs) and gen_sdfs is not None:
                    start, stop = start_stop_idxs[i*n_cores+core]
                    # Get ref mol
                    if gen_sdfs[start] is not None:
                        frag_smi = gen_sdfs[start].GetProp("_StartingPoint")
                        # idx_direct = ref_names.index(ref_smiles[ref_starting_points.index(frag_smi)])
                        # ref_idx_sdf = ref_idxs[idx_direct]
                        # ref_mol = Chem.Mol(ref_sdfs[ref_idx_sdf])
                        # Prepare jobs
                        gen_mols = [(Chem.Mol(gen_sdfs[idx]), Chem.Mol(ref_mol)) for idx in range(start, stop) if gen_sdfs[idx] is not None] # Test addition
                        jobs.append(gen_mols)

            # Get SC_RDKit scores
            set_scores = parallel((delayed(frag_utils.SC_RDKit_full_scores)(gen_mols) for gen_mols in jobs)) # Multiprocessing step
            for core, scores in enumerate(set_scores):
                start, stop = start_stop_idxs[i*n_cores+core]
                names_full.append(gen_sdfs[start].GetProp("_Name"))
                names_start_pts_full.append(gen_sdfs[start].GetProp("_StartingPoint"))
                best_scores_full.append(max(scores))
                idx_best_poses_full.append((np.argmax(scores)+start, 0))

    # Save output
    np.save(save_path+version+"/names_full_"+version, names_full)
    np.save(save_path+version+"/best_scores_full_"+version, best_scores_full)
    np.save(save_path+version+"/idx_best_poses_full_"+version, idx_best_poses_full)
    np.save(save_path+version+"/names_start_pts_full_"+version, names_start_pts_full)

    if verbose:
        print("Number of molecules to assess: \t%d" % len(start_stop_idxs))
        print("Number of scores: \t\t%d" % len(best_scores_full))
        print("Number of failed mols: \t\t%d" % len([s for s in best_scores_full if s<0]))

    if len(best_scores_full) != 0:
        # Print SC_RDKit Full results
        print("Average SC_RDKit Full score: %.3f +- %.3f\n" % (np.mean(best_scores_full), np.std(best_scores_full)))

        thresholds_SC_RDKit = [0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
        for thresh in thresholds_SC_RDKit:
            print("SC_RDKit Full - Molecules above %.2f: %.2f%%" % (thresh, len([score for score in best_scores_full if score >= thresh]) / len(best_scores_full)*100))

    # Calculate SC_RDKit fragments scores
    names_frags = []
    best_scores_frags = []
    idx_best_poses_frags = []
    names_frags_start_pts = []

    with Parallel(n_jobs=n_cores, backend='multiprocessing') as parallel:
        for i in range(-(-len(start_stop_idxs)//n_cores)):
            jobs = []
            for core in range(n_cores):
                if i*n_cores+core < len(start_stop_idxs) and gen_sdfs is not None:
                    start, stop = start_stop_idxs[i*n_cores+core]
                    # Get ref mol
                    if gen_sdfs[start] is not None:
                        frag_smi = gen_sdfs[start].GetProp("_StartingPoint")
                        # idx_direct = ref_names.index(ref_smiles[ref_starting_points.index(frag_smi)])
                        # ref_idx_sdf = ref_idxs[idx_direct]
                        # ref_mol = Chem.Mol(ref_sdfs[ref_idx_sdf])
                        # Prepare jobs
                        gen_mols = [(Chem.Mol(gen_sdfs[idx]), Chem.Mol(ref_mol), str(frag_smi)) for idx in range(start, stop) if gen_sdfs[idx] is not None] # Test addition
                        jobs.append(gen_mols)
                    
            # Get SC_RDKit scores
            set_scores = parallel((delayed(frag_utils.SC_RDKit_frag_scores)(gen_mols) for gen_mols in jobs)) # Multiprocessing step
            for core, scores in enumerate(set_scores):
                start, stop = start_stop_idxs[i*n_cores+core]
                names_frags.append(gen_sdfs[start].GetProp("_Name"))
                names_frags_start_pts.append(gen_sdfs[start].GetProp("_StartingPoint"))
                best_scores_frags.append(max(scores))
                idx_best_poses_frags.append((np.argmax(scores)+start, 0))
                
    # Save output
    np.save(save_path+version+"/names_frags_"+version, names_frags)
    np.save(save_path+version+"/best_scores_frags_"+version, best_scores_frags)
    np.save(save_path+version+"/idx_best_poses_frags_"+version, idx_best_poses_frags)
    np.save(save_path+version+"/names_frags_start_pts_"+version, names_frags_start_pts)
                
    if verbose:
        print("Number of molecules to assess: \t%d" % len(start_stop_idxs))
        print("Number of scores: \t\t%d" % len(best_scores_frags))
        print("Number of failed mols: \t\t%d" % len([s for s in best_scores_frags if s<0]))

    if len(best_scores_frags) != 0:
        # Print SC_RDKit Fragments results
        print("Average SC_RDKit Fragments score: %.3f +- %.3f\n" % (np.mean(best_scores_frags), np.std(best_scores_frags)))

        thresholds_SC_RDKit = [0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
        for thresh in thresholds_SC_RDKit:
            print("SC_RDKit Fragments - Molecules above %.2f: %.2f%%" % (thresh, len([score for score in best_scores_frags if score >= thresh]) / len(best_scores_frags)*100))

    # Calculate fragments RMSDs
    names_rmsd_frags = []
    best_rmsd_frags = []
    idx_best_rmsd_poses_frags = []
    names_rmsd_frags_start_pts = []

    with Parallel(n_jobs=n_cores, backend='multiprocessing') as parallel:
        for i in range(-(-len(start_stop_idxs)//n_cores)):
            jobs = []
            for core in range(n_cores):
                if i*n_cores+core < len(start_stop_idxs) and gen_sdfs is not None:
                    start, stop = start_stop_idxs[i*n_cores+core]
                    # Get ref mol
                    if gen_sdfs[start] is not None:
                        frag_smi = gen_sdfs[start].GetProp("_StartingPoint")
                        # idx_direct = ref_names.index(ref_smiles[ref_starting_points.index(frag_smi)])
                        # ref_idx_sdf = ref_idxs[idx_direct]
                        # ref_mol = Chem.Mol(ref_sdfs[ref_idx_sdf])
                        # Prepare jobs
                        gen_mols = [(Chem.Mol(gen_sdfs[idx]), Chem.Mol(ref_mol), str(frag_smi)) for idx in range(start, stop) if gen_sdfs[idx] is not None] # Test addition
                        jobs.append(gen_mols)

            # Calculate RMSDs
            set_scores = parallel((delayed(frag_utils.rmsd_frag_scores)(gen_mols) for gen_mols in jobs)) # Multiprocessing step
            for core, scores in enumerate(set_scores):
                start, stop = start_stop_idxs[i*n_cores+core]
                names_rmsd_frags.append(gen_sdfs[start].GetProp("_Name"))
                names_rmsd_frags_start_pts.append(gen_sdfs[start].GetProp("_StartingPoint"))
                best_rmsd_frags.append(min(scores))
                idx_best_rmsd_poses_frags.append((np.argmin(scores)+start, 0))

    # Save output
    np.save(save_path+version+"/names_rmsd_frags_"+version, names_rmsd_frags)
    np.save(save_path+version+"/best_rmsd_frags_"+version, best_rmsd_frags)
    np.save(save_path+version+"/idx_best_rmsd_poses_frags_"+version, idx_best_rmsd_poses_frags)
    np.save(save_path+version+"/names_rmsd_frags_start_pts_"+version, names_rmsd_frags_start_pts)
                
    if verbose:
        print("Number of molecules to assess: \t%d" % len(start_stop_idxs))
        print("Number of scores: \t\t%d" % len(best_rmsd_frags))
        print("Number of failed mols: \t\t%d" % len([s for s in best_rmsd_frags if s>90]))

    if len(best_rmsd_frags) != 0:
        # Print RMSD Fragments results
        print("Average Fragments RMSD: %.3f +- %.3f\n" % (np.mean(best_rmsd_frags), np.std(best_rmsd_frags)))

        thresholds_rmsd = [1.0, 0.75, 0.6, 0.5]
        for thresh in thresholds_rmsd:
            print("RMSD Fragments - Molecules below %.2f: %.2f%%" % (thresh, len([score for score in best_rmsd_frags if score <= thresh]) / len(best_rmsd_frags)*100))

    # Best by SC_RDKit Full
    new_best_mols_SC_RDKit_Full = []
    new_names_full = []
    best_mols_SC_RDKit_Full = sorted(list(zip(names_full, best_scores_full, idx_best_poses_full)), key=lambda x: x[1], reverse=True)
    for best_mol_SC_RDKit_Full in best_mols_SC_RDKit_Full:
        if best_mol_SC_RDKit_Full[0] not in new_names_full:
            new_names_full.append(best_mol_SC_RDKit_Full[0])
            new_best_mols_SC_RDKit_Full.append(best_mol_SC_RDKit_Full)
    if len(new_names_full) >= 100:
        new_best_mols_SC_RDKit_Full = new_best_mols_SC_RDKit_Full[:100]
    else:
        pass
    # print('best_mols_SC_RDKit_Full: ', best_mols_SC_RDKit_Full)
    idx_list_1 = []
    for i in new_best_mols_SC_RDKit_Full:
        idx_list_1.append(i[2][0])
    print('best_idx_SC_RDKit_Full: ', idx_list_1)

    # Best by SC_RDKit Fragments
    new_best_mols_SC_RDKit_Frag = []
    new_names_frags = []
    best_mols_SC_RDKit_Frag = sorted(list(zip(names_frags, best_scores_frags, idx_best_poses_frags)), key=lambda x: x[1], reverse=True)
    for best_mol_SC_RDKit_Frag in best_mols_SC_RDKit_Frag:
        if best_mol_SC_RDKit_Frag[0] not in new_names_frags:
            new_names_frags.append(best_mol_SC_RDKit_Frag[0])
            new_best_mols_SC_RDKit_Frag.append(best_mol_SC_RDKit_Frag)
    if len(new_names_frags) >= 100:
        new_best_mols_SC_RDKit_Frag = new_best_mols_SC_RDKit_Frag[:100]
    else:
        pass
    # print('best_mols_SC_RDKit_Frag: ', best_mols_SC_RDKit_Frag)
    idx_list_2 = []
    for i in new_best_mols_SC_RDKit_Frag:
        idx_list_2.append(i[2][0])
    print('best_idx_SC_RDKit_Frag: ', idx_list_2)

    # Best by RMSD Fragments
    new_best_mols_RMSD_Frag = []
    new_names_rmsd_frags = []
    best_mols_RMSD_Frag = sorted(list(zip(names_rmsd_frags, best_rmsd_frags, idx_best_rmsd_poses_frags)), key=lambda x: x[1])
    for best_mol_RMSD_Frag in best_mols_RMSD_Frag:
        if best_mol_RMSD_Frag[0] not in new_names_rmsd_frags:
            new_names_rmsd_frags.append(best_mol_RMSD_Frag[0])
            new_best_mols_RMSD_Frag.append(best_mol_RMSD_Frag)
    if len(new_names_rmsd_frags) >= 100:
        new_best_mols_RMSD_Frag = new_best_mols_RMSD_Frag[:100]
    else:
        pass
    # print('best_mols_RMSD_Frag: ', best_mols_RMSD_Frag)
    idx_list_3 = []
    for i in new_best_mols_RMSD_Frag:
        idx_list_3.append(i[2][0])
    print('best_idx_RMSD_Frag: ', idx_list_3)
    print('best_mols_RMSD_Frag: ', new_best_mols_RMSD_Frag)
