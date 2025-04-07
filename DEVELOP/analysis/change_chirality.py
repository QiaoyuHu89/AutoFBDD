import rdkit
from rdkit import Chem
import sys,os

def change_chirality(smi_file):
    frags = []
    ori_smis = []
    gen_smis = []
    gen_mols = []

    with open(smi_file, 'r') as fr:
        smis_list = fr.readlines()
        for i in smis_list:
            frags.append(i.strip().split(' ')[0])
            ori_smis.append(i.strip().split(' ')[1])
            gen_smis.append(i.strip().split(' ')[2])

    for i in gen_smis:
        gen_mol = Chem.MolFromSmiles(i)
        gen_mols.append(gen_mol)

    def spam(n):
        out=[]
        for perm in getPerms(n):
            elem = [ int(i) for i in list(perm) ]
            out.append(elem)
        return out

    def getPerms(n):
        from itertools import permutations
        for i in getCandidates(n):
            for perm in set(permutations(i)):
                yield ''.join(perm)

    def getCandidates(n):
        for i in range(0, n+1):
            res = "1" * i + "0" * (n - i)
            yield res

    def GetStereoIsomers(mol):
        from rdkit import Chem
        from copy import copy
        out = []
        smis = []

        chiralCentres = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        #return the molecule object when no chiral centres where identified
        if chiralCentres == []:
            return [mol]

        #All bit permutations with number of bits equals number of chiralCentres
        elements = spam(len(chiralCentres))

        for isoId,element in enumerate(elements):       
            for centreId,i in enumerate(element):
                atomId = chiralCentres[centreId][0]
                if i == 0:
                    mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
                elif i == 1:
                    mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            outmol = copy(mol        )
            out.append(outmol)
            smis.append(Chem.MolToSmiles(mol,isomericSmiles=True))
        return out, smis

    brick1_smi = frags[0].split('.')[0].replace('*', '')
    brick2_smi = frags[0].split('.')[1].replace('*', '')
    brick1 = Chem.MolFromSmiles(brick1_smi)
    brick2 = Chem.MolFromSmiles(brick2_smi)
    gen_chiral_smis = []
    for mol in gen_mols:
        if len(GetStereoIsomers(mol)) > 1:
            for i in GetStereoIsomers(mol)[1]:
                gen_mol = Chem.MolFromSmiles(i)
                brick1_list = gen_mol.GetSubstructMatch(brick1, useChirality=True)
                brick2_list = gen_mol.GetSubstructMatch(brick2, useChirality=True)
                if len(brick1_list) != 0 and len(brick2_list) != 0:
                    gen_chiral_smis.append(Chem.MolToSmiles(gen_mol))
    print(gen_chiral_smis)

    with open(smi_file, 'w') as fw:
        for gen_chiral_smi in gen_chiral_smis:
            line = frags[0] + ' ' + ori_smis[0] + ' ' + gen_chiral_smi + '\n'
            print(line)
            fw.write(line)

if __name__ == "__main__":
    change_chirality(sys.argv[1])