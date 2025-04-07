import os
import argparse

def main(input_name):
    with open('addC_template.py', 'r') as fr:
        lines = fr.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace('ligand', input_name)
    chimera_ligand = input_name + '.py'
    with open(chimera_ligand, 'w') as fw:
        fw.writelines(lines)

    os.system('$chimera/chimera --nogui --script ./' + chimera_ligand)
    os.system('rm ' + chimera_ligand)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pass parameters!')
    parser.add_argument('-i', type=str, help='input file')
    args = parser.parse_args()
    input_file = args.i
    input_name = input_file.split('.')[0]
   
    main(input_name)
