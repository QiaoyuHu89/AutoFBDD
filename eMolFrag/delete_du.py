import os
import glob
from pymol import cmd

if cmd._COb is None:
    import pymol2
    import pymol.invocation
    pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
    pymol2.SingletonPyMOL().start()
    
bricks = glob.glob('test_out/output-brick/*.sdf')
for brick in bricks:
    cmd.load(brick)
    cmd.select('organic_atoms', 'name H+B+C+N+O+S+P+F+Cl+Br+I')
    cmd.select('du', 'not organic_atoms')
    cmd.remove('du')
    brick_pdb = brick.replace('.sdf', '.pdb')
    brick_newpdb = brick.replace('.sdf', 'new.pdb')
    cmd.save(brick_pdb)
    cmd.delete('all')
    with open(brick_pdb, 'r') as fr:
        lines = fr.readlines()
    with open(brick_newpdb, 'w') as fw:
        for line in lines:
            if line.startswith('HETATM'):
                fw.write(line)
    cmd.load(brick_newpdb)
    cmd.save(brick)
    cmd.delete('all')
    os.system('rm ' + brick_pdb)
    os.system('rm ' + brick_newpdb)
    print(brick + ' du is deleted!')
print('du deletion is done!')