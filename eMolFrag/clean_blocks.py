import os
import shutil
from pymol import cmd

if cmd._COb is None:
    import pymol2
    import pymol.invocation
    pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
    pymol2.SingletonPyMOL().start()

name_list = os.listdir('KRAS_G12D_out/clean-brick')
os.chdir('KRAS_G12D_out/')
if not os.path.exists('clean-brick-sdf'):
    os.mkdir('clean-brick-sdf')
os.chdir('clean-brick')
for name in name_list:
    name_sdf = name.split('.')[0] + '.' + name.split('.')[1] + '.sdf'
    # cmd.load(name)
    # cmd.save(name_pdb)
    # cmd.delete('all')

    # cmd.load(name_pdb)
    # cmd.select('dummies', 'name *')
    # cmd.remove('dummies')
    # cmd.save(name_pdb)
    # cmd.delete('all')
    os.system('babel ' + name + ' ' + name_sdf)
    shutil.move(name_sdf, '../clean-brick-sdf')
    print('Clean ' + name + ' successfully!')
