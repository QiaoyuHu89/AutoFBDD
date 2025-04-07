#!/bin/python
######################################################################
### Usage: $chimera/chimera --nogui --script ./ifitpre.py ###
######################################################################
import sys
import os
from chimera import runCommand

runCommand('open ligand.sdf')
runCommand('addh')
runCommand('write format pdb atomTypes gaff 0 ligand.pdb')
