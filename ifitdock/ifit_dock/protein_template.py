#!/bin/python
######################################################################
### Usage: $chimera/chimera --nogui --script ./ifitpre.py ###
######################################################################
import sys
import os
from chimera import runCommand

runCommand('open target.pdb')
runCommand('select solvent')
runCommand('del sel')
runCommand('select ligand')
runCommand('del sel')
runCommand('select ions')
runCommand('del sel')
runCommand('select H')
runCommand('del sel')
runCommand('write format pdb 0 target_noH.pdb')
runCommand('addh')
runCommand('addcharge ')
runCommand('write format mol2 resnum atomTypes sybyl 0 target.mol2')
runCommand('write format pdb 0 target.pdb')
runCommand('select HC')
runCommand('del sel')
runCommand('write format pdb 0 target_gbsa.pdb')
runCommand('close all')
