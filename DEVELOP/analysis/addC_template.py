#!/bin/python
######################################################################
### Usage: $chimera/chimera --nogui --script ./ifitpre.py ###
######################################################################
import sys
import os
from chimera import runCommand

runCommand('open ligand.sdf')
runCommand('addh')
runCommand('addcharge ')
runCommand('write format mol2 atomTypes sybyl 0 ligand.mol2')
