import time
import os
from os.path import (join, exists)

from scm.plams.core.functions import read_molecules

from .attachment.dye import bob, monosubstitution


##################################          input             #####################################


# Time measuring
start = time.time()

# Path to the working folder where are prepared molecules and where folder with new coordinares
# will be made with the specific name
path = os.getcwd()
input_ligands = read_molecules(join(path, 'LIGANDS'))
input_cores = read_molecules(join(path, 'CORES'))

# Bob does what Bob has to do.
for lig in input_ligands:
    lig.guess_bonds()
    bob(lig)

# As it is written so shall it be.
for core in input_cores:
    core.guess_bonds()
    bob(core)


##############################          output folder             #################################


# Makes new folder
new_dir = 'new_molecules'
if not exists(join(path, new_dir)):
    os.makedirs(join(path, new_dir))


############################ Monosubstitution and disubstitution ##################################


# Generate structures by combining ligands and cores
monosub_molecules = monosubstitution(input_ligands, input_cores)

# Export molecules to if the minimum core/ligand distance is smaller than min_dist
min_dist = 0.0
for mol in monosub_molecules:
    if mol.properties.min_distance > min_dist:
        mol.write(join(path, new_dir, mol.properties.name + '.xyz'))
    else:
        mol.write(join(path, new_dir, 'err_' + mol.properties.name + '.xyz'))

# The End
end = time.time()
print("Elapsed wall-clock time:  ", end - start)




"""
# mono_subs contains sublists with plams molecule ('molecule') and binary ('check')
# depending on binary value, False or True, name of the coordinate file is with or without 'err'.
for item in new_molecules:
    n = ['mono', 'di']
    for molecule, check in item:
        if check:
            molecule.write(join(path, new_dir, n[new_molecules.index(item)] + '_' + molecule.properties.name + '.xyz'))
        else:
            molecule.write(join(path, new_dir, 'err_' + n[new_molecules.index(item)] + molecule.properties.name + '.xyz'))
"""


