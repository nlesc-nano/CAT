import time
import os
from itertools import chain
from os.path import (join, exists)

from scm.plams.core.functions import read_molecules

from CAT.attachment.dye import bob_ligand, bob_core, substitution

### ovo je cat.01 !!!

##################################          input             #####################################


# Time measuring
start = time.time()

# Path to the working folder where are prepared molecules and where folder with new coordinares
# will be made with the specific name
path = os.getcwd()
input_ligands = read_molecules(join(path, 'LIGANDS'))
input_ligands = list(input_ligands.values())
input_cores = read_molecules(join(path, 'CORES'))
input_cores = list(input_cores.values())

# Bob does what Bob has to do.
for lig in input_ligands:
    lig.guess_bonds()
    bob_ligand(lig)

# As it is written so shall it be.
for core in input_cores:
    core.guess_bonds()
    bob_core(core)

##############################          output folder             #################################


# Makes new folder
new_dir = 'new_molecules'
if not exists(join(path, new_dir)):
    os.makedirs(join(path, new_dir))

min_dist = 1.3

############################ Monosubstitution and disubstitution ##################################

# Generate structures by combining ligands and cores
mono = substitution(input_ligands, input_cores)
di = substitution(input_ligands, mono)
tri = substitution(input_ligands, di)
tetra = substitution(input_ligands, tri)

# Combine and flatten all new molecules into a generator
new_molecules = chain.from_iterable([mono, di, tri, tetra])


# Export molecules to if the minimum core/ligand distance is smaller than min_dist
for mol in new_molecules:
    if mol.properties.min_distance > min_dist:
        mol.write(join(path, new_dir, mol.properties.name + '.xyz'))
    else:
        mol.write(join(path, new_dir, 'err_' + mol.properties.name + '.xyz'))
        print (mol.properties.name + ": \tDistance between ligand and core atoms is smaler than %f" % min_dist)

# The End
end = time.time()
print("Elapsed wall-clock time:  ", end - start)
