import time
import os
from itertools import chain
from os.path import (join, exists)

from scm.plams.core.functions import read_molecules

from CAT.attachment.dye import bob_ligand, bob_core, substitution
from CAT.attachment.substitution_symmetry import del_equiv_structures

##################################          input             #####################################


# Time measuring
start = time.time()

# Path to the working folder where are prepared molecules and where folder with new coordinares
# will be made with the specific name
path = os.getcwd()
input_ligands = read_molecules(join(path, 'LIGANDStest'))
input_ligands = list(input_ligands.values())
input_cores = read_molecules(join(path, 'COREStest'))
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
new_dir = ['new_molecules', 'err_molecules']

if not exists(join(path, new_dir[0])):
    os.makedirs(join(path, new_dir[0]))

min_dist = 1.2

############################            Substitution              ##################################

# Generate structures by combining ligands and cores
mono = substitution(input_ligands, input_cores, min_dist)
di = substitution(input_ligands, mono, min_dist)
di_unique = del_equiv_structures(di)
tri = substitution(input_ligands, di, min_dist)
tetra = substitution(input_ligands, tri, min_dist)
tetra_unique = del_equiv_structures(tetra)

# Combine and flatten all new molecules into a generator
new_molecules = chain.from_iterable([mono, di_unique, tri, tetra_unique])


# Export molecules to if the minimum core/ligand distance is smaller than min_dist
for mol in new_molecules:
    if mol.properties.min_distance > min_dist:
        mol.write(join(path, new_dir[0], mol.properties.name + '.xyz'))
    else:
        if not exists(join(path, new_dir[1])):
            os.makedirs(join(path, new_dir[1]))
        mol.write(join(path, new_dir[1], 'err_' + mol.properties.name + '.xyz'))
        print (mol.properties.name + ": \tDistance between ligand and core atoms is smaler than %f" % min_dist)

# The End
end = time.time()
print("Elapsed wall-clock time:  ", end - start)
