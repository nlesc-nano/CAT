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

# Preparing molecules = in the comment section of xyz file, write serial numbers of atoms that will be substituted

# Path to the working folder where are prepared molecules and where folder with new coordinares
# will be made with the specific name
path = os.getcwd()
input_ligands = read_molecules(join(path, 'Ligands'))
input_ligands = list(input_ligands.values())

input_cores = read_molecules(join(path, 'Core'))
input_cores = list(input_cores.values())

# Bob does what Bob has to do (numbering the ligands)
bob_ligand(input_ligands)

# As it is written so shall it be (reading the comment section)
for core in input_cores:
    core.guess_bonds()
    bob_core(core)

##############################          output folder             #################################


# Makes new folders
new_dir = ['The_new_PDI_structures', 'err_molecules']

if not exists(join(path, new_dir[0])):
    os.makedirs(join(path, new_dir[0]))

# Criteria for the minimal interatomic distances
min_dist = 1.2

############################            Substitution              ##################################

# Generate structures by combining ligands and cores
# If substitution is symmetric, produces equivalent structures, those should be deleted

mono = substitution(input_ligands, input_cores, min_dist)

di = substitution(input_ligands, mono, min_dist)
di_unique = del_equiv_structures(di, 'linear')

#tri = substitution(input_ligands, di, min_dist)

#tetra = substitution(input_ligands, tri, min_dist)
#tetra_unique = del_equiv_structures(tetra, 'D2h')

# Combine and flatten all new molecules into a generator
new_molecules = chain.from_iterable([mono, di_unique])#, tri, tetra_unique])


# Export molecules to folders depending on the minimum core/ligand distance
for mol in new_molecules:
    if mol.properties.min_distance > min_dist:
        mol.write(join(path, new_dir[0], mol.properties.name + '.xyz'))
    else:
        # if minimal distance is smaler than given criteria, place structures in different folder
        if not exists(join(path, new_dir[1])):
            os.makedirs(join(path, new_dir[1]))
        mol.write(join(path, new_dir[1], 'err_' + mol.properties.name + '.xyz'))
        print (mol.properties.name + ": \n Minimal distance %f A is smaler than %f A" %
               (mol.properties.min_distance, min_dist))

# The End
end = time.time()
print("Elapsed wall-clock time:  ", end - start)
