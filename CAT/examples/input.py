""" An example input file. """

import os
import yaml
import CAT.base as CAT


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = '/Users/basvanbeek/Documents/CdSe/Week_5'

# The input cores from path/core/
input_cores = yaml.load("""
-   - Cd68Se55.xyz
    - guess_bonds: False
    - core_indices: [130, 149, 137, 143]
""")

# Cd176Se147.xyz    - core_indices: [338, 375, 327, 372, 341, 355, 374, 356, 340]
# Cd68Se55.xyz    - core_indices: [130, 149, 137, 143]


# The input ligands from path/ligand/
input_ligands = yaml.load("""
- OC(C)=O
""")

# Optional arguments: these can safely be left to their default values
argument_dict = {
    'dir_name_list': ('core', 'ligand', 'QD'),
    'dummy': 'Cl',
    'use_database': False,
    'ligand_opt': True,
    'ligand_crs': False,
    'qd_opt': True,
    'qd_int': False,
    'qd_dissociate': {
        'job1': None, 's1': None, 
        'job2': None, 's2': None
    },
    'maxiter': 1000,
    'split': True,
}

# Runs the script: the ligand, core and quantum dot lists are returned
qd_list, core_list, ligand_list = CAT.prep(input_ligands, input_cores, path, argument_dict)
