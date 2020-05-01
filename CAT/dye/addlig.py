"""
CAT.dye.addlig
==============

TODO: Insert a short description of the module

Index
-----
.. currentmodule:: CAT.dye.addlig
.. autosummary::
    add_ligands
    export_dyes

API
---
.. autofunction:: get_lig_charge
.. autofunction:: export_dyes

"""

import time
import os
from itertools import chain
from typing import Iterator, Iterable
from os.path import join, exists

from scm.plams import read_molecules, Molecule

import CAT
from CAT.attachment.dye import bob_ligand, bob_core, substitution
from CAT.attachment.substitution_symmetry import del_equiv_structures

__all__ = ['add_ligands', 'export_dyes']

"""
path = join(CAT.__path__[0], 'data')
core_dir = join(path, 'CORES')
ligand_dir = join(path, 'LIGANDS')
"""

def add_ligands(core_dir: str, ligand_dir: str, min_dist: float = 1.2) -> Iterator[Molecule]:
    # Time measuring
    start = time.time()

    # Preparing molecules = in the comment section of xyz file, write serial numbers of atoms that will be substituted

    # Path to the working folder where are prepared molecules and where folder with new coordinares
    # will be made with the specific name
    input_ligands = list(read_molecules(ligand_dir).values())
    input_cores = list(read_molecules(core_dir).values())

    # Bob does what Bob has to do (numbering the ligands)
    bob_ligand(input_ligands)

    # As it is written so shall it be (reading the comment section)
    for core in input_cores:
        core.guess_bonds()
        bob_core(core)

    # Criteria for the minimal interatomic distances
    min_dist = 1.2

    ############################            Substitution              ##################################

    # Generate structures by combining ligands and cores
    # If substitution is symmetric, produces equivalent structures, those should be deleted

    mono = substitution(input_ligands, input_cores, min_dist)

    di = substitution(input_ligands, mono, min_dist)
    di_unique = del_equiv_structures(di,'linear')

    tri = substitution(input_ligands, di, min_dist)

    tetra = substitution(input_ligands, tri, min_dist)

    tetra_unique = del_equiv_structures(tetra, 'D2h')

    # Combine and flatten all new molecules into a generator
    ret = chain.from_iterable([mono, di_unique, tri, tetra_unique])

    end = time.time()
    print("Elapsed wall-clock time:  ", end - start)
    return ret


def export_dyes(mol_list: Iterable[Molecule],
                new_dir: str = 'new_molecules',
                err_dir: str = 'err_molecules',
                min_dist: float = 1.2) -> None:
    if not exists(new_dir):
        os.makedirs(new_dir)
    if not exists(err_dir):
        os.makedirs(err_dir)

    # Export molecules to folders depending on the minimum core/ligand distance
    for mol in mol_list:
        name = mol.properties.name
        mol_distance = mol.properties.min_distance

        if mol_distance > min_dist:
            mol.write(join(new_dir, f'{name}.xyz'))
        else:
            # if minimal distance is smaler than given criteria, place structures in different folder
            mol.write(join(err_dir, f'err_{name}.xyz'))
            print(f"{name}: \n Minimal distance {mol_distance} A is smaller than {min_dist} A")
