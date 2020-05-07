"""
CAT.dye.addlig
==============

A module for multiple compound attachment and export of the xyz files.

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

import os
from itertools import chain
from typing import Iterator, Iterable
from os.path import join, exists

from scm.plams import read_molecules, Molecule

import CAT
from CAT.attachment.dye import label_lig, label_core, substitution
from CAT.attachment.substitution_symmetry import del_equiv_structures

__all__ = ['add_ligands', 'export_dyes']


def add_ligands(core_dir: str, ligand_dir: str, min_dist: float = 1.2, n: int = 1, symmetry: Tuple[str] = []) -> Iterator[Molecule]:
    """Add ligand(s) to one core.

    Parameters
    ----------
    core_dir: str
        Name of directory where core coordinates are located
    ligand_dir: str
        Name of directory where ligands coordinates are located
    min_dist: float
        Criterion for the minimal interatomic distances
    n: int
        Number of substitutions
    symmetry: tulpe[str]
        Keywords for substitution symmetry for deleting equivalent molecules

    Returns
    -------
    Iterator[Molecule]
        New structures that are containg core and lingad fragments

    """
    input_ligands = list(read_molecules(ligand_dir).values())
    input_cores = list(read_molecules(core_dir).values())

    # Reading the comment section and labeling substitution atom and numbering the ligands
    label_lig(input_ligands)

    # Reading the comment section and labeling substitution atoms
    for core in input_cores:
        core.guess_bonds()
        label_core(core)

    # Generate structures by combining ligands and cores
    mono = substitution(input_ligands, input_cores, min_dist)
    new_mols = []
    new_mols.append(mono)
    if n > 1:
        for i in range(1,n):
            poli = substitution(input_ligands, mono, min_dist)
            new_mols.append(poli)
            mono = poli

    if symmetry == []:
        pass
    elif 'linear' in symmetry:
        new_mols[1] = del_equiv_structures(new_mols[1], 'linear')
    elif 'triangle' in symmetry:
        new_mols[2] = del_equiv_structures(new_mols[2], 'triangle')
    elif 'D2h' in symmetry:
        new_mols[3] = del_equiv_structures(new_mols[2], 'D2h')

    # Combine and flatten all new molecules into a generator
    ret = chain.from_iterable(new_mols)

    return ret


def export_dyes(mol_list: Iterable[Molecule],
                new_dir: str = 'new_molecules',
                err_dir: str = 'err_molecules',
                min_dist: float = 1.2) -> None:
    """Exports molecular coordinates to .xyz files

    Parameters
    ----------
    mol_list: Iterable[Molecule]
        List of molecules
    new_dir: str
        Name of the durectory to place new structures
    err_dir: str
        Name of the directory for molecules that do not fullfill min_dist criterion
    min_dist: float
        Criterion for the minimal interatomic distances

    """
    if not exists(new_dir):
        os.makedirs(new_dir)
    if not exists(err_dir):
        os.makedirs(err_dir)

    # Export molecules to folders depending on the minimum core/ligand distance
    print(mol_list)
    for mol in mol_list:
        name = mol.properties.name
        mol_distance = mol.properties.min_distance

        if mol_distance > min_dist:
            mol.write(join(new_dir, f'{name}.xyz'))
        else:
            mol.write(join(err_dir, f'err_{name}.xyz'))
            print(f"{name}: \n Minimal distance {mol_distance} A is smaller than {min_dist} A")
