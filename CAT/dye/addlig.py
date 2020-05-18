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
.. autofunction:: add_ligands
.. autofunction:: export_dyes

"""

import os
from itertools import chain
from typing import Iterator, Iterable, Collection
from os.path import join, exists
import numpy as np
import pickle
import gzip
import math

from scm.plams import read_molecules, Molecule
from scm.plams.interfaces.molecule.rdkit import to_rdmol

from CAT.logger import logger

from rdkit.Chem import AllChem, rdMolDescriptors, FindMolChiralCenters

import CAT
from CAT.attachment.dye import label_lig, label_core, substitution
from CAT.attachment.substitution_symmetry import del_equiv_structures

__all__ = ['add_ligands', 'export_dyes', 'SA_scores']


def add_ligands(core_dir: str,
                ligand_dir: str,
                min_dist: float = 1.2,
                n: int = 1,
                symmetry: Collection[str] = ()) -> Iterator[Molecule]:
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
    symmetry: tuple[str]
        Keywords for substitution symmetry for deleting equivalent molecules

    Returns
    -------
    Iterator[Molecule]
        New structures that are containg core and lingad fragments

    """
    # Validate the symmetry
    valid_sym = {'linear', 'triangle', 'D2h'}
    if not valid_sym.issuperset(symmetry):
        raise ValueError(f"Invalid 'symmetry' value: {symmetry!r}")

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
        for i in range(1, n):
            poli = substitution(input_ligands, mono, min_dist)
            new_mols.append(poli)
            mono = poli

    if not symmetry:
        pass
    elif 'linear' in symmetry:
        new_mols[1] = del_equiv_structures(new_mols[1], 'linear')
    elif 'triangle' in symmetry:
        new_mols[2] = del_equiv_structures(new_mols[2], 'triangle')
    elif 'D2h' in symmetry:
        new_mols[3] = del_equiv_structures(new_mols[3], 'D2h')

    # Combine and flatten all new molecules into a generator
    return chain.from_iterable(new_mols)


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
    for mol in mol_list:
        name = mol.properties.name
        mol_distance = mol.properties.min_distance

        if mol_distance > min_dist:
            mol.write(join(new_dir, f'{name}.xyz'))
        else:
            mol.write(join(err_dir, f'err_{name}.xyz'))
            logger.warning(
                f"{name}: \n Minimal distance {mol_distance} A is smaller than {min_dist} A"
            )

# The functions '_compute_SAS' and 'SA_scores' as well as data set 'SA_score.pkl.gz'
# are copied and adapted from:
###########################################################################
#    Title: MolGAN: An implicit generative model for small molecular graphs
#    Authors: De Cao, Nicola and Kipf, Thomas
#    Date: 2018
#    Availability: https://github.com/nicola-decao/MolGAN
###########################################################################

SA_model = {i[j]: float(i[0])
         for i in pickle.load(gzip.open('SA_score.pkl.gz')) for j in range(1, len(i))}


def _compute_SAS(mol):
    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0
    # for bitId, v in fps.items():
    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += SA_model.get(sfp, -4) * v
    score1 /= nf

    # features score
    nAtoms = mol.GetNumAtoms()
    nChiralCenters = len(FindMolChiralCenters(
        mol, includeUnassigned=True))
    ri = mol.GetRingInfo()
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgeheads = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms ** 1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = 0.

    # ---------------------------------------
    # This differs from the paper, which defines:
    # macrocyclePenalty = math.log10(nMacrocycles+1)
    # This form generates better results when 2 or more macrocycles are present
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0. - sizePenalty - stereoPenalty - \
           spiroPenalty - bridgePenalty - macrocyclePenalty

    # correction for the fingerprint density
    # not in the original publication, added in version 1.1
    # to make highly symmetrical molecules easier to synthetise
    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.
    # smooth the 10-end
    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    return sascore


def SA_scores(mols, norm=False):
    mols = [to_rdmol(mol) for mol in mols]
    scores = [_compute_SAS(mol) if mol is not None else None for mol in mols]
    scores = np.array(list(map(lambda x: 10 if x is None else x, scores)))
    scores = np.clip(remap(scores, 5, 1.5), 0.0, 1.0) if norm else scores

    return scores
