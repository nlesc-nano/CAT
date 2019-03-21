""" A module designed for finding ligand functional groups. """

__all__ = ['init_ligand_anchoring']

from itertools import chain

import numpy as np
import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import (get_time, get_template)
from ..mol_utils import separate_mod
from ..data_handling.input_sanitizer import santize_smiles


def init_ligand_anchoring(ligand_df, arg):
    """ Initialize the ligand functional group searcher.

    :parameter ligand_df: A dataframe of valid ligands.
    :type ligand_df: |pd.DataFrame|_ (columns: |str|_, index=|int|_, values=|plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_)
    :return: A dataframe of ligands with functional groups that can serve as valid anchor points.
    :rtype: |pd.DataFrame|_ (columns: |str|_, index=|str|_, values=|plams.Molecule|_)
    """
    # Find all functional groups; return a copy of each mol for each functional group
    mol_list = []
    for lig in ligand_df['mol']:
        if not lig.properties.dummies:  # Functional group search
            mol_list += find_substructure(lig, split=arg.optional.ligand.split)
        else:  # Manual specification of a functional group
            if len(lig.properties.dummies) == 1:  # optional.ligand.split = False
                lig.properties.dummies = lig.properties.dummies[0] - 1
                split = False
            elif len(lig.properties.dummies) == 2:  # optional.ligand.split = True
                lig.properties.dummies = [i - 1 for i in lig.properties.dummies]
                split = True
            mol_list += [substructure_split(lig, lig.properties.dummies, split=split)]

    # Convert the results into a dataframe
    return _get_df(mol_list)


def _get_df(mol_list):
    """ Create and return a new ligand dataframe.

    :parameter mol_list: A list of PLAMS molecules.
    :type mol_list: |list|_ [|plams.Molecule|_]
    :return: A dataframe of ligands with functional groups that can serve as valid anchor points.
    :rtype: |pd.DataFrame|_ (columns: |str|_, index=|str|_, values=|plams.Molecule|_)
    """
    # Create the dataframe index and columns
    idx_tuples = [(mol.properties.smiles, mol.properties.anchor) for mol in mol_list]
    idx = pd.MultiIndex.from_tuples(idx_tuples, names=['smiles', 'anchor'])
    columns_tuples = [('mol', ''), ('formula', ''), ('hdf5_index', '')]
    columns = pd.MultiIndex.from_tuples(columns_tuples, names=['index', 'sub index'])

    # Create, fill and return the dataframe
    df = pd.DataFrame(-1, index=idx, columns=columns)
    df['mol'] = mol_list
    df['formula'] = [lig.get_formula() for lig in df['mol']]
    return df


def find_substructure(ligand, split=True):
    """ Identify interesting functional groups within the ligand.

    :parameter ligand: The ligand molecule.
    :type ligand: |plams.Molecule|_
    :parameter bool split: If a functional group should be split from **ligand** (*True*)
        or not (*False*).
    :return: A list of ligands. A single copy of **ligand** is created for each identified
        functional group, removing parts of the functional group if required (see **split**).
    :rtype: |list|_ [|plams.Molecule|_].
    """
    rdmol = molkit.to_rdmol(ligand)

    # Creates a list containing predefined functional groups, each saved as an rdkit molecule
    # IMPORTANT: The first atom should ALWAYS be the atom that should attach to the core
    if split:
        func_groups = get_template('smiles.yaml').split
    else:
        func_groups = get_template('smiles.yaml').no_split
    func_groups = chain.from_iterable(func_groups.values())
    func_groups = [Chem.MolFromSmarts(smarts) for smarts in func_groups]

    # Searches for functional groups (defined by functional_group_list) within the ligand
    # Duplicates are removed
    rdmatch = rdmol.GetSubstructMatches
    matches = chain(*[rdmatch(smarts) for smarts in func_groups])

    # Remove all duplicate matches, each heteroatom (match[0]) should have <= 1 entry
    ligand_indices = []
    ref = []
    for match in matches:
        if match[0] not in ref:
            ligand_indices.append(match)
            ref.append(match[0])

    if ligand_indices:
        ligand_list = [substructure_split(ligand, tup, split) for tup in ligand_indices]
    else:
        print(get_time() + 'No functional groups were found (optional.ligand.split = ' + str(split)
              + ') for ligand: ' + ligand.properties.smiles)
        ligand_list = []

    return ligand_list


def substructure_split(ligand, idx, split=True):
    """
    Delete the hydrogen or mono-/polyatomic counterion attached to the functional group.
    Sets the charge of the remaining heteroatom to -1 if split=True.

    :parameter ligand: The ligand molecule.
    :type ligand: |plams.Molecule|_
    :parameter idx: A list of 2 atomic indices associated with a functional group.
    :type idx: 2 |list|_ [|int|_]
    :parameter bool split: If a functional group should be split from **ligand** (*True*)
        or not (*False*).
    :return: A copy of **ligand**, with part of its functional group removed (see **split**).
    :rtype: |plams.Molecule|_
    """
    lig = ligand.copy()
    at1 = lig[idx[0] + 1]
    at2 = lig[idx[-1] + 1]

    if split:
        if len(lig.separate()) == 1:
            lig.delete_atom(at2)
        else:
            mol1, mol2 = lig.separate_mod()
            if at1 in mol1:
                lig = mol1
            else:
                lig = mol2

        # Check if the ligand heteroatom has a charge assigned, assigns a charge if not
        if not at1.properties.charge or at1.properties.charge == 0:
            at1.properties.charge = -1

    # Update ligand properties
    lig.properties.dummies = at1
    lig.properties.anchor = at1.symbol + str(lig.atoms.index(at1) + 1)
    lig.properties.charge = sum(atom.properties.charge for atom in lig if atom.properties.charge)

    # Update the ligand smiles string
    rdmol = molkit.to_rdmol(lig)
    smiles = Chem.MolToSmiles(rdmol)
    lig.properties.smiles = Chem.CanonSmiles(smiles)
    lig.properties.name = santize_smiles(lig.properties.smiles) + '@' + lig.properties.anchor

    return lig
