"""A module designed for finding ligand functional groups."""

from itertools import chain
from typing import (Sequence, List, Tuple)

import pandas as pd

from scm.plams import Molecule
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..properties_dataframe import PropertiesDataFrame
from ..utils import (get_time, get_template)
from ..mol_utils import separate_mod
from ..data_handling.input_sanitizer import santize_smiles

__all__ = ['init_ligand_anchoring']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
FORMULA = ('mol', '')
HDF5_INDEX = ('hdf5 index', '')
OPT = ('opt', '')


def init_ligand_anchoring(ligand_df: PropertiesDataFrame):
    """Initialize the ligand functional group searcher.

    Parameters
    ----------
    ligand_df : |CAT.PropertiesDataFrame|_
        A dataframe of valid ligands.

    Returns
    -------
    |CAT.PropertiesDataFrame|_
        A dataframe of ligands with functional groups that can serve as valid anchor points.

    """
    # Find all functional groups; return a copy of each mol for each functional group
    mol_list = []
    for lig in ligand_df[MOL]:
        # Functional group search
        if not lig.properties.dummies:
            split = ligand_df.properties.optional.ligand.split
            mol_list += find_substructure(lig, split=split)

        # Manual specification of a functional group
        else:
            if len(lig.properties.dummies) == 1:  # optional.ligand.split = False
                lig.properties.dummies = lig.properties.dummies[0] - 1
                split = False
            elif len(lig.properties.dummies) == 2:  # optional.ligand.split = True
                lig.properties.dummies = [i - 1 for i in lig.properties.dummies]
                split = True
            mol_list += [substructure_split(lig, lig.properties.dummies, split=split)]

    # Convert the results into a dataframe
    return _get_df(mol_list)


def _get_df(mol_list: Sequence[Molecule]) -> PropertiesDataFrame:
    """Create and return a new ligand dataframe.

    Parameters
    ----------
    mol_list : |list|_ [|plams.Molecule|_]
        A list of PLAMS molecules.

    Returns
    -------
    |CAT.PropertiesDataFrame|_
        A dataframe of ligands with functional groups that can serve as valid anchor points.

    """
    # Create the dataframe index and columns
    idx_tuples = [(mol.properties.smiles, mol.properties.anchor) for mol in mol_list]
    idx = pd.MultiIndex.from_tuples(idx_tuples, names=['smiles', 'anchor'])
    columns_tuples = [MOL, FORMULA, HDF5_INDEX, OPT]
    columns = pd.MultiIndex.from_tuples(columns_tuples, names=['index', 'sub index'])

    # Create, fill and return the dataframe
    df = pd.DataFrame(-1, index=idx, columns=columns)
    df[MOL] = mol_list
    df[FORMULA] = [lig.get_formula() for lig in df[MOL]]
    df[OPT] = False
    return df


def find_substructure(ligand: Molecule,
                      split: bool = True) -> List[Molecule]:
    """Identify interesting functional groups within the ligand.

    Parameters
    ----------
    ligand : |plams.Molecule|_
        The ligand molecule.

    split : bool
        If a functional group should be split from **ligand** (``True``) or not (``False``).

    Returns
    -------
    |list|_ [|plams.Molecule|_]
        A list of ligands.
        A single copy of **ligand** is created for each identified functional group,
        removing parts of the functional group if required (see **split**).
        An empty list is returned if no valid functional groups are found.

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
        msg = 'No functional groups were found (optional.ligand.split = {}) for ligand: {}'
        print(get_time() + msg.format(split, ligand.properties.smiles))
        ligand_list = []

    return ligand_list


def substructure_split(ligand: Molecule,
                       idx: Tuple[int, int],
                       split: bool = True) -> Molecule:
    """Delete the hydrogen or mono-/polyatomic counterion attached to the functional group.

    Sets the charge of the remaining heteroatom to -1 if ``split=True``.

    Parameters
    ----------
    ligand: |plams.Molecule|_
        The ligand molecule.

    idx : |tuple|_ [|int|_]
        A tuple with 2 atomic indices associated with a functional group.

    split : bool
        If a functional group should be split from **ligand** (``True``) or not (``False``).

    Returns
    -------
    |plams.Molecule|_
        A copy of **ligand**, with part of its functional group removed (see **split**).

    """
    lig = ligand.copy()
    at1 = lig[idx[0] + 1]
    at2 = lig[idx[-1] + 1]

    if split:
        lig.delete_atom(at2)
        mol_list = lig.separate_mod()
        for mol in mol_list:
            if at1 in mol:
                lig = mol
                break

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
    lig.properties.path = ligand.properties.path

    return lig
