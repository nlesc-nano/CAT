""" A module designed for finding ligand functional groups. """

__all__ = ['init_ligand_anchoring']

from itertools import chain

import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import (get_time, get_template)
from ..mol_utils import separate_mod
from ..data_handling.input_sanitizer import santize_smiles


def init_ligand_anchoring(ligand_list, arg):
    """ """
    ret = []
    for lig in ligand_list:
        if not lig.properties.dummies:
            ret += find_substructure(lig, split=arg.optional.ligand.split)
        else:
            if len(lig.properties.dummies) == 1:
                lig.properties.dummies = lig.properties.dummies[0] - 1
                split = False
            elif len(lig.properties.dummies) == 2:
                lig.properties.dummies = [i - 1 for i in lig.properties.dummies]
                split = True
            ret += [substructure_split(lig, lig.properties.dummies, split=split)]
    return ret


def find_substructure(ligand, split=True):
    """ Identify interesting functional groups within the ligand.

    ligand <plams.Molecule>: The ligand molecule.
    split <bool>: If a functional group should be split from the ligand (True) or not (False).

    return <list>[<plams.Molecule>]: A copy of the ligand for each identified functional group.
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

    ligand <plams.Molecule>: The ligand molecule.
    idx <list>[<int>, <int>]: A list of atomic indices associated with a functional group.
    split <bool>: If a functional group should be split from the ligand (True) or not (False).

    return <plams.Molecule>: The ligand molecule.
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
