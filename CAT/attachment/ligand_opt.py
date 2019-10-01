"""
CAT.attachment.ligand_opt
=========================

A module designed for optimizing the geometry of ligands.

Index
-----
.. currentmodule:: CAT.attachment.ligand_opt
.. autosummary::
    init_ligand_opt
    _parse_overwrite
    read_data
    start_ligand_jobs
    optimize_ligand
    _ligand_to_db
    split_mol
    get_frag_size
    get_dihed
    set_dihed

API
---
.. autofunction:: init_ligand_opt
.. autofunction:: _parse_overwrite
.. autofunction:: read_data
.. autofunction:: start_ligand_jobs
.. autofunction:: optimize_ligand
.. autofunction:: _ligand_to_db
.. autofunction:: get_frag_size
.. autofunction:: split_bond
.. autofunction:: get_dihed
.. autofunction:: set_dihed

"""

import itertools
from typing import (List, Tuple, Dict, Any, Union)

import numpy as np
import pandas as pd

from scm.plams import (Molecule, Atom, Bond, Settings, MoleculeError, add_to_class, Units)
from scm.plams.recipes.global_minimum import global_minimum_scan_rdkit
import scm.plams.interfaces.molecule.rdkit as molkit

import rdkit
from rdkit.Chem import AllChem

from .mol_split_cm import SplitMol
from .remove_atoms_cm import RemoveAtoms
from ..logger import logger
from ..mol_utils import (to_symbol, fix_carboxyl, get_index, round_coords, from_rdmol, to_atnum)
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.mol_to_file import mol_to_file

__all__ = ['init_ligand_opt']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
OPT = ('opt', '')
FORMULA = ('formula', '')
HDF5_INDEX = ('hdf5 index', '')
SETTINGS1 = ('settings', '1')


def init_ligand_opt(ligand_df: SettingsDataFrame) -> None:
    """Initialize the ligand optimization procedure.

    Performs an inplace update of **ligand_df**.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of valid ligands.

    """
    settings = ligand_df.settings.optional
    db = settings.database.db
    overwrite = db and 'ligand' in settings.database.overwrite
    read = db and 'ligand' in settings.database.read
    write = db and 'ligand' in settings.database.write
    optimize = settings.ligand.optimize
    lig_path = settings.ligand.dirname
    mol_format = settings.database.mol_format

    # Searches for matches between the input ligand and the database; imports the structure
    if read:
        read_data(ligand_df, read)

    if write:
        _ligand_to_db(ligand_df, opt=False)

    # Identify the to be optimized ligands and optimize them
    if optimize:
        idx = _parse_overwrite(ligand_df, overwrite)
        start_ligand_jobs(ligand_df, idx)

    # Write newly optimized structures to the database
    if write and optimize:
        _ligand_to_db(ligand_df)

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    if 'ligand' in settings.database.write and optimize and mol_format:
        mol_to_file(ligand_df[MOL], lig_path, mol_format=mol_format)


def _parse_overwrite(ligand_df: SettingsDataFrame, overwrite: bool) -> Tuple[pd.Series, str]:
    """Return a series for dataframe slicing and a to-be printer message."""
    if overwrite:
        return pd.Series(True, index=ligand_df.index, name=MOL)
    return np.invert(ligand_df[OPT])


def read_data(ligand_df: SettingsDataFrame, read: bool) -> None:
    """Read ligands from the database if **read** = ``True``."""
    db = ligand_df.settings.optional.database.db
    logger.info('Pulling ligands from database')
    db.from_csv(ligand_df, database='ligand')

    for i, mol in zip(ligand_df[OPT], ligand_df[MOL]):
        if not i:
            continue
        logger.info(f'{mol.properties.name} has been pulled from the database')
    ligand_df[OPT] = ligand_df[OPT].astype(bool, copy=False)


def start_ligand_jobs(ligand_df: SettingsDataFrame, idx: pd.Series) -> None:
    """Loop over all molecules in ``ligand_df.loc[idx]`` and perform geometry optimizations."""
    if not idx.any():
        logger.info(f'No new to-be optimized ligands found\n')
        return None
    else:
        logger.info(f'Starting ligand optimization')

    for ligand in ligand_df.loc[idx, MOL]:
        logger.info(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization has started')
        try:
            optimize_ligand(ligand)
        except Exception as ex:
            logger.error(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization '
                         'has failed')
            logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)
        else:
            logger.info(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization '
                        'is successful')

    logger.info('Finishing ligand optimization\n')
    return None


def optimize_ligand(ligand: Molecule) -> None:
    """Optimize a ligand molecule."""
    # Split the branched ligand into linear fragments and optimize them individually
    bonds = split_mol(ligand)
    with SplitMol(ligand, bonds) as mol_frags:
        for mol in mol_frags:
            mol.set_dihed(180.0)

    # Find the optimal dihedrals angle between the fragments
    for bond in bonds:
        global_minimum_scan_rdkit(ligand, ligand.get_index(bond))

    # RDKit UFF can sometimes mess up the geometries of carboxylates: fix them
    fix_carboxyl(ligand)
    ligand.round_coords()
    return None


def _ligand_to_db(ligand_df: SettingsDataFrame, opt: bool = True):
    """Export ligand optimziation results to the database."""
    # Extract arguments
    settings = ligand_df.settings.optional
    db = settings.database.db
    overwrite = 'ligand' in settings.database.overwrite

    kwargs: Dict[str, Any] = {'overwrite': overwrite}
    if opt:
        kwargs['job_recipe'] = Settings({
            '1': {'key': 'RDKit_' + rdkit.__version__, 'value': 'UFF'}
        })
        kwargs['columns'] = [FORMULA, HDF5_INDEX, SETTINGS1]
        kwargs['database'] = 'ligand'
        kwargs['opt'] = True
    else:
        kwargs['columns'] = [FORMULA, HDF5_INDEX]
        kwargs['database'] = 'ligand_no_opt'

    db.update_csv(ligand_df, **kwargs)


def split_mol(plams_mol: Molecule) -> List[Bond]:
    """Split a molecule into multiple smaller fragments; returning the bonds that have to be broken.

    One fragment is created for every branch within **plams_mol**.

    Parameters
    ----------
    plams_mol : |plams.Molecule|
        The input molecule with the properties.dummies attribute.

    Returns
    -------
    :class:`list` [|plams.Bond|]
        A list of plams bonds.

    """
    def _in_ring(bond: Bond) -> bool:
        """Check if one of the atoms in **bond** is part of a ring system."""
        return (plams_mol.in_ring(bond.atom1) or plams_mol.in_ring(bond.atom2))

    def _is_valid_bond(bond: Bond) -> bool:
        """Check if one atom in **bond** has at least 3 neighbours and the other at least 2."""
        n1, n2 = plams_mol.neighbors(bond.atom1), plams_mol.neighbors(bond.atom2)
        return (len(n1) >= 3 and len(n2) >= 2) or (len(n1) >= 2 and len(n2) >= 3)

    def _get_frag_size(bond: Bond) -> int:
        """Return the size of the largest fragment were **plams_mol** to be split along **bond**."""
        return plams_mol.get_frag_size(bond, plams_mol.properties.dummies)

    # Temporary remove hydrogen atoms
    atom_gen = (at for at in plams_mol if at.atnum == 1)
    with RemoveAtoms(plams_mol, atom_gen):
        # Remove undesired bonds
        bond_list = [bond for bond in plams_mol.bonds if not _in_ring(bond)]

        # Remove even more undesired bonds
        for bond in reversed(bond_list):
            if not _is_valid_bond(bond):
                bond_list.remove(bond)

    atom_list = list(itertools.chain.from_iterable((bond.atom1, bond.atom2) for bond in bond_list))
    atom_set = {atom for atom in atom_list if atom_list.count(atom) >= 3}
    atom_dict = {atom: [bond for bond in atom.bonds if bond in bond_list] for atom in atom_set}

    # Fragment the molecule such that the functional group is on the largest fragment
    ret = []
    for at, bond_list in atom_dict.items():
        for _ in bond_list[2:]:
            frag_size = [_get_frag_size(bond) for bond in bond_list]
            idx = np.argmax(frag_size)  # The index of the largest fragment
            bond = bond_list[idx]
            ret.append(bond)
    return ret


@add_to_class(Molecule)
def get_frag_size(self, bond: Bond, atom: Atom) -> int:
    """Return the size of the fragment containing **atom** if **self** was split into two
    molecules by the breaking of **bond**.

    Parameters
    ----------
    bond : |plams.Bond|
        A PLAMS bond.

    atom : |plams.Atom|
        A PLAMS atom. The size of the fragment containg this atom will be returned.

    Returns
    -------
    :class:`int`
        The number of atoms in the fragment containing **atom**.

    """
    if bond not in self.bonds:
        raise MoleculeError('get_frag_size: The argument bond should be of type plams.Bond and '
                            'be part of the Molecule')
    elif atom not in self.atoms:
        raise MoleculeError('get_frag_size: The argument atom should be of type plams.Atom and '
                            'be part of the Molecule')

    for at in self:
        at._visited = False

    def dfs(at1, len_at=0, has_atom=False, atom=atom):
        at1._visited = True
        len_at += 1
        if at1 is atom:
            has_atom = True
        for bond in at1.bonds:
            at2 = bond.other_end(at1)
            if at2._visited:
                continue
            i, j = dfs(at2)
            len_at += i
            has_atom = has_atom or j
        return len_at, has_atom

    bond.atom1._visited = bond.atom2._visited = True
    size1, has_atom1 = dfs(bond.atom1)
    size2, _ = dfs(bond.atom2)
    for at in self.atoms:
        del at._visited

    if has_atom1:
        return size1
    return size2


def get_dihed(atoms: Tuple[Atom, Atom, Atom, Atom], unit: str = 'degree') -> float:
    """Return the dihedral angle defined by four atoms.

    Parameters
    ----------
    atoms : |Iterable| [|plams.atoms|]
        An iterable consisting of 4 PLAMS atoms

    unit : :class:`str`
        The output unit.

    Returns
    -------
    :class:`float`
        A dihedral angle expressed in **unit**.

    """
    at1, at2, at3, at4 = atoms
    vec1 = -np.array(at1.vector_to(at2))
    vec2 = np.array(at2.vector_to(at3))
    vec3 = np.array(at3.vector_to(at4))

    v1v2, v2v3 = np.cross(vec1, vec2), np.cross(vec3, vec2)
    v1v2_v2v3 = np.cross(v1v2, v2v3)
    v2_norm_v2 = vec2 / np.linalg.norm(vec2)
    epsilon = np.arctan2(v1v2_v2v3@v2_norm_v2, v1v2@v2v3)

    return Units.convert(epsilon, 'radian', unit)


@add_to_class(Molecule)
def neighbors_mod(self, atom: Atom, exclude: Union[int, str] = 1) -> List[Atom]:
    """A modified PLAMS function: Allows the exlucison of specific elements from the return list.

    Return a list of neighbors of **atom** within the molecule.
    Atoms with **atom** has to belong to the molecule.
    Returned list follows the same order as the **atom.bond** attribute.

    Parameters
    ----------
    atom : |plams.Atom|_
        The plams atom whose neighbours will be returned.

    exclude : |str|_ or |int|_
        Exclude all neighbours with a specific atomic number or symbol.

    Returns
    -------
    |list|_ [|plams.Atom|_]
        A list of all neighbours of **atom**.

    """
    exclude = to_atnum(exclude)
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom).atnum != exclude]


@add_to_class(Molecule)
def set_dihed(self, angle: float, opt: bool = True, unit: str = 'degree') -> None:
    """Change a dihedral angle into a specific value.

    Performs an inplace update of this instance.

    Parameters
    ----------
    angle : :class:`float`
        The desired dihedral angle.

    opt : :class:`bool`
        Whether or not the dihedral adjustment should be followed up by an RDKit UFF optimization.

    unit : :class:`str`
        The input unit.

    """
    angle = Units.convert(angle, unit, 'degree')
    bond_list = [bond for bond in self.bonds if bond.atom1.atnum != 1 and bond.atom2.atnum != 1
                 and bond.order == 1 and not self.in_ring(bond)]

    anchor = self.properties.dummies
    for bond in bond_list:
        n1, n2 = self.neighbors_mod(bond.atom1), self.neighbors_mod(bond.atom2)
        n1 = [atom for atom in n1 if atom != bond.atom2]
        n2 = [atom for atom in n2 if atom != bond.atom1]
        if len(n1) > 1:
            n1 = [atom for atom in n1 if (len(self.neighbors_mod(atom)) > 1 or atom is anchor)]
        if len(n2) > 1:
            n2 = [atom for atom in n2 if (len(self.neighbors_mod(atom)) > 1 or atom is anchor)]

        if n1 and n2:
            dihed = get_dihed((n1[0], bond.atom1, bond.atom2, n2[0]))
            if anchor not in bond:
                self.rotate_bond(bond, bond.atom1, angle - dihed, unit='degree')
            else:
                self.rotate_bond(bond, bond.atom1, -dihed, unit='degree')

    if opt:
        rdmol = molkit.to_rdmol(self)
        AllChem.UFFGetMoleculeForceField(rdmol).Minimize()
        self.from_rdmol(rdmol)
