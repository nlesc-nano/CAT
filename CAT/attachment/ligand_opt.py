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
    _ligand_to_db
    remove_duplicates
    split_bond
    neighbors_mod
    split_mol
    get_frag_size
    recombine_mol
    get_dihed
    set_dihed

API
---
.. autofunction:: init_ligand_opt
.. autofunction:: _parse_overwrite
.. autofunction:: read_data
.. autofunction:: start_ligand_jobs
.. autofunction:: _ligand_to_db
.. autofunction:: remove_duplicates
.. autofunction:: split_bond
.. autofunction:: neighbors_mod
.. autofunction:: get_frag_size
.. autofunction:: recombine_mol
.. autofunction:: get_dihed
.. autofunction:: set_dihed

"""

import itertools
from typing import (Union, Sequence, List, Tuple, Dict, Any)

import numpy as np
import pandas as pd

from scm.plams import (Molecule, Atom, Bond, Settings)
from scm.plams.core.errors import MoleculeError
from scm.plams.core.functions import add_to_class
from scm.plams.tools.units import Units
from scm.plams.recipes.global_minimum import global_minimum_scan_rdkit
import scm.plams.interfaces.molecule.rdkit as molkit

import rdkit
from rdkit.Chem import AllChem

from .ligand_attach import (rot_mol_angle, sanitize_dim_2)
from ..logger import logger
from ..settings_dataframe import SettingsDataFrame
from ..mol_utils import (to_symbol, fix_carboxyl, get_index, round_coords,
                         from_mol_other, from_rdmol, separate_mod)
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
    ligand_df[OPT] = ligand_df[OPT].astype(bool, copy=False)

    if write:
        _ligand_to_db(ligand_df, opt=False)

    # Optimize all new ligands
    if optimize:
        # Identify the to be optimized ligands
        idx = _parse_overwrite(ligand_df, overwrite)

        # Optimize the ligands
        lig_new = start_ligand_jobs(ligand_df, idx)

        # Update the ligand dataframe
        if lig_new:
            if len(lig_new) == 1:  # pd.DataFrame.loc has serious issues when assigning 1 molecue
                idx, _ = next(ligand_df[idx].iterrows())
                ligand_df.at[idx, MOL] = lig_new[0]
            else:
                ligand_df.loc[idx, MOL] = lig_new

    remove_duplicates(ligand_df)

    # Write newly optimized structures to the database
    if write and optimize:
        _ligand_to_db(ligand_df)

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    if 'ligand' in settings.database.write and optimize and mol_format:
        mol_to_file(ligand_df[MOL], lig_path, mol_format=mol_format)


def _parse_overwrite(ligand_df: SettingsDataFrame,
                     overwrite: bool) -> Tuple[pd.Series, str]:
    """Return a series for dataframe slicing and a to-be printer message."""
    if overwrite:
        return pd.Series(True, index=ligand_df.index, name=MOL)
    else:
        return np.invert(ligand_df[OPT])


def read_data(ligand_df: SettingsDataFrame,
              read: bool) -> None:
    """Read ligands from the database if **read** = ``True``."""
    db = ligand_df.settings.optional.database.db
    logger.info('Pulling ligands from database')
    db.from_csv(ligand_df, database='ligand')

    for i, mol in zip(ligand_df[OPT], ligand_df[MOL]):
        if i == -1:
            continue
        logger.info(f'{mol.properties.name} has been pulled from the database')
    ligand_df[OPT] = ligand_df[OPT].astype(bool, copy=False)


def start_ligand_jobs(ligand_df: SettingsDataFrame,
                      idx: pd.Series) -> List[Molecule]:
    """Loop over all molecules in ``ligand_df.loc[idx]`` and perform geometry optimizations."""
    if not idx.any():
        logger.info(f'No new to-be optimized ligands found\n')
        return []
    else:
        logger.info(f'Starting ligand optimization')

    lig_new = []
    for ligand in ligand_df[MOL][idx]:
        logger.info(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization has started')
        try:
            mol_list = split_mol(ligand)
            for mol in mol_list:
                mol.set_dihed(180.0)
            ligand_tmp = recombine_mol(mol_list)
            fix_carboxyl(ligand_tmp)
            ligand_tmp.round_coords()
            lig_new.append(ligand_tmp)
            logger.info(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization '
                        'is successful')
        except Exception as ex:
            logger.error(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization '
                         'has failed')
            logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)

    logger.info('Finishing ligand optimization\n')
    return lig_new


def _ligand_to_db(ligand_df: SettingsDataFrame,
                  opt: bool = True):
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


def remove_duplicates(df: pd.DataFrame) -> None:
    """Remove duplicate rows from a dataframe.

    Duplicates are identified based on their index.
    Performs an inplace update of **ligand_df**.
    """
    # Remove duplicate ligands and sort
    if df.index.duplicated().any():
        idx_name = df.index.names
        df.reset_index(inplace=True)
        i, j = idx_name
        df.drop_duplicates(subset=((i, ''), (j, '')), inplace=True)
        df.set_index(idx_name, inplace=True)
        df.index.names = idx_name
    df.sort_index(inplace=True)


@add_to_class(Molecule)
def split_bond(self, bond: Sequence[Atom],
               atom_type: Union[str, int] = 'H') -> None:
    """Delete a bond and cap the resulting fragments.

    A link to the two atoms previously defining the bond & the two capping atoms is stored under
    self.properties.mark in a list of 4-tuples.
    Performs in inplace update of **self**.

    Parameters
    ----------
    bond : |plams.Bond|_
        A PLAMS bond.

    atom_type : |str|_ or |int|_
        The atomic symbol or number of the two to be created capping atoms.

    """
    atom_type = to_symbol(atom_type)
    at1, at2 = bond.atom1, bond.atom2
    at3, at4 = Atom(symbol=atom_type, coords=at1.coords), Atom(symbol=atom_type, coords=at2.coords)
    bond_length = at3.radius + at2.radius

    self.add_atom(at3, adjacent=[at2])
    self.add_atom(at4, adjacent=[at1])
    self.bonds[-1].resize(at1, bond_length)
    self.bonds[-2].resize(at2, bond_length)
    if self.properties.mark:
        self.properties.mark.append((at1, at4, at2, at3))
    else:
        self.properties.mark = [(at1, at4, at2, at3)]
    self.delete_bond(bond)


@add_to_class(Molecule)
def neighbors_mod(self, atom: Atom,
                  exclude: Union[int, str] = 1) -> List[Atom]:
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
    exclude = to_symbol(exclude)
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom).atnum != exclude]


def split_mol(plams_mol: Molecule) -> List[Molecule]:
    """Split a molecule into multiple smaller fragments.

    One fragment is created for every branch within **plams_mol**.

    Parameters
    ----------
    plams_mol : |plams.Molecule|_
        The input molecule with the properties.dummies attribute.

    Returns
    -------
    |list|_ [|plams.Molecule|_]
        A list of one or more plams molecules.

    """
    # Temporary remove hydrogen atoms
    h_atoms = []
    h_bonds = []
    for atom in reversed(plams_mol.atoms):
        if atom.atnum == 1:
            h_atoms.append(atom)
            h_bonds.append(atom.bonds[0])
            plams_mol.delete_atom(atom)

    # Remove undesired bonds
    bond_list = [bond for bond in plams_mol.bonds if not plams_mol.in_ring(bond.atom1) and not
                 plams_mol.in_ring(bond.atom2)]

    # Remove even more undesired bonds
    for bond in reversed(bond_list):
        n1, n2 = plams_mol.neighbors_mod(bond.atom1), plams_mol.neighbors_mod(bond.atom2)
        if not (len(n1) >= 3 and len(n2) >= 2) and not (len(n1) >= 2 and len(n2) >= 3):
            bond_list.remove(bond)

    # Add the hydrogen atoms and bonds back to the molecule
    for atom, bond in zip(reversed(h_atoms), reversed(h_bonds)):
        plams_mol.add_atom(atom)
        plams_mol.add_bond(bond)

    atom_list = itertools.chain.from_iterable((bond.atom1, bond.atom2) for bond in bond_list)
    atom_set = {atom for atom in atom_list if atom_list.count(atom) >= 3}
    atom_dict = {atom: [bond for bond in atom.bonds if bond in bond_list] for atom in atom_set}

    # Fragment the molecule such that the functional group is on the largest fragment
    for at in atom_dict:
        for i in atom_dict[at][2:]:
            len_atom = [plams_mol.get_frag_size(bond, plams_mol.properties.dummies) for
                        bond in atom_dict[at]]
            idx = len_atom.index(max(len_atom))
            bond = atom_dict[at][idx]
            plams_mol.split_bond(bond)
            atom_dict[at].remove(bond)

    # Copy the properties attribute to all fragment molecules
    properties = plams_mol.properties
    mol_list = plams_mol.separate_mod()
    for mol in mol_list:
        mol.properties = properties

    return mol_list


@add_to_class(Molecule)
def get_frag_size(self, bond: Bond,
                  atom: Atom) -> int:
    """Return the size of the fragment containing **atom** if **self** was split into two
    molecules by the breaking of **bond**.

    Parameters
    ----------
    bond : |plams.Bond|_
        A PLAMS bond.

    atom : |plams.Atom|_
        A PLAMS atom. The size of the fragment containg this atom will be returned.

    Returns
    -------
    |int|_
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


def recombine_mol(mol_list: Sequence[Molecule]) -> Molecule:
    """Recombine a list of molecules into a single molecule.

    A list of 4-tuples of plams.Atoms will be read from mol_list[0].properties.mark.
    A bond will be created between tuple[0] & tuple[2]; tuple[1] and tuple[3] will be deleted.

    Parameters
    ----------
    mol_list : |list|_ [|plams.Molecule|_]
        A list of on or more plams molecules with the properties.mark atribute.

    Returns
    -------
    |plams.Molecule|_
        The (re-)merged PLAMS molecule.

    """
    if len(mol_list) == 1:
        return mol_list[0]
    tup_list = mol_list[0].properties.mark
    if not tup_list:
        raise IndexError('No PLAMS atoms specified in mol_list[0].properties.mark, '
                         'aborting recombine_mol()')

    for tup in tup_list:
        # Allign mol1 & mol2
        mol1, mol2 = tup[0].mol, tup[2].mol
        vec1 = sanitize_dim_2(tup[3]) - sanitize_dim_2(tup[2])
        vec2 = sanitize_dim_2(tup[0]) - sanitize_dim_2(tup[1])
        idx = tup[2].get_atom_index() - 1
        mol_array = rot_mol_angle(mol2, vec1, vec2, atoms_other=tup[0], idx=idx, bond_length=1.5)
        mol2.from_array(mol_array)

        # Merge mol1 & mol2
        mol1.merge_mol(mol2)
        mol1.delete_atom(tup[1])
        mol1.delete_atom(tup[3])
        mol1.add_bond(tup[0], tup[2])
        bond_tup = mol1.get_bond_index(mol1.bonds[-1])
        mol1.from_mol_other(global_minimum_scan_rdkit(mol1, bond_tup))

    del mol1.properties.mark
    return mol1


def get_dihed(atoms: Tuple[Atom, Atom, Atom, Atom],
              unit: str = 'degree') -> float:
    """Return the dihedral angle defined by four atoms.

    Parameters
    ----------
    atoms : |tuple|_ [|plams.atoms|_]
        An iterable consisting of 4 PLAMS atoms

    unit : str
        The output unit.

    Returns
    -------
    |float|_
        A dihedral angle expressed in **unit**.

    """
    vec1 = -np.array(atoms[0].vector_to(atoms[1]))
    vec2 = np.array(atoms[1].vector_to(atoms[2]))
    vec3 = np.array(atoms[2].vector_to(atoms[3]))

    v1v2, v2v3 = np.cross(vec1, vec2), np.cross(vec3, vec2)
    v1v2_v2v3 = np.cross(v1v2, v2v3)
    v2_norm_v2 = vec2 / np.linalg.norm(vec2)
    epsilon = np.arctan2(v1v2_v2v3@v2_norm_v2, v1v2@v2v3)

    return Units.convert(epsilon, 'radian', unit)


@add_to_class(Molecule)
def set_dihed(self, angle: float,
              opt: bool = True,
              unit: str = 'degree') -> None:
    """Change a dihedral angle into a specific value.

    Performs an inplace update of this instance.

    Parameters
    ----------
    angle : float
        The desired dihedral angle.

    opt : bool
        Whether or not the dihedral adjustment should be followed up by an RDKit UFF optimization.

    unit : str
        The input unit.

    """
    angle = Units.convert(angle, unit, 'degree')
    bond_list = [bond for bond in self.bonds if bond.atom1.atnum != 1 and bond.atom2.atnum != 1
                 and bond.order == 1 and not self.in_ring(bond)]

    for bond in bond_list:
        n1, n2 = self.neighbors_mod(bond.atom1), self.neighbors_mod(bond.atom2)
        n1 = [atom for atom in n1 if atom != bond.atom2]
        n2 = [atom for atom in n2 if atom != bond.atom1]
        if len(n1) > 1:
            n1 = [atom for atom in n1 if len(self.neighbors_mod(atom)) > 1]
        if len(n2) > 1:
            n2 = [atom for atom in n2 if len(self.neighbors_mod(atom)) > 1]
        if n1 and n2:
            dihed = get_dihed((n1[0], bond.atom1, bond.atom2, n2[0]))
            if self.properties.dummies not in bond:
                self.rotate_bond(bond, bond.atom1, angle - dihed, unit='degree')
            else:
                self.rotate_bond(bond, bond.atom1, -dihed, unit='degree')

    if opt:
        rdmol = molkit.to_rdmol(self)
        AllChem.UFFGetMoleculeForceField(rdmol).Minimize()
        self.from_rdmol(rdmol)
