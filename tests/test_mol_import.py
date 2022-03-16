"""Tests for :mod:`CAT.data_handling.mol_import`."""

import random
from os.path import join

import numpy as np

from scm.plams import (Settings, Molecule)
import scm.plams.interfaces.molecule.rdkit as molkit
from assertionlib import assertion

from CAT.utils import get_formula
from CAT.data_handling.mol_import import (
    read_mol_xyz, read_mol_pdb, read_mol_mol, read_mol_smiles, read_mol_plams, read_mol_rdkit,
    read_mol_folder, read_mol_txt, get_charge_dict, set_mol_prop, canonicalize_mol
)

PATH = join('tests', 'test_files')
REF_MOL = Molecule(join(PATH, 'Methanol.xyz'))
REF_MOL.guess_bonds()
canonicalize_mol(REF_MOL)


def test_read_mol_xyz() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_xyz`."""
    xyz = join(PATH, 'Methanol.xyz')
    mol_dict = Settings({'mol': xyz, 'guess_bonds': True})
    mol = read_mol_xyz(mol_dict)

    assertion.isinstance(mol, Molecule)
    np.testing.assert_allclose(mol.as_array(), REF_MOL.as_array())
    assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_pdb() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_pdb`."""
    pdb = join(PATH, 'Methanol.pdb')
    mol_dict = Settings({'mol': pdb, 'guess_bonds': False})
    mol = read_mol_pdb(mol_dict)

    assertion.isinstance(mol, Molecule)
    assertion.lt(mol.as_array().sum() - REF_MOL.as_array().sum(), 0.01)
    assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_mol() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_mol`."""
    mol_file = join(PATH, 'Methanol.mol')
    mol_dict = Settings({'mol': mol_file, 'guess_bonds': False})
    mol = read_mol_mol(mol_dict)

    assertion.isinstance(mol, Molecule)
    assertion.lt(mol.as_array().sum() - REF_MOL.as_array().sum(), 0.01)
    assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_smiles() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_smiles`."""
    smiles = 'CO'
    mol_dict = Settings({'mol': smiles, 'guess_bonds': False})
    mol = read_mol_smiles(mol_dict)

    assertion.isinstance(mol, Molecule)
    assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_plams() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_smiles`."""
    mol = REF_MOL.copy()
    random.shuffle(mol.atoms)
    mol_dict = Settings({'mol': mol, 'guess_bonds': False})
    mol = read_mol_plams(mol_dict)

    assertion.isinstance(mol, Molecule)
    assertion.lt(mol.as_array().sum() - REF_MOL.as_array().sum(), 0.01)
    assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_rdkit() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_rdkit`."""
    mol = REF_MOL.copy()
    random.shuffle(mol.atoms)
    rdmol = molkit.to_rdmol(mol)
    mol_dict = Settings({'mol': rdmol, 'guess_bonds': False})
    mol = read_mol_rdkit(mol_dict)

    assertion.isinstance(mol, Molecule)
    assertion.lt(mol.as_array().sum() - REF_MOL.as_array().sum(), 0.01)
    assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_folder() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_folder`."""
    mol_dict = Settings({'mol': PATH, 'path': PATH, 'guess_bonds': True, 'is_core': False})
    _mol_list = read_mol_folder(mol_dict)
    mol_list = [mol for mol in _mol_list if get_formula(mol) == 'C1H4O1']

    for mol in mol_list:
        assertion.isinstance(mol, Molecule)
        assertion.lt(mol.as_array().sum() - REF_MOL.as_array().sum(), 0.01)
        assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])


def test_read_mol_txt() -> None:
    """Test :func:`CAT.data_handling.validate_input.read_mol_txt`."""
    txt = join(PATH, 'Methanol.txt')
    mol_dict = Settings({'mol': txt, 'path': PATH, 'guess_bonds': True, 'is_core': False})
    mol_list = read_mol_txt(mol_dict)

    for mol in mol_list[:-1]:
        assertion.isinstance(mol, Molecule)
        assertion.lt(mol.as_array().sum() - REF_MOL.as_array().sum(), 0.01)
        assertion.eq([at.symbol for at in mol], [at.symbol for at in REF_MOL])

    assertion.isinstance(mol_list[-1], Molecule)
    assertion.eq([at.symbol for at in mol_list[-1]], [at.symbol for at in REF_MOL])


def test_get_charge_dict() -> None:
    """Test :func:`CAT.data_handling.validate_input.get_charge_dict`."""
    charge_dict = get_charge_dict()
    ref = {
        'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
        'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
        'N': -3, 'P': -3, 'As': -3, 'Sb': -3, 'Bi': -3,
        'O': -2, 'S': -2, 'Se': -2, 'Te': -2, 'Po': -2,
        'H': -1, 'F': -1, 'Cl': -1, 'Br': -1, 'I': -1, 'At': -1,
        'Cd': 2, 'Pb': 2,
        'In': 3,
    }

    assertion.eq(charge_dict, ref)


def test_set_mol_prop() -> None:
    """Test :func:`CAT.data_handling.validate_input.set_mol_prop`."""
    mol = REF_MOL.copy()
    mol.properties = Settings()
    mol_dict = Settings({'is_core': False, 'path': PATH, 'name': 'CO'})

    set_mol_prop(mol, mol_dict)
    ref = {'name': 'CO', 'dummies': {}, 'path': PATH, 'job_path': [], 'smiles': 'CO'}
    assertion.eq(mol.properties, ref)

    ref1 = Settings({
        'stereo': {}, 'charge': 0,
        'pdb_info': {'ResidueName': 'LIG', 'Occupancy': 1.0, 'TempFactor': 0.0,
                     'ResidueNumber': 1, 'ChainId': 'A', 'IsHeteroAtom': False},
    })
    ref2 = Settings({
        'stereo': {}, 'charge': 0,
        'pdb_info': {'ResidueName': 'LIG', 'Occupancy': 1.0, 'TempFactor': 0.0,
                     'ResidueNumber': 1, 'ChainId': 'A', 'IsHeteroAtom': True},
    })
    for at in mol:
        del at.properties.pdb_info.Name
        if at.symbol == 'O':
            assertion.eq(at.properties, ref2)
        else:
            assertion.eq(at.properties, ref1)
