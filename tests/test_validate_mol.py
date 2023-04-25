"""Tests for :mod:`CAT.data_handling.validate_mol`."""

from os.path import join

from scm.plams import (Settings, Molecule)
import scm.plams.interfaces.molecule.rdkit as molkit
from assertionlib import assertion

from CAT.data_handling.validate_mol import (
    validate_mol, santize_smiles, _parse_name_type, _check_core
)

PATH = join('tests', 'test_files')
MOL_PATH = join(PATH, 'Methanol.xyz')
MOL = Molecule(MOL_PATH)
MOL.guess_bonds()


def test_santize_smiles() -> None:
    """Test :func:`CAT.data_handling.validate_mol.santize_smiles`."""
    assertion.eq(santize_smiles('CO'), 'CO')
    assertion.eq(santize_smiles('C[H]O'), 'C[H]O')
    assertion.eq(santize_smiles('C(C)O'), 'C[C]O')
    assertion.eq(santize_smiles('CC=CC'), 'CC=CC')
    assertion.eq(santize_smiles('C/C=C/C'), 'trans-CC=CC')
    assertion.eq(santize_smiles(r'C/C=C\C'), 'cis-CC=CC')
    assertion.eq(santize_smiles('C/C=C/C/C=C/C'), 'trans-trans-CC=CCC=CC')
    assertion.eq(santize_smiles(r'C/C=C\C/C=C/C'), 'cis-trans-CC=CCC=CC')


def test_parse_mol_type() -> None:
    """Test :func:`CAT.data_handling.validate_mol._parse_mol_type`."""
    assertion.eq(_check_core('input_cores'), True)
    assertion.eq(_check_core('input_ligands'), False)
    assertion.assert_(_check_core, 'bob', exception=ValueError)
    assertion.assert_(_check_core, 1, exception=AttributeError)


def test_parse_name_type() -> None:
    """Test :func:`CAT.data_handling.validate_mol._parse_name_type`."""
    mol_dict = Settings({'mol': 'CCCO'})
    _parse_name_type(mol_dict)
    assertion.eq(mol_dict, {'mol': 'CCCO', 'type': 'smiles', 'name': 'CCCO'})

    mol_dict = Settings({'mol': MOL})
    _parse_name_type(mol_dict)
    assertion.eq(mol_dict, {'mol': MOL, 'type': 'plams_mol', 'name': 'Methanol'})

    mol = MOL.copy()
    mol.properties = Settings()
    mol_dict = Settings({'mol': mol})
    _parse_name_type(mol_dict)
    assertion.eq(mol_dict, {'mol': mol, 'type': 'plams_mol', 'name': 'CO'})

    rdmol = molkit.to_rdmol(MOL)
    mol_dict = Settings({'mol': rdmol})
    _parse_name_type(mol_dict)
    assertion.eq(mol_dict, {'mol': rdmol, 'type': 'rdmol', 'name': 'CO'})

    mol_dict = Settings({'mol': PATH})
    _parse_name_type(mol_dict)
    assertion.eq(mol_dict, {'mol': PATH, 'type': 'folder', 'name': 'test_files'})

    mol_dict = Settings({'mol': MOL_PATH})
    _parse_name_type(mol_dict)
    assertion.eq(mol_dict, {'mol': MOL_PATH, 'type': 'xyz', 'name': 'Methanol'})

    mol_dict = Settings({'mol': 1})  # Excception: Invalid type
    assertion.assert_(_parse_name_type, mol_dict, exception=TypeError)


def test_validate_mol() -> None:
    """Test :func:`CAT.data_handling.validate_mol.validate_mol`."""
    args1 = ['Methanol.xyz', 'Ethylene.xyz']
    args2 = [
        Settings({'Acetate.xyz': {'guess_bonds': False}}),
        Settings({'Methanol_rotate.xyz': {'guess_bonds': False}})
    ]

    ref1 = [
        {
            'path': PATH,
            'is_core': True,
            'is_qd': False,
            'mol': join(PATH, 'Methanol.xyz'),
            'type': 'xyz',
            'name': 'Methanol',
            'parsed': True,
            'canonicalize': False,
        },
        {
            'path': PATH,
            'is_core': True,
            'is_qd': False,
            'mol': join(PATH, 'Ethylene.xyz'),
            'type': 'xyz',
            'name': 'Ethylene',
            'parsed': True,
            'canonicalize': False,
        },
    ]

    ref2 = [
        {
            'guess_bonds': False,
            'is_core': False,
            'is_qd': False,
            'path': PATH,
            'mol': join(PATH, 'Acetate.xyz'),
            'type': 'xyz',
            'name': 'Acetate',
            'parsed': True,
            'canonicalize': True,
        },
        {
            'guess_bonds': False,
            'is_core': False,
            'is_qd': False,
            'path': PATH,
            'mol': join(PATH, 'Methanol_rotate.xyz'),
            'type': 'xyz',
            'name': 'Methanol_rotate',
            'parsed': True,
            'canonicalize': True,
        }
    ]

    validate_mol(args1, 'input_cores', PATH)
    validate_mol(args2, 'input_ligands', PATH)
    assertion.eq(args1, ref1)
    assertion.eq(args2, ref2)
