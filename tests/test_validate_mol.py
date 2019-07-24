"""Tests for :mod:`CAT.data_handling.validate_mol`."""

from os.path import join

from scm.plams import (Settings, Molecule)
import scm.plams.interfaces.molecule.rdkit as molkit

from CAT.assertion_functions import (assert_eq, assert_exception)
from CAT.data_handling.validate_mol import (
    validate_mol, santize_smiles, _parse_name_type, _parse_mol_type
)

PATH = join('tests', 'test_files')
MOL_PATH = join(PATH, 'Methanol.xyz')
MOL = Molecule(MOL_PATH)
MOL.guess_bonds()


def test_santize_smiles() -> None:
    """Test :func:`CAT.data_handling.validate_mol.santize_smiles`."""
    assert_eq(santize_smiles('CO'), 'CO')
    assert_eq(santize_smiles('C[H]O'), 'C[H]O')
    assert_eq(santize_smiles('C(C)O'), 'C[C]O')
    assert_eq(santize_smiles('CC=CC'), 'CC=CC')
    assert_eq(santize_smiles('C/C=C/C'), 'trans-CC=CC')
    assert_eq(santize_smiles(r'C/C=C\C'), 'cis-CC=CC')
    assert_eq(santize_smiles('C/C=C/C/C=C/C'), 'trans-trans-CC=CCC=CC')
    assert_eq(santize_smiles(r'C/C=C\C/C=C/C'), 'cis-trans-CC=CCC=CC')


def test_parse_mol_type() -> None:
    """Test :func:`CAT.data_handling.validate_mol._parse_mol_type`."""
    assert_eq(_parse_mol_type('input_cores'), True)
    assert_eq(_parse_mol_type('input_ligands'), False)
    assert_exception(ValueError, _parse_mol_type, 'bob')
    assert_exception(AttributeError, _parse_mol_type, 1)


def test_parse_name_type() -> None:
    """Test :func:`CAT.data_handling.validate_mol._parse_name_type`."""
    mol_dict = Settings({'mol': 'CCCO'})
    _parse_name_type(mol_dict)
    assert_eq(mol_dict, {'mol': 'CCCO', 'type': 'smiles', 'name': 'CCCO'})

    mol_dict = Settings({'mol': MOL})
    _parse_name_type(mol_dict)
    assert_eq(mol_dict, {'mol': MOL, 'type': 'plams_mol', 'name': 'Methanol'})

    mol = MOL.copy()
    mol.properties = Settings()
    mol_dict = Settings({'mol': mol})
    _parse_name_type(mol_dict)
    assert_eq(mol_dict, {'mol': mol, 'type': 'plams_mol', 'name': 'CO'})

    rdmol = molkit.to_rdmol(MOL)
    mol_dict = Settings({'mol': rdmol})
    _parse_name_type(mol_dict)
    assert_eq(mol_dict, {'mol': rdmol, 'type': 'rdmol', 'name': 'CO'})

    mol_dict = Settings({'mol': PATH})
    _parse_name_type(mol_dict)
    assert_eq(mol_dict, {'mol': PATH, 'type': 'folder', 'name': 'test_files'})

    mol_dict = Settings({'mol': MOL_PATH})
    _parse_name_type(mol_dict)
    assert_eq(mol_dict, {'mol': MOL_PATH, 'type': 'xyz', 'name': 'Methanol'})

    mol_dict = Settings({'mol': 1})  # Excception: Invalid type
    assert_exception(TypeError, _parse_name_type, mol_dict)


def test_validate_mol() -> None:
    """Test :func:`CAT.data_handling.validate_mol.validate_mol`."""
    args1 = ['Methanol.xyz', 'Ethylene.xyz']
    args2 = [
        Settings({'Acetate.xyz': {'guess_bonds': False}}),
        Settings({'Methanol_rotate.xyz': {'guess_bonds': False}})
    ]

    ref1 = [
        {'path': PATH,
         'is_core': True,
         'mol': join(PATH, 'Methanol.xyz'),
         'type': 'xyz',
         'name': 'Methanol'},
        {'path': PATH,
         'is_core': True,
         'mol': join(PATH, 'Ethylene.xyz'),
         'type': 'xyz',
         'name': 'Ethylene'}
    ]

    ref2 = [
        {'guess_bonds': False,
         'is_core': False,
         'path': PATH,
         'mol': join(PATH, 'Acetate.xyz'),
         'type': 'xyz',
         'name': 'Acetate'},
        {'guess_bonds': False,
         'is_core': False,
         'path': PATH,
         'mol': join(PATH, 'Methanol_rotate.xyz'),
         'type': 'xyz',
         'name': 'Methanol_rotate'}
    ]

    validate_mol(args1, 'input_cores', PATH)
    validate_mol(args2, 'input_ligands', PATH)
    assert_eq(args1, ref1)
    assert_eq(args2, ref2)
