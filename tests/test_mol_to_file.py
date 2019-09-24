"""Tests for :mod:`CAT.data_handling.mol_to_file`."""

from os import mkdir
from os.path import (join, isdir)
from shutil import rmtree

from CAT.assertion.assertion_manager import assertion
from CAT.data_handling.mol_to_file import mol_to_file

from scm.plams import Molecule

_PATH = join('tests', 'test_files')
PATH = join(_PATH, 'mol_to_file')
MOL = Molecule(join(_PATH, 'Methanol.pdb'))


def test_mol_to_file() -> None:
    """Test :func:`CAT.data_handling.mol_to_file.mol_to_file`."""
    if not isdir(PATH):
        mkdir(PATH)

    _ref = ('mol.xyz', 'mol.pdb', 'mol.mol', 'mol.mol2')
    ref = [join(PATH, item) for item in _ref]

    mol = MOL.copy()
    mol.properties.name = 'mol'

    mol_list = [mol]
    mol_format = ('xyz', 'pdb', 'mol', 'mol2')
    kwargs = {'path': PATH, 'mol_format': mol_format}

    try:
        mol_to_file(mol_list, **kwargs)
        for file in ref:
            assertion.isfile(file)

        rmtree(PATH)
        mkdir(PATH)
        kwargs['mol_format'] = kwargs['mol_format'][0:2]
        mol_to_file(mol_list, **kwargs)
        for file in ref[:2]:
            assertion.isfile(file)
        for file in ref[2:]:
            assertion.isfile(file, invert=True)

        kwargs['path'] = join(PATH, 'bob')
        assertion.exception(FileNotFoundError, mol_to_file, mol_list, **kwargs)

        kwargs['path'] = join(PATH, 'mol.xyz')
        assertion.exception(NotADirectoryError, mol_to_file, mol_list, **kwargs)

        kwargs['path'] = PATH
        kwargs['mol_format'] = ('bob')
        assertion.exception(ValueError, mol_to_file, mol_list, **kwargs)

    finally:
        rmtree(PATH)
