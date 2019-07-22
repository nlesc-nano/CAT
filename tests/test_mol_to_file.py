"""Tests for :mod:`CAT.data_handling.mol_to_file`."""

from os import mkdir
from os.path import join
from shutil import rmtree

from CAT.assertion_functions import (assert_isfile, assert_exception, Invert)
from CAT.data_handling.mol_to_file import mol_to_file

from scm.plams import Molecule

_PATH = 'tests/test_files'
PATH = join(_PATH, 'mol_to_file')
MOL = Molecule(join(_PATH, 'Methanol.pdb'))


def test_mol_to_file() -> None:
    """Test :func:`CAT.data_handling.mol_to_file.mol_to_file`."""
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
            assert_isfile(file)

        rmtree(PATH)
        mkdir(PATH)
        kwargs['mol_format'] = kwargs['mol_format'][0:2]
        mol_to_file(mol_list, **kwargs)
        for file in ref[:2]:
            assert_isfile(file)
        for file in ref[2:]:
            with Invert(assert_isfile) as func:
                func(file)

        kwargs['path'] = join(PATH, 'bob')
        assert_exception(FileNotFoundError, mol_to_file, mol_list, **kwargs)

        kwargs['path'] = join(PATH, 'mol.xyz')
        assert_exception(NotADirectoryError, mol_to_file, mol_list, **kwargs)

        kwargs['path'] = PATH
        kwargs['mol_format'] = ('bob')
        assert_exception(ValueError, mol_to_file, mol_list, **kwargs)

    finally:
        rmtree(PATH)
        mkdir(PATH)
