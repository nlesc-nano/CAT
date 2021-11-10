"""Tests for :mod:`CAT.data_handling.entry_points`."""

import os
import sys
import shutil
import warnings
from pathlib import Path
from typing import Generator

import pytest
import numpy as np

from scm.plams import Molecule

from CAT.data_handling.warn_map import MoleculeWarning
from CAT.data_handling.entry_points import main

if sys.version_info >= (3, 7):
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict

PATH = Path('tests') / 'test_files'

LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'
PATH_REF = PATH / 'ligand_ref'


class TestMain:
    @pytest.fixture(scope="class", autouse=True)
    def run_cat(self) -> Generator[None, None, None]:
        # Setup
        filename = str(PATH / 'input2.yaml')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MoleculeWarning)
            main([filename])

        yield None

        # Teardown
        files = [LIG_PATH, QD_PATH, DB_PATH]
        for f in files:
            if os.path.isdir(f):
                shutil.rmtree(f)

    PATH_DICT = OrderedDict({
        os.path.splitext(f)[0]: f for f in os.listdir(PATH_REF) if f.endswith('xyz')
    })

    def test_raise(self) -> None:
        filename = str(PATH / 'input2.yaml')
        with pytest.raises(FileNotFoundError):
            main([f'{filename}bob'])

    @pytest.mark.parametrize("filename", PATH_DICT.values(), ids=PATH_DICT.keys())
    @pytest.mark.xfail(
        sys.version_info >= (3, 9),
        reason="Geometries must be updated for RDKit >2019.09.2",
    )
    def test_mol(self, filename: str) -> None:
        # Coordinates
        mol = Molecule(LIG_PATH / filename)
        mol_ref = Molecule(PATH_REF / filename)
        np.testing.assert_allclose(mol, mol_ref, rtol=0, atol=10**-2)

        # Atomic symbols
        symbols = [at.symbol for at in mol]
        symbol_ref = [at.symbol for at in mol_ref]
        np.testing.assert_array_equal(symbols, symbol_ref)
