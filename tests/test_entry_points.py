"""Tests for :mod:`CAT.data_handling.entry_points`."""

import os
import sys
import shutil
import warnings
from pathlib import Path
from typing import Generator

import rdkit
import pytest
from scm.plams import Molecule
from packaging.version import Version

from CAT.data_handling.entry_points import main
from CAT.test_utils import assert_mol_allclose

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
        with pytest.raises(FileNotFoundError), warnings.catch_warnings():
            warnings.filterwarnings(
                action="ignore",
                message="divide by zero encountered in double_scalars",
                category=RuntimeWarning,
            )
            main([f'{filename}bob'])

    @pytest.mark.parametrize("filename", PATH_DICT.values(), ids=PATH_DICT.keys())
    @pytest.mark.xfail(
        Version(rdkit.__version__) < Version("2021.03.4"),
        reason="requires rdkit >= 2021.03.4",
    )
    def test_mol(self, filename: str) -> None:
        mol = Molecule(LIG_PATH / filename)
        mol_ref = Molecule(PATH_REF / filename)

        xfail_linux = {
            "CCCCCCCCP[CCCCCCCC]CCCCCCCC@P25.xyz",
            "CCCCCCCCP[=O][CCCCCCCC]CCCCCCCC@O26.xyz",
        }
        xfail_windows = xfail_linux | {
            "O=C[COc1ccc[Cl]cc1Cl]Nc1cc[C[=O][O-]]ccc1F@O21.xyz",
        }
        xfail_unix = {
            "O=C[[O-]]c1ccc[S[=O][=O]CCCC[F][F]F]o1@O15.xyz",
            "O=C[COc1ccc[Cl]cc1Cl]Nc1cc[C[=O][O-]]ccc1F@O21.xyz",
        }
        if filename in xfail_linux and sys.platform.startswith("linux"):
            pytest.xfail("platform dependant geometry")
        elif filename in xfail_windows and sys.platform == "win32":
            pytest.xfail("platform dependant geometry")
        elif filename in xfail_unix and (
            sys.platform.startswith("linux") or sys.platform == "darwin"
        ):
            pytest.xfail("platform dependant geometry")

        # Coordinates
        assert_mol_allclose(mol, mol_ref, rtol=0, atol=10**-2)
