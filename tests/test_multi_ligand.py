"""Tests for :mod:`CAT.multi_ligand`."""

from pathlib import Path

import os
import yaml
from assertionlib import assertion

import pytest
import rdkit
from nanoutils import delete_finally
from scm.plams import Settings
from packaging.version import Version
from CAT import base

PATH = Path('tests') / 'test_files'
PDB_REF = PATH / 'multi_ligand.pdb'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'

with open(PATH / 'input3.yaml', 'r') as f:
    SETTINGS = Settings(yaml.load(f, Loader=yaml.FullLoader))


def _iterfiles(file1, file2):
    for i, j in zip(file1, file2):
        yield (
            i.rstrip(os.linesep).rstrip(),
            j.rstrip(os.linesep).rstrip()
        )


@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
@pytest.mark.xfail(
    Version(rdkit.__version__) < Version("2021.03.4"),
    reason="requires rdkit >= 2021.03.4",
)
def test_multi_ligand() -> None:
    """Test :mod:`CAT.multi_ligand`."""
    base.prep(SETTINGS.copy())
    pdb = PATH / 'qd' / 'Br25Cd360Cl25F25I27Se309__25_CCCC[O-]@O5__50_CCCCCCCCCCCCC[O-]@O14__77_C[O-]@O2.pdb'  # noqa
    with open(pdb, 'r') as f1, open(PDB_REF, 'r') as f2:
        iterator = enumerate(_iterfiles(f1, f2), 1)
        for i, (line1, line2) in iterator:
            assertion.eq(line1, line2, message=f'atom {i}')
