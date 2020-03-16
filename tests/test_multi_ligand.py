"""Tests for :mod:`CAT.multi_ligand`."""

import os
import shutil
from pathlib import Path

import yaml
from assertionlib import assertion

from scm.plams import Settings
from CAT import base

PATH = Path('tests') / 'test_files'
PDB_REF = PATH / 'mutli_ligand.pdb'

with open(PATH / 'input3.yaml', 'r') as f:
    SETTINGS = Settings(yaml.load(f, Loader=yaml.FullLoader))


def test_multi_ligand() -> None:
    """Test :mod:`CAT.multi_ligand`."""
    try:
        base.prep(SETTINGS.copy())
        pdb = PATH / 'qd' / 'Br25Cd360Cl25F25I27Se309__25_CCCC[O-]@O5__50_CCCCCCCCCCCCC[O-]@O14__77_C[O-]@O2.pdb'  # noqa
        with open(pdb, 'r') as f1, open(PDB_REF, 'r') as f2:
            for i, j in zip(f1, f2):
                assertion.eq(i, j)
    finally:
        shutil.rmtree(PATH / 'ligand') if os.path.isdir(PATH / 'ligand') else None
        shutil.rmtree(PATH / 'qd') if os.path.isdir(PATH / 'qd') else None
        shutil.rmtree(PATH / 'database') if os.path.isdir(PATH / 'database') else None
