"""Tests for the various workflows within the CAT package."""

import shutil
from pathlib import Path
from os.path import isdir

import yaml
import pytest
from scm.plams import Settings

from CAT.base import prep

PATH = Path('tests') / 'test_files'


@pytest.mark.slow
def test_CAT() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))
    arg.path = PATH

    try:
        qd_df, core_df, ligand_df = prep(arg)
    finally:
        lig_dir = PATH / 'ligand'
        qd_dir = PATH / 'qd'
        db_dir = PATH / 'database'

        shutil.rmtree(lig_dir) if isdir(lig_dir) else None
        shutil.rmtree(qd_dir) if isdir(qd_dir) else None
        shutil.rmtree(db_dir) if isdir(db_dir) else None
