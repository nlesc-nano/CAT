"""Tests for the various workflows within the CAT package."""

from pathlib import Path

import yaml
import pytest
from scm.plams import Settings
from nanoutils import delete_finally

from CAT.base import prep

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@pytest.mark.slow
@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_cat() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

    arg.path = PATH
    prep(arg)
