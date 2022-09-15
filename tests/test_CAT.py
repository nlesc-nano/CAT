"""Tests for the various workflows within the CAT package."""

from pathlib import Path
from typing import Optional

import yaml
import pytest
from scm.plams import Settings
from nanoutils import delete_finally

from CAT.base import prep

try:
    import nanoCAT  # noqa: F401
    NANOCAT_EX: Optional[ImportError] = None
except ImportError as ex:
    NANOCAT_EX = ex

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
@pytest.mark.skipif(NANOCAT_EX is not None, reason=str(NANOCAT_EX))
def test_cat_fail() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

    arg.path = PATH
    del arg.optional.core.anchor
    arg.optional.ligand["cosmo-rs"] = False
    arg.optional.qd.activation_strain = False
    arg.optional.qd.optimize = False
    del arg.input_cores[0]["Cd68Se55.xyz"]
    arg.input_cores[0]["CdCl.xyz"] = {"indices": [2]}
    with pytest.raises(ValueError):
        prep(arg)


@pytest.mark.slow
@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
@pytest.mark.skipif(NANOCAT_EX is not None, reason=str(NANOCAT_EX))
def test_cat() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

    arg.path = PATH
    prep(arg)
