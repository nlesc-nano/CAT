"""A pytest ``conftest.py`` file."""

import pytest
import rdkit

import warnings
from nanoutils import VersionInfo


@pytest.fixture(scope="session", autouse=True)
def rdkit_version() -> None:
    rdkit_version = VersionInfo.from_str(rdkit.__version__)
    if rdkit_version < (2021, 3, 4):
        warnings.warn(
            "The CAT test suite requires RDKit >= 2021.03.4 to properly run; "
            f"observed version: {rdkit.__version__}"
        )
