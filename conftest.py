"""A pytest ``conftest.py`` file."""

import pytest
import rdkit

from nanoutils import VersionInfo


@pytest.fixture(scope="session", autouse=True)
def rdkit_version() -> None:
    rdkit_version = VersionInfo.from_str(rdkit.__version__)
    if rdkit_version < (2021, 3, 4):
        pytest.fail(
            "The CAT test suite requires RDKit >= 2021.03.4; "
            f"observed version: {rdkit.__version__}",
            pytrace=True,
        )
