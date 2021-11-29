"""A pytest ``conftest.py`` file."""

import logging
from typing import Any

import pytest
import rdkit

from nanoutils import VersionInfo


def pytest_configure(config: Any) -> None:
    logging.getLogger("filelock").setLevel(logging.WARNING)
    logging.getLogger("CAT").setLevel(logging.WARNING)


@pytest.fixture(scope="session", autouse=True)
def rdkit_version() -> None:
    rdkit_version = VersionInfo.from_str(rdkit.__version__)
    if rdkit_version < (2021, 3, 4):
        pytest.fail(
            "The CAT test suite requires RDKit >= 2021.03.4; "
            f"observed version: {rdkit.__version__}",
            pytrace=True,
        )
