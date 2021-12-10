"""A pytest ``conftest.py`` file."""

import re
import logging
import warnings

import pytest
import rdkit
from nanoutils import VersionInfo

from CAT import logger


@pytest.fixture(scope="session", autouse=True)
def rdkit_version() -> None:
    """Validate the rdkit version."""
    rdkit_version = VersionInfo.from_str(rdkit.__version__)
    if rdkit_version < (2021, 3, 4):
        warnings.warn(
            "The CAT test suite requires RDKit >= 2021.03.4 to properly run; "
            f"observed version: {rdkit.__version__}"
        )


class TestFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        if re.fullmatch(r"The optional [\-a-zA-Z]+ package was not found", record.msg):
            return False
        return True


@pytest.fixture(scope="session", autouse=True)
def modify_logger() -> None:
    """Suppres the loggers ``stderr`` output and filter out all optional-package warnings."""
    for handler in logger.handlers.copy():
        logger.removeHandler(handler)

    filter = TestFilter(logger.name)
    logger.addFilter(filter)
