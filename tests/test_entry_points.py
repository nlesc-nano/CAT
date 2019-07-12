"""Tests for :mod:`CAT.data_handling.entry_points`."""

from os.path import join
from shutil import rmtree

from CAT.assertion_functions import assert_exception
from CAT.data_handling.entry_points import main

PATH = 'tests/test_files'


def test_main() -> None:
    """Test :func:`CAT.data_handling.entry_points.main`."""
    filename = join(PATH, 'input2.yaml')
    try:
        main([filename])
    finally:
        rmtree(join(PATH, 'ligand'))
        rmtree(join(PATH, 'qd'))
        rmtree(join(PATH, 'database'))

    assert_exception(FileNotFoundError, main, [filename + 'bob'])
