"""Tests for :mod:`CAT.attachment.ligand_attach`."""

import sys
import shutil
from pathlib import Path
from typing import Generator, NamedTuple, TYPE_CHECKING

import yaml
import pytest
import numpy as np
from assertionlib import assertion
from scm.plams import Settings, Molecule

from CAT.base import prep
from CAT.attachment.ligand_attach import (_get_rotmat1, _get_rotmat2)
from CAT.workflows import MOL

if TYPE_CHECKING:
    import _pytest

PATH = Path('tests') / 'test_files'

LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


def test_get_rotmat1() -> None:
    """Test :func:`CAT.attachment.ligand_attach._get_rotmat1`."""
    vec1 = np.array([1, 0, 0], dtype=float)
    vec2 = np.array([0, 1, 0], dtype=float)
    vec3 = np.array([[0, 1, 0], [0, 0, 1]], dtype=float)
    vec4 = np.zeros(3, dtype=float)

    ref1 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], ndmin=3, dtype=float)
    _ref2 = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]], ndmin=3, dtype=float)
    ref2 = np.vstack([ref1, _ref2])
    ref3 = np.zeros((2, 3, 3))
    ref3[:] = np.identity(3)
    ref4 = np.array([[[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                     [[0, 0, -1], [0, 1, 0], [1, 0, 0]]], dtype=float)
    ref5 = np.identity(3)[None, ...]

    rotmat1 = _get_rotmat1(vec1, vec2)
    np.testing.assert_allclose(rotmat1, ref1)

    rotmat2 = _get_rotmat1(vec1, vec3)
    np.testing.assert_allclose(rotmat2, ref2)

    rotmat3 = _get_rotmat1(vec3, vec3)
    np.testing.assert_allclose(rotmat3, ref3)

    rotmat4 = _get_rotmat1(vec3, vec1)
    np.testing.assert_allclose(rotmat4, ref4)

    rotmat5 = _get_rotmat1(vec4, vec1)
    np.testing.assert_allclose(rotmat5, ref5)

    assertion.assert_(_get_rotmat1, vec3[..., None], vec3, exception=ValueError)


def test_get_rotmat2() -> None:
    """Test :func:`CAT.attachment.ligand_attach._get_rotmat2`."""
    vec1 = np.array([1, 0, 0], dtype=float)
    vec2 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)

    ref1 = np.load(PATH / 'rotmat2_1.npy')
    ref2 = np.load(PATH / 'rotmat2_2.npy')

    rotmat1 = _get_rotmat2(vec1)
    np.testing.assert_allclose(rotmat1, ref1)

    rotmat2 = _get_rotmat2(vec2)
    np.testing.assert_allclose(rotmat2, ref2)


class DihedTup(NamedTuple):
    mol: Molecule
    ref: np.recarray


class TestDihedral:
    PARAMS = {
        "dihed_45": 45,
        "dihed_45_deg": "45 deg",
        "dihed_180": 180.0,
    }

    @pytest.fixture(scope="class", autouse=True, name="output", params=PARAMS.items(), ids=PARAMS)
    def run_cat(self, request: "_pytest.fixtures.SubRequest") -> Generator[DihedTup, None, None]:
        # Setup
        name, dihed = request.param  # type: str, str | float
        yaml_path = PATH / 'CAT_dihedral.yaml'
        with open(yaml_path, 'r') as f:
            arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

        arg.path = PATH
        arg.optional.ligand.anchor.dihedral = dihed
        qd_df, _, _ = prep(arg)
        qd = qd_df[MOL].iloc[0]

        ref = np.load(PATH / f"test_dihedral_{name}.npy").view(np.recarray)
        yield DihedTup(qd, ref)

        # Teardown
        files = [LIG_PATH, QD_PATH, DB_PATH]
        for f in files:
            shutil.rmtree(f, ignore_errors=True)

    def test_atoms(self, output: DihedTup) -> None:
        dtype = [("symbols", "U2"), ("coords", "f8", 3)]
        atoms = np.fromiter(
            [(at.symbol, at.coords) for at in output.mol], dtype=dtype
        ).view(np.recarray)

        assertion.eq(atoms.dtype, output.ref.dtype)
        np.testing.assert_array_equal(atoms.symbols, output.ref.symbols)

        if sys.version_info >= (3, 9):
            pytest.xfail("Geometries must be updated for RDKit >2019.09.2")
        np.testing.assert_allclose(atoms.coords, output.ref.coords)
