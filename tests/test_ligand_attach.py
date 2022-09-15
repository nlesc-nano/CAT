"""Tests for :mod:`CAT.attachment.ligand_attach`."""

import shutil
from pathlib import Path
from typing import Generator, NamedTuple, TYPE_CHECKING, Any

import rdkit
import yaml
import pytest
import h5py
import numpy as np
from assertionlib import assertion
from packaging.version import Version
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
    name: str


class AllignmentTup(NamedTuple):
    mol: Molecule
    atoms_ref: np.recarray
    bonds_ref: np.recarray
    name: str


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
        yield DihedTup(qd, ref, name)

        # Teardown
        files = [LIG_PATH, QD_PATH, DB_PATH]
        for file in files:
            shutil.rmtree(file, ignore_errors=True)

    @pytest.mark.xfail(
        Version(rdkit.__version__) < Version("2021.03.4"),
        reason="requires rdkit >= 2021.03.4",
    )
    def test_atoms(self, output: DihedTup) -> None:
        dtype = [("symbols", "U2"), ("coords", "f8", 3)]
        atoms = np.fromiter(
            [(at.symbol, at.coords) for at in output.mol], dtype=dtype
        ).view(np.recarray)

        assertion.eq(atoms.dtype, output.ref.dtype)
        np.testing.assert_array_equal(atoms.symbols, output.ref.symbols)
        np.testing.assert_allclose(atoms.coords, output.ref.coords)


class TestAllignment:
    PARAMS = ("sphere", "surface", "sphere_invert", "surface_invert")

    @pytest.fixture(scope="class", name="output", params=PARAMS)
    def run_cat(
        self, request: "_pytest.fixtures.SubRequest"
    ) -> Generator[AllignmentTup, None, None]:
        # Setup
        allignment: str = request.param
        yaml_path = PATH / 'CAT_allignment.yaml'
        with open(yaml_path, 'r') as f1:
            arg = Settings(yaml.load(f1, Loader=yaml.FullLoader))

        arg.path = PATH
        arg.optional.core.allignment = allignment
        qd_df, _, _ = prep(arg)
        qd = qd_df[MOL].iloc[0]

        with h5py.File(PATH / "test_allignment.hdf5", "r") as f2:
            atoms_ref = f2[f"TestAllignment/{allignment}/atoms"][...].view(np.recarray)
            bonds_ref = f2[f"TestAllignment/{allignment}/bonds"][...].view(np.recarray)
        yield AllignmentTup(qd, atoms_ref, bonds_ref, allignment)

        # Teardown
        files = [LIG_PATH, QD_PATH, DB_PATH]
        for file in files:
            shutil.rmtree(file, ignore_errors=True)

    @pytest.mark.xfail(
        Version(rdkit.__version__) < Version("2021.03.4"),
        reason="requires rdkit >= 2021.03.4",
    )
    def test_atoms(self, output: AllignmentTup) -> None:
        dtype = [("symbols", "S2"), ("coords", "f8", 3)]
        iterator = ((at.symbol, at.coords) for at in output.mol)
        atoms = np.fromiter(iterator, dtype=dtype).view(np.recarray)

        assertion.eq(atoms.dtype, output.atoms_ref.dtype)
        np.testing.assert_array_equal(atoms.symbols, output.atoms_ref.symbols)
        np.testing.assert_allclose(atoms.coords, output.atoms_ref.coords)

    def test_bonds(self, output: AllignmentTup) -> None:
        dtype = [("atom1", "i8"), ("atom2", "i8"), ("order", "f8")]
        try:
            output.mol.set_atoms_id()
            iterator = ((b.atom1.id, b.atom2.id, b.order) for b in output.mol.bonds)
            bonds = np.fromiter(iterator, dtype=dtype).view(np.recarray)
        finally:
            output.mol.unset_atoms_id()

        assertion.eq(bonds.dtype, output.bonds_ref.dtype)
        np.testing.assert_array_equal(bonds.atom1, output.bonds_ref.atom1)
        np.testing.assert_array_equal(bonds.atom2, output.bonds_ref.atom2)
        np.testing.assert_allclose(bonds.order, output.bonds_ref.order)


class TestCoreAnchor:
    PARAMS = {
        "HCl": ("Cd68Se55_HCl.pdb", {"group": "[H]Cl", "group_idx": 0, "remove": 0}),
        "formate": ("Cd68Cl26Se55__26_O=C[O-]@O2O3.pdb", {
            "group": "[H]C([O-])=O",
            "group_idx": [2, 3],
            "remove": [0, 1, 2, 3],
            "kind": "mean",
        }),
        "ethoxide": ("Cd68Cl26Se55__26_CC[O-]@O3.pdb", {
            "group": "[O-]C([H])([H])C([H])([H])[H]",
            "group_idx": 0,
            "remove": [0, 1, 2, 3, 4, 5, 6, 7],
            "kind": "first",
        }),
    }

    @pytest.fixture(scope="class", name="output", params=PARAMS.items(), ids=PARAMS)
    def run_cat(
        self, request: "_pytest.fixtures.SubRequest"
    ) -> Generator[AllignmentTup, None, None]:
        # Setup
        name, (core, kwargs) = request.param  # type: str, tuple[str, dict[str, Any]]
        yaml_path = PATH / 'CAT_allignment.yaml'
        with open(yaml_path, 'r') as f1:
            arg = Settings(yaml.load(f1, Loader=yaml.FullLoader))

        arg.path = PATH
        arg.input_cores = [core]
        arg.optional.core.anchor = kwargs
        if name == "formate":
            arg.optional.core.allignment = "anchor"
        qd_df, _, _ = prep(arg)
        qd = qd_df[MOL].iloc[0]

        with h5py.File(PATH / "test_allignment.hdf5", "r") as f2:
            atoms_ref = f2[f"TestCoreAnchor/{name}/atoms"][...].view(np.recarray)
            bonds_ref = f2[f"TestCoreAnchor/{name}/bonds"][...].view(np.recarray)
        yield AllignmentTup(qd, atoms_ref, bonds_ref, name)

        # Teardown
        files = [LIG_PATH, QD_PATH, DB_PATH]
        for file in files:
            shutil.rmtree(file, ignore_errors=True)

    @pytest.mark.xfail(
        Version(rdkit.__version__) < Version("2021.03.4"),
        reason="requires rdkit >= 2021.03.4",
    )
    def test_atoms(self, output: AllignmentTup) -> None:
        dtype = [("symbols", "S2"), ("coords", "f8", 3)]
        iterator = ((at.symbol, at.coords) for at in output.mol)
        atoms = np.fromiter(iterator, dtype=dtype).view(np.recarray)

        assertion.eq(atoms.dtype, output.atoms_ref.dtype)
        np.testing.assert_array_equal(atoms.symbols, output.atoms_ref.symbols)
        np.testing.assert_allclose(atoms.coords, output.atoms_ref.coords)

    def test_bonds(self, output: AllignmentTup) -> None:
        dtype = [("atom1", "i8"), ("atom2", "i8"), ("order", "f8")]
        try:
            output.mol.set_atoms_id()
            iterator = ((b.atom1.id, b.atom2.id, b.order) for b in output.mol.bonds)
            bonds = np.fromiter(iterator, dtype=dtype).view(np.recarray)
        finally:
            output.mol.unset_atoms_id()

        assertion.eq(bonds.dtype, output.bonds_ref.dtype)
        np.testing.assert_array_equal(bonds.atom1, output.bonds_ref.atom1)
        np.testing.assert_array_equal(bonds.atom2, output.bonds_ref.atom2)
        np.testing.assert_allclose(bonds.order, output.bonds_ref.order)
