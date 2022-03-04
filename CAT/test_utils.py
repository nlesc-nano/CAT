"""Utility functions for the CAT tests.

Index
-----
.. currentmodule:: CAT.test_utils
.. autosummary::
    assert_mol_allclose

API
---
.. autofunction:: assert_mol_allclose

"""

from typing import TYPE_CHECKING, Any

import numpy as np
from scm.plams import Molecule

if TYPE_CHECKING:
    from numpy.typing import NDArray, ArrayLike
    _RecArray = np.recarray[Any, np.dtype[np.record]]

__all__ = ["assert_mol_allclose"]


def _is_within_tol(
    a: "NDArray[np.float64]",
    b: "NDArray[np.float64]",
    rtol: float,
    atol: float,
) -> "tuple[bool, int | np.intp, float | np.float64, float | np.float64]":
    """Check if ``a`` is equal (within a given tolerance) to any of the element in ``b``."""
    isclose = np.isclose(a, b, rtol, atol).all(axis=-1)
    if isclose.any():
        return True, isclose.nonzero()[0][0], 0.0, 0.0
    else:
        error = abs(a - b)
        i = np.nanargmin(error.sum(axis=-1))
        max_error = error[i].max()
        rel_error = (error[i] / abs(b[i])).max()
        return False, -1, max_error, rel_error


def _overlap_coordinates(
    xyz1: "NDArray[np.float64]",
    xyz2: "NDArray[np.float64]",
) -> "NDArray[np.float64]":
    """Remove all translations/rotations of ``xyz2`` w.r.t. ``xyz1``.

    Rotations are removed via a partial Procrustes superimposition.

    """
    ret = xyz2 - (xyz2.mean(axis=0) - xyz1.mean(axis=0))

    # Peform a singular value decomposition on the covariance matrix
    H = ret.T @ xyz1
    U, _, Vt = np.linalg.svd(H)
    V, Ut = Vt.T, U.T

    # Construct the rotation matrix
    rotmat = np.ones_like(U)
    rotmat[2, 2] = np.linalg.det(V @ Ut)
    rotmat *= V @ Ut

    # Apply the rotation matrix
    return np.matmul(ret, rotmat.T, out=ret)


def _compare_atoms(
    a: "NDArray[np.float64]",
    b: "NDArray[np.float64]",
    rtol: float,
    atol: float,
) -> "_RecArray":
    """Compare the coordinates of two Cartesian coordinate arrays."""
    b = _overlap_coordinates(a, b)

    dtype = np.dtype([
        ("isclose", "?"),
        ("idx", "u8"),
        ("max_err", "f8"),
        ("rel_err", "f8"),
    ])
    ret = np.zeros(len(a), dtype=dtype).view(np.recarray)

    # Iterate through all atoms and mask any successfully matched atoms with `nan`
    for i, item in enumerate(a):
        ret[i] = isclose, idx, *_ = _is_within_tol(item, b, rtol, atol)
        if isclose:
            b[idx] = np.nan
    return ret


def _compare_bonds(
    actual: Molecule,
    desired: Molecule,
    argsort_idx: "NDArray[np.integer]",
) -> "NDArray[np.bool_]":
    """Compare the bonds of two Molecules."""
    a_bonds = np.triu(actual.bond_matrix())[..., argsort_idx]
    b_bonds = np.triu(desired.bond_matrix())
    idx_tup = b_bonds.nonzero()
    return np.isclose(a_bonds, b_bonds, equal_nan=True)[idx_tup]


def _compare_atnums(
    actual: Molecule,
    desired: Molecule,
    argsort_idx: "NDArray[np.integer]",
) -> "NDArray[np.bool_]":
    """Compare the atomic numbers of two Molecules."""
    a = np.fromiter([at.atnum for at in actual], dtype=np.int64, count=len(actual))[argsort_idx]
    b = np.fromiter([at.atnum for at in desired], dtype=np.int64, count=len(actual))
    return a == b


def _compare_lattice(actual: Molecule, desired: Molecule, rtol: float, atol: float) -> np.bool_:
    """Compare the lattice vectors of two Molecules."""
    a_lat: "None | ArrayLike" = actual.lattice
    b_lat: "None | ArrayLike" = desired.lattice
    if a_lat is b_lat is None:
        return np.True_
    elif (a_lat is None and b_lat is not None) or (b_lat is None and a_lat is not None):
        return np.False_

    a = np.asarray(a_lat, dtype=np.float64)
    b = np.asarray(a_lat, dtype=np.float64)
    if a.shape != b.shape:
        return np.False_
    else:
        return np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True).all()


def assert_mol_allclose(
    actual: Molecule,
    desired: Molecule,
    *,
    rtol: float = 1e-05,
    atol: float = 1e-08,
    err_msg: str = "",
    verbose: bool = True,
) -> None:
    """Raises an AssertionError if two molecules are not equal up to desired tolerance.

    The test is equivalent to performing ``np.isclose(actual[i], desired, rtol, atol).any(axis=1)``
    for all atoms in ``actual``, making sure that each atoms is only allowed to match once.

    Parameters
    ----------
    actual : plams.Molecule
        Molecule obtained.
    desired : plams.Molecule
        Molecule desired.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.
    err_msg : str
        The error message to be printed in case of failure.
    verbose : bool
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
        If actual and desired are not equal up to specified precision.

    """
    __tracebackhide__ = True

    a = np.array(actual, dtype=np.float64)
    b = np.array(desired, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"Molecules must have the same number of atoms: {len(a)} {len(b)}")
    elif len(actual.bonds) != len(desired.bonds):
        raise ValueError(f"Molecules must have the same number of bonds: {len(actual.bonds)} {len(desired.bonds)}")  # noqa
    elif not (np.isfinite(a) & np.isfinite(b)).all():
        raise ValueError("Moleculair coordinates must not contain NaN and/or Inf")

    # Compare atomic coordinates, bonds and lattice vectors
    atoms_match = _compare_atoms(a, b, rtol=rtol, atol=atol)
    atnum_match = _compare_atnums(actual, desired, argsort_idx=atoms_match.idx)
    bonds_match = _compare_bonds(actual, desired, argsort_idx=atoms_match.idx)
    lattice_match = _compare_lattice(actual, desired, rtol=rtol, atol=atol)

    if atoms_match.isclose.all() & bonds_match.all() & lattice_match:
        return None

    at_mismatch = (~atoms_match.isclose).sum()
    at_percentage = round(100 * at_mismatch / len(a))
    atnum_mismatch = (~atnum_match).sum()
    atnum_percentage = round(100 * atnum_mismatch / len(a))
    bond_mismatch = (~bonds_match).sum()
    bond_percentage = round(100 * bond_mismatch / len(a))

    header = f"Not equal to tolerance rtol={rtol:g}, atol={atol:g}"
    err_msg += f"\nMismatched lattice:            {~lattice_match}"
    err_msg += f"\nMismatched bonds:              {bond_mismatch} / {len(a)} ({bond_percentage:.3g}%)"
    err_msg += f"\nMismatched atoms:              {at_mismatch} / {len(a)} ({at_percentage:.3g}%)"
    err_msg += f"\nMismatched atomic symbols:     {atnum_mismatch} / {len(a)} ({atnum_percentage:.3g}%)"
    err_msg += f"\nAtoms max absolute difference: {atoms_match.max_err.max()}"
    err_msg += f"\nAtoms max relative difference: {atoms_match.rel_err.max()}"
    raise AssertionError(np.testing.build_err_msg(
        [a, b], err_msg=err_msg, names=("actual", "desired"),
        verbose=verbose, header=header
    ))
