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

from typing import TYPE_CHECKING

import numpy as np
from scm.plams import Molecule

if TYPE_CHECKING:
    from numpy.typing import NDArray


def _is_within_tol(
    a: "NDArray[np.float64]",
    b: "NDArray[np.float64]",
    rtol: float,
    atol: float,
) -> "tuple[bool, int | np.intp, float | np.float64, float | np.float64]":
    """Check if ``a`` is equal (within a given tolerance) to any of the element in ``b``."""
    nan_mask = ~np.isnan(b)
    isclose = abs(a - b) <= (atol + rtol * abs(b))
    ret = np.zeros(isclose.shape[:-1], dtype=np.bool_)
    ret = np.logical_and.reduce(isclose, where=nan_mask, out=ret, axis=-1)
    if ret.any():
        return True, np.where(ret)[0][0], 0.0, 0.0
    else:
        error = abs(a - b)
        i = np.nanargmin(error.sum(axis=-1))
        max_error = error[i].max()
        rel_error = (error[i] / abs(b[i])).max()
        return False, -1, max_error, rel_error


def _overlap_coordinates(xyz1: "NDArray[np.float64]", xyz2: "NDArray[np.float64]") -> None:
    """Remove all translations/rotations of ``xyz2`` w.r.t. ``xyz1``.

    Rotations are removed via a partial Procrustes superimposition.

    Performs an inplace update of the coordinates in ``xyz2``.

    """
    xyz2 -= (xyz2.mean(axis=0) - xyz1.mean(axis=0))

    # Peform a singular value decomposition on the covariance matrix
    H = xyz2.T @ xyz1
    U, _, Vt = np.linalg.svd(H)
    V, Ut = Vt.T, U.T

    # Construct the rotation matrix
    rotmat = np.ones_like(U)
    rotmat[2, 2] = np.linalg.det(V @ Ut)
    rotmat *= V @ Ut

    # Apply the rotation matrix inplace
    np.matmul(xyz2, rotmat.T, out=xyz2)


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
        raise ValueError(f"Operands should have the same shapes: {a.shape} {b.shape}")
    elif not (np.isfinite(a).all() & np.isfinite(b).all()):
        raise NotImplementedError
    _overlap_coordinates(a, b)

    dtype = np.dtype([("isclose", "?"), ("max_err", "f8"), ("rel_err", "f8")])
    ret = np.zeros(len(a), dtype=dtype).view(np.recarray)
    for i, item in enumerate(a):
        isclose, idx, max_err, rel_err = _is_within_tol(item, b, rtol, atol)
        if isclose:
            b[idx] = np.nan
        ret[i] = (isclose, max_err, rel_err)

    if ret.isclose.all():
        return None

    n_mismatch = (~ret.isclose).sum()
    percentage = round(100 * n_mismatch / len(b))

    header = f"Not equal to tolerance rtol={rtol:g}, atol={atol:g}"
    err_msg += f"\nMismatched atoms: {n_mismatch} / {len(b)} ({percentage:.3g}%)"
    err_msg += f"\nMax absolute difference: {ret.max_err.max()}"
    err_msg += f"\nMax relative difference: {ret.rel_err.max()}"
    raise AssertionError(np.testing.build_err_msg(
        [a, b], err_msg=err_msg, names=("actual", "desired"),
        verbose=verbose, header=header
    ))
