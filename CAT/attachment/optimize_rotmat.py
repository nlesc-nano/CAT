"""A module for optimizing rottion matrices.

Index
-----
.. currentmodule:: CAT.attachment.optimize_rotmat
.. autosummary::
    optimize_rotmat
API
---
.. autofunction:: optimize_rotmat

"""

from __future__ import annotations

from typing import Callable, TYPE_CHECKING

import numpy as np
import scipy
from scipy.optimize import minimize
from packaging.version import Version

from scm.plams import rotation_matrix, Atom, Molecule

if Version(scipy.__version__) >= Version("1.1.0"):
    from scipy.optimize import Bounds
    BOUNDS = Bounds(-1, 1, keep_feasible=True)
else:
    BOUNDS = None

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from numpy import float64 as f8

__all__ = ['optimize_rotmat']


@np.errstate(invalid="ignore")
def _get_angle(
    xyz: "NDArray[f8]",
    anchor_idx: int,
    ignore_idx: "None | NDArray[np.integer]",
) -> "NDArray[f8] | f8":
    """Return the maximum angle in ``xyz`` w.r.t. to the X-axis."""
    vecs = xyz / np.linalg.norm(xyz, axis=-1)[..., None]
    if ignore_idx is not None:
        vecs[ignore_idx] = [1, 0, 0]
    vecs[anchor_idx] = [1, 0, 0]
    angles = np.arccos(vecs @ [1, 0, 0])
    return np.nanmax(angles, axis=-1)


def _minimize_func2(
    vec1: "NDArray[f8]",
    vec2: "NDArray[f8]",
    xyz: "NDArray[f8]",
    anchor: int,
    ignore_idx: "None | NDArray[np.integer]",
) -> f8:
    """Rotate the X-axis in ``xyz`` to ``vec`` and \
    compute the maximum angle w.r.t. to the X-axis."""
    rotmat = rotation_matrix(vec1, vec2)
    xyz_rot = xyz @ rotmat.T
    xyz_rot -= xyz_rot[anchor]
    return _get_angle(xyz_rot, anchor, ignore_idx)


def _minimize_func(vec1: np.ndarray, vec2: np.ndarray, xyz: np.ndarray, anchor: int) -> float:
    """The function whose output is to-be minimized by :func:`scipy.optimize.minimize`."""
    # Rotate and translate
    rotmat = rotation_matrix(vec1, vec2)
    xyz_new = xyz@rotmat
    xyz_new -= xyz_new[anchor]

    # Apply the cost function: e^(||yz||)
    yz = xyz_new[:, 1:]
    distance = np.linalg.norm(yz, axis=1)
    return np.exp(distance).sum()


def optimize_rotmat(
    mol: "NDArray[np.float64] | Molecule",
    anchor: int = 0,
    ignore_idx: "None | NDArray[np.integer]" = None,
    func: Callable = _minimize_func,
) -> np.ndarray:
    r"""Find the rotation matrix for **xyz** that minimizes its deviation from the Cartesian X-axis.

    A set of vectors, :math:`v`, is constructed for all atoms in **xyz**,
    denoting their deviation from the Cartesian X-axis.
    Subsequently, the rotation matrix that minimizes :math:`\sum_{i}^{n} {e^{||v_{i}||}}`
    is returned.

    Parameters
    ----------
    mol : :class:`numpy.ndarray` [:class:`float`], shape :math:`(n,3)`
        An array-like object of Cartesian coordinates.
        *e.g.* :class:`list`, :class:`numpy.ndarray` or |plams.Molecule|.

    anchor : |plams.Atom| or :class:`int`
        The index (0-based) of the anchor atom in **xyz**.
        Used for defining the origin in **xyz** (*i.e.* :code:`xyz[i] == (0, 0, 0)`).
        Alternativelly, a PLAMS atom can be passed.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`float`], shape :math:`(3,3)` or :math:`(3,)`
        An optimized rotation matrix.
        If ``as_vec=True``, instead return the initial vector
        used for constructing the rotation matrix.

    See Also
    --------
    :func:`minimize()<scipy.optimize.minimize>`
        Minimization of scalar function of one or more variables.

    """
    if hasattr(anchor, '__int__'):  # This encompasses both int and np.integer instances
        i = int(anchor)
    elif isinstance(anchor, Atom):
        i = anchor.mol.atoms.index(anchor)
    else:
        raise TypeError("The passed anchor is neither an 'int' nor 'Atom'; "
                        f"observed type: {repr(type(anchor))}")
    xyz = np.array(mol, dtype=float, ndmin=2, copy=False)

    # Create a first guess for the starting
    vec1 = xyz.mean(axis=0) - xyz[i]
    vec2 = np.array([1, 0, 0], dtype=float)
    with np.errstate(invalid='raise'):
        try:
            rotmat1 = rotation_matrix(vec1, vec2)
        except FloatingPointError:
            return np.eye(3)
    xyz_new = xyz@rotmat1.T

    # Optimize the trial vector; return the matching rotation matrix
    trial_vec = xyz_new.mean(axis=0) - xyz_new[i]
    output = minimize(
        _minimize_func2,
        trial_vec,
        args=(vec2, xyz_new, i, ignore_idx),
        bounds=BOUNDS,
    )
    rotmat2 = rotation_matrix(output.x, vec2)
    return (rotmat1.T@rotmat2.T).T
