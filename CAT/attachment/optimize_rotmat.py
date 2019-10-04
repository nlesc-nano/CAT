"""
nanoCAT.mol_bulk
================

A module for calculating the bulkiness of ligands.

Index
-----
.. currentmodule:: nanoCAT.mol_bulk
.. autosummary::
    init_lig_bulkiness
    get_cone_angles
    get_V
    _export_to_db
    _get_anchor_idx

API
---
.. autofunction:: init_lig_bulkiness
.. autofunction:: get_cone_angles
.. autofunction:: get_V
.. autofunction:: _export_to_db
.. autofunction:: _get_anchor_idx

"""

import numpy as np
from scipy.optimize import minimize

from scm.plams import rotation_matrix, Atom

__all__ = ['optimize_rotmat']


def optimize_rotmat(mol: np.ndarray, anchor: int = 0) -> np.ndarray:
    r"""Find the rotation matrix for **xyz** that minimizes its deviation from the Cartesian X-axis.

    A set of vectors, :math:`v`, is constructed for all atoms in **xyz**,
    denoting their deviation from the Cartesian X-axis.
    Subsequently, the rotation matrix that minimizes :math:`\sum_{i}^{n} {e^{||v_{i}||}}`
    is returned.

    Parameters
    ----------
    mol : :math:`n*3` :class:`numpy.ndarray` [:class:`float`]
        An array-like object of Cartesian coordinates.
        *e.g.* :class:`list`, :class:numpy.ndarray: or |plams.Molecule|.

    anchor : |plams.Atom| or :class:`int`
        The index (0-based) of the anchor atom in **xyz**.
        Used for defining the origin in **xyz** (*i.e.* :code:`xyz[i] == (0, 0, 0)`).
        Alternativelly, a PLAMS atom can be passed.

    Returns
    -------
    :math:`3*3` or :math:`3` :class:`numpy.ndarray` [:class:`float`]
        An optimized rotation matrix.
        If ``as_vec=True``, instead return the initial vector
        used for constructing the rotation matrix.

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
    rotmat1 = rotation_matrix(vec1, vec2)
    xyz_new = xyz@rotmat1.T

    # Optimize the trial vector; return the matching rotation matrix
    trial_vec = xyz_new.mean(axis=0) - xyz_new[i]
    output = minimize(_minimize_func, trial_vec, args=(vec2, xyz_new, i))
    rotmat2 = rotation_matrix(output.x, vec2)
    return (rotmat1.T@rotmat2.T).T


def _minimize_func(vec1: np.ndarray, vec2: np.ndarray, xyz: np.ndarray, anchor: int) -> float:
    """The function whose output is to-be minimized by :func:`scipy.optimize.minimize`."""
    # Rotate and translate
    rotmat = rotation_matrix(vec1, vec2)
    xyz_new = xyz@rotmat
    xyz_new -= xyz_new[anchor]

    # Apply the cost function: e^(|yz|)
    yz = xyz_new[:, 1:]
    distance = np.linalg.norm(yz, axis=1)
    return np.exp(distance).sum()
