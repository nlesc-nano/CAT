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


def optimize_rotmat(mol: np.ndarray, anchor: int = 0, as_vec: bool = False) -> np.ndarray:
    r"""Find the rotation matrix for **xyz** that minimizes its deviation from the Cartesian X-axis.

    A set of vectors, :math:`v`, is constructed for all atoms in **xyz**,
    denoting their deviation from the Cartesian X-axis.
    Subsequently, the rotation matrix that minimizes :math:`\sum_{i}^{n} {e^{v_{i}}}` is returned.

    Parameters
    ----------
    : : :math:`n*3` :class:`numpy.ndarray` [:class:`float`]
        An array-like object of Cartesian coordinates.
        *e.g.* :class:`list`, :class:numpy.ndarray: or |plams.Molecule|.

    anchor : |plams.Atom| or :class:`int`
        The index (0-based) of the anchor atom in **xyz**.
        Used for defining the origin in **xyz** (*i.e.* :code:`xyz[i] == (0, 0, 0)`).
        Alternativelly, a PLAMS atom can be passed.

    as_vec : :class:`bool`
        If ``True``, return the initial vector used for constructing the rotation matrix.

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
        i = anchor.mol.index(anchor)
    else:
        raise TypeError("The passed anchor is neither an 'int' nor 'Atom'; "
                        f"observed type: {repr(type(anchor))}")
    xyz = np.array(mol, dtype=float, ndmin=2, copy=False)

    # Construct the initial vectors
    vec1_trial = xyz.mean(axis=0) - xyz[i]
    vec2 = np.array([1, 0, 0], dtype=float)

    # Optimize the trial vector; return the matching rotation matrix
    output = minimize(_minimize_func, vec1_trial, args=(vec2, xyz, i))
    vec1 = output.x
    if as_vec:
        return vec1
    return rotation_matrix(vec1, vec2)


def _minimize_func(vec1: np.ndarray, vec2: np.ndarray, xyz: np.ndarray, anchor: int) -> float:
    """The function whose output is to-be minimized by :func:`scipy.optimize.minimize`."""
    # Rotate and translate
    rotmat = rotation_matrix(vec1, vec2)
    xyz = xyz@rotmat
    xyz -= xyz[anchor]

    # Apply the cost function: e^(|yz|)
    yz = xyz[:, 1:]
    distance = np.linalg.norm(yz, axis=1)
    return np.exp(distance).sum()
