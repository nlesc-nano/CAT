"""
CAT.attachment.perp_surface
===========================

A module for constructing vectors perpendicular to a molecules surface.

Index
-----
.. currentmodule:: CAT.attachment.perp_surface
.. autosummary::
    get_surface_vec

API
---
.. autofunction:: get_surface_vec

"""

from typing import Union

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

from scm.plams import Molecule

__all__ = ['get_surface_vec']


def get_surface_vec(mol: Union[Molecule, np.ndarray],
                    anchor: Union[Molecule, np.ndarray]) -> np.ndarray:
    """Construct a set of vectors perpendicular to the surface of **mol** and assign one to each atom in **anchor**.

    Utilizes a convex hull algorithm for identifying and partitioning the surface.

    Parameters
    ----------
    mol : array-like [:class:`float`], shape :math:`(n, 3)`
        A 2D array-like object representing the Cartesian coordinates of a polyhedron.

    anchor : array-like [:class:`float`], shape :math:`(m, 3)`
        A 2D array-like object representing the Cartesian coordinates of a polyhedron.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`float`], shape :math:`(m, 3)`
        An array with vectors perpendicular to the surface of **mol**,
        one for each atom in **anchor**.

    See Also
    --------
    :class:`ConvexHull<scipy.spatial.ConvexHull>`
        Convex hulls in N dimensions.

    """  # noqa
    xyz = np.array(mol, dtype=float, ndmin=2, copy=False)
    anchor = np.array(anchor, dtype=float, ndmin=2, copy=False)

    hull = ConvexHull(xyz)
    simplice = np.swapaxes(xyz[hull.simplices], 0, 1)
    simplice_center = simplice.mean(axis=0)

    vec = _get_perp_vecs(*simplice)
    _flip_vec(simplice_center, vec)
    return vec[_find_nearest_center(anchor, simplice_center)]


def _find_nearest_center(anchor: np.ndarray, center: np.ndarray) -> np.ndarray:
    """Find the points in **center** closest to those in **anchor**."""
    tree = cKDTree(center)
    _, idx = tree.query(anchor, k=1)
    return idx


def _get_perp_vecs(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Construct a unit vector perpendicular to a set of triangular polygons."""
    v1 = p2 - p1
    v2 = p3 - p1
    vec = np.cross(v1, v2)
    vec /= np.linalg.norm(vec, axis=-1)[..., None]
    return vec


def _flip_vec(simplice_centra: np.ndarray, perp_vec: np.ndarray) -> None:
    """Ensure that **perp_vec** is pointing away from the surface defined by **simplice_centra**."""
    center = simplice_centra.mean(axis=0)
    dot = np.einsum('ij,ij->i', simplice_centra - center, perp_vec)
    perp_vec[dot < 0] *= -1
