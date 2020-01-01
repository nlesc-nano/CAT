"""
CAT.attachment.edge_distance
============================

Functions for creating distance matrices from

Index
-----
.. currentmodule:: CAT.attachment.edge_distance
.. autosummary::
    edge_dist
    array_combinations
    to_convex

API
---
.. autofunction:: edge_dist
.. autofunction:: array_combinations
.. autofunction:: to_convex

"""

import reprlib
from itertools import combinations
from math import factorial

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import ConvexHull

__all__ = ['edge_dist']


def array_combinations(array: np.ndarray, r: int = 2) -> np.ndarray:
    """Construct an array with all :func:`combinations<itertools.combinations>` of **ar** along axis 1.

    Parameters
    ----------
    array : array-like
        A 2D array-like object.

    r : :class:`int`
        The length of each combination.

    Returns
    -------
    ``(n, len(ar), r)`` :class:`numpy.ndarray`
        A 3D array with all **ar** combinations (of length ``e``) along axis 1.
        ``n`` represents the number of combinations: :math:`n! / r! / (n-r)!`.

    """  # noqa
    ar = np.asarray(array)
    if ar.ndim != 2:
        raise ValueError(f"'array' excpected a 2D array; observed dimensionality: {ar.ndim}")
    n = ar.shape[1]
    combinations_len = int(factorial(n) / factorial(r) / factorial(n - r))

    shape = combinations_len, len(ar), r
    ret = np.empty(shape, dtype=ar.dtype)
    for i, jk in enumerate(combinations(range(ar.shape[1]), r=r)):
        ret[i] = ar[:, jk]
    return ret


def to_convex(xyz: np.ndarray, n: float = 1.0) -> np.ndarray:
    r"""Round all edges in **xyz** by a factor **n**: (:math:`0 < n \le 1`)."""
    if not (0.0 < n <= 1.0):
        raise ValueError(f"Expected '0 < n <= 1'; observed: {reprlib.repr(n)}")
    xyz = np.asarray(xyz, dtype=float)
    xyz = xyz - xyz.mean(axis=0)

    r = np.linalg.norm(xyz, axis=-1)
    dr = n * (r.mean() - r)

    elongate = (dr + r) / r
    xyz *= elongate[..., None]
    return xyz


def edge_dist(xyz: np.ndarray, n: float = 1.0) -> np.ndarray:
    r"""Calculate all shortest paths in the polyhedron **xyz** by traversing its edges.

    After converting **xyz** into a polyhedron with triangular faces,
    the shortest paths between all possible point pairs is calculated.
    Paths are constrained to the edges of the polyhedron,
    effectively constraining all movement to a 2D surface
    as opposed to the 3D volume of an ordinary distance matrix.

    Given the matrix of Cartesian coordinates :math:`X \in \mathbb{R}^{n, 3}`,
    the matching distance matrix :math:`D \in \mathbb{R}^{n, n}` and
    the matching edge-distance matrix :math:`D^{\text{edge}} \in \mathbb{R}^{n, n}`,
    then element :math:`D_{i, j}` is defined as following:

    .. math::

        D_{i, j} = ||X_{i,:} - X_{j,:}||_{2}

    while :math:`D_{i, j}^{\text{edge}}` is defined as:

    .. math::

        D_{i, j}^{\text{edge}} = \min_{\boldsymbol{k}} \Bigl{(}
            ||X_{i,:} - X_{\boldsymbol{k}_{0},:}||_{2} + ... +
            ||X_{\boldsymbol{k}_{m},:} - X_{j,:}||_{2} \Bigr{)}

        \quad \text{with} \quad
        \boldsymbol{k}_{0} \ne i
        \quad \text{and} \quad
        \boldsymbol{k}_{m} \ne j

    The vector :math:`\boldsymbol{k} \in \mathbb{Z}^{m}` is a path,
    represented by the indices of neighbouring vertices in :math:`X`,
    whose length (*i.e.* Euclidean distance) is to-be minimized.

    Parameters
    ----------
    xyz : :class:`numpy.ndarray`
        A 2D array-like object of Cartesian coordinates representing a polyhedron.
        The supplied polyhedron should be convex in shape.

    n : :class:`float`
        Smoothing factor for constructing a convex hull.
        Should obey :math:`0 <= n <= 1`.

    Returns
    -------
    :class:`numpy.ndarray`
        A 2D array containing all possible (Euclidean) distance-pairs in **xyz**.
        Distances are calculated by traversing the (minimum-length) edges of **xyz**,
        rather than moving directly through space.

    See Also
    --------
    :class:`ConvexHull<scipy.spatial.ConvexHull>`
        Convex hulls in N dimensions.

    :func:`dijkstra<scipy.sparse.csgraph.dijkstra>`
        Dijkstra algorithm using Fibonacci Heaps.

    :func:`cdist<scipy.spatial.distance.cdist>`
        Compute distance between each pair of the two collections of inputs.

    """
    xyz = np.array(xyz, dtype=float, ndmin=2, copy=False)
    xyz_convex = to_convex(xyz, n=n) if not np.allclose(n, 0.0) else xyz
    xyz_len = len(xyz)

    # Create an array with all index-pairs forming the polyhedron edges
    hull = ConvexHull(xyz_convex)
    idx = array_combinations(hull.simplices, r=2)
    idx.shape = -1, 2
    i, j = idx.T

    # Create a distancea matrix containing all edge distances
    dist = np.zeros((xyz_len, xyz_len), dtype=float)
    dist[i, j] = dist[j, i] = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    dist_sparse = csr_matrix(dist)

    # Traverse the edges
    ret = np.zeros_like(dist)
    for k in range(xyz_len):
        ret[k] = dijkstra(dist_sparse, indices=k, return_predecessors=False, min_only=True)
    return ret


def _plot_polyhedron(xyz: np.ndarray, show: bool = True) -> 'matplotlib.pyplot.Figure':
    """Plot a polyhedron, represented by an array of Cartesian coordinates, with matplotlib."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    xyz = np.asarray(xyz, dtype=float)
    hull = ConvexHull(xyz)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_trisurf(*xyz.T, triangles=hull.simplices, cmap=plt.cm.Spectral)

    if show:
        plt.show(block=True)
    return fig
