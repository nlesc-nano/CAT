"""
CAT.attachment.distribution
===========================

Functions for creating distributions of atomic indices (*i.e.* core anchor atoms).

Index
-----
.. currentmodule:: CAT.attachment.distribution
.. autosummary::
    distribute_idx
    uniform_idx
    random_idx

API
---
.. autofunction:: distribute_idx
.. autofunction:: uniform_idx
.. autofunction:: random_idx

"""

import reprlib
from itertools import islice
from typing import Generator, Optional, Iterable, FrozenSet, Any, Union

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import Molecule

from .edge_distance import edge_dist

__all__ = ['distribute_idx']

#: A set of allowed values for the **mode** parameter in :func:`get_distribution`.
MODE_SET: FrozenSet[str] = frozenset({'uniform', 'random', 'cluster'})


def distribute_idx(core: Union[Molecule, np.ndarray], idx: Union[int, Iterable[int]], p: float,
                   mode: str = 'uniform', **kwargs: Any) -> np.ndarray:
    r"""Create a new distribution of atomic indices from **idx** of length :code:`p * len(idx)`.

    Parameters
    ----------
    core : :math:`(m, 3)` array-like [:class:`float`]
        A 2D array-like object (such as a :class:`Molecule` instance) consisting
        of Cartesian coordinates.

    idx : :class:`int` or :math:`(i,)` :class:`Iterable<collections.abc.Iterable>` [:class:`int`]
        An integer or iterable of unique integers representing the 0-based indices of
        all anchor atoms in **core**.

    p : :class:`float`
        A float obeying the following condition: :math:`0.0 < p <= 1.0`.
        Represents the fraction of **idx** that will be returned.

    mode : :class:`str`
        How the subset of to-be returned indices will be generated.
        Accepts one of the following values:

        * ``"random"``: A random distribution.
        * ``"uniform"``: A uniform distribution; the distance between each successive atom and
          all previous points is maximized.
        * ``"cluster"``: A clustered distribution; the distance between each successive atom and
          all previous points is minmized.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for the **mode**-specific functions.

    Returns
    -------
    :math:`(p*i,)` :class:`numpy.ndarray` [:class:`int`]
        A 1D array of atomic indices.
        If **idx** has :math:`i` elements,
        then the length of the returned list is equal to :math:`\max(1, p*i)`.

    See Also
    --------
    :func:`uniform_idx`
        Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    :func:`cluster_idx`
        Return the column-indices of **dist** which yield a clustered distribution.

    """
    # Convert **idx** into an array
    try:
        idx_ar = np.array(idx, dtype=int, ndmin=1, copy=False)
    except TypeError:  # A Collection or Iterator
        idx_ar = np.fromiter(idx, dtype=int)

    # Validate the input
    if mode not in MODE_SET:
        raise ValueError(f"Invalid value for 'mode' ({reprlib.repr(mode)}); "
                         f"accepted values: {reprlib.repr(tuple(MODE_SET))}")
    elif not (0.0 < p <= 1.0):
        raise ValueError("'p' should be larger than 0.0 and smaller than or equal to 1.0; "
                         f"observed value: {reprlib.repr(p)}")
    elif p == 1.0:  # Ensure that **idx** is always returned as copy
        return idx_ar.copy() if idx_ar is idx else idx_ar

    # Create an array of indices
    stop = max(1, round(p * len(idx_ar)))
    if mode in ('uniform', 'cluster'):
        xyz = np.array(core, dtype=float, ndmin=2, copy=False)[idx_ar]
        dist = edge_dist(xyz) if kwargs.get('follow_edge', False) else cdist(xyz, xyz)
        if mode == 'uniform':
            generator1 = uniform_idx(dist, 'max', p=p, start=kwargs.get('start', None))
            generator2 = islice(generator1, stop)
            ret = idx_ar[np.fromiter(generator2, count=stop, dtype=int)]
        else:
            ret = idx_ar[cluster_idx(dist, start=kwargs.get('start', None))]

    elif mode == 'random':
        ret = np.random.permutation(idx_ar)

    # Return a list of `p * len(idx)` atomic indices
    return ret[:stop]


def uniform_idx(dist: np.ndarray, operation: str = 'max', p: Optional[float] = 0.5,
                start: Optional[int] = None) -> Generator[int, None, None]:
    r"""Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    Given the (symmetric) distance matrix :math:`D \in \mathbb{R}^{n,n}` and
    the vector :math:`\boldsymbol{d} \in \mathbb{Z}^{m}`
    (representing a subset of indices in :math:`D`),
    then the :math:`i`'th element :math:`\boldsymbol{d}_{i}` is
    defined as following:

    .. math::

        \DeclareMathOperator*{\argmax}{\arg\!\max}
        \boldsymbol{d}_{i} = \argmax_{k}
        \sqrt{ \sum_{0 \le j < i} {D_{k, \boldsymbol{d}_{j}}}^2 }
        \quad \text{with} \quad
        k \notin \boldsymbol{d}[0, ..., i-1]

    The row in :math:`D` corresponding to :math:`\boldsymbol{d}_{i=0}`
    can be specified by **start**.

    THe distance matrix can be truncated with the **p** parameter.

    The :math:`\text{argmax}` operation can be exchanged for :math:`\text{argmin}` by settings
    **operation** to ``"min"``, thus yielding a clustered- rather than uniform-distribution.

    Parameters
    ----------
    dist : :math:`(m, m)` :class:`numpy.ndarray` [:class:`float`]
        A symmetric 2D NumPy array (:code:`(dist == dist.T).all()`) representing the
        distance matrix :math:`D`.

    operation : :class:`str`
        Whether to minimize or maximize the distance between points.
        Accepted values are ``"min"`` and ``"max"``.

    p : :class:`float`, optional
        A float obeying the following condition: :math:`0.0 < p <= 1.0`.
        Represents the fraction of **dist** which is of interest to the user.
        If not ``None``, used for truncating the distance matrix :math:`D`:

        .. math::

            r_{truncate} = \max(2, r_{nn} * \log_{2} p)
            \quad \text{with} \quad
            r_{nn} = \frac{1}{N} \sum_{i=0}^{N} D_{i,:}

    start : :class:`int`, optional
        The index of the starting row in **dist**.
        If ``None``, start in whichever row contains the global minimum/maximum (see **operation**).

    Yields
    ------
    :class:`int`
        Column-indices specified in :math:`\boldsymbol{d}`.

    """  # noqa
    if operation not in ('min', 'max'):
        raise ValueError(f"Invalid value for 'mode' ({reprlib.repr(operation)}); "
                         f"accepted values: ('min', 'max')")

    # Truncate and square the distance matrix
    dist_sqr = np.array(dist, dtype=float, copy=True)
    if p is not None:
        np.fill_diagonal(dist_sqr, np.inf)
        n = max(2, -np.log2(p))
        trunc = n * dist_sqr.min(axis=0).mean()
        dist_sqr[dist_sqr > trunc] = trunc
    dist_sqr **= 2

    # Use either argmin or argmax
    if operation == 'min':
        arg_func = np.nanargmin
    else:
        arg_func = np.nanargmax
    start = arg_func(np.linalg.norm(dist, axis=1)) if start is None else start
    np.fill_diagonal(dist_sqr, np.nan)

    # Yield indices
    dist_1d_sqr = dist_sqr[start].copy()
    dist_1d_sqr[start] = np.nan
    yield start

    for _ in range(len(dist_1d_sqr)-1):
        dist_1d = dist_1d_sqr**0.5
        i = arg_func(dist_1d)
        dist_1d_sqr[i] = np.nan
        dist_1d_sqr += dist_sqr[i]
        yield i


def cluster_idx(dist: np.ndarray, start: Optional[int] = None) -> np.ndarray:
    r"""Return the column-indices of **dist** which yield a clustered distribution.

    Given the (symmetric) distance matrix :math:`D \in \mathbb{R}^{n, n}` and the row index
    :math:`\DeclareMathOperator*{\argmin}{\arg\!\min} i = \argmin_{i} ||D_{i, :}||_{2}`,
    return the column-indices of :math:`D_{i, :}` sorted in order of ascending distance.

    Parameters
    ----------
    dist : :math:`(m, m)` array-like [:class:`float`]
        A symmetric 2D NumPy array (:code:`(dist == dist.T).all()`) representing the
        distance matrix :math:`D`.

    start : :class:`int`, optional
        The index of the starting row in **dist**.
        If ``None``, start in row:
        :math:`\DeclareMathOperator*{\argmin}{\arg\!\min} i = \argmin_{i} ||D_{i, :}||_{2}`.

    Returns
    -------
    :math:`(m,)` :class:`numpy.ndarray` [:class:`int`]
        A 1D array of indices.

    """
    dist = np.asarray(dist, dtype=float)
    start = np.linalg.norm(dist, axis=1).argmin() if start is None else start

    r = dist[start]
    r_arg = r.argsort()
    idx = np.arange(len(r))
    return idx[r_arg]


def truncate_dist(dist: np.ndarray, p: float) -> None:
    np.fill_diagonal(dist, np.inf)
    nn_dist = dist.min(axis=0).mean()
    base = get_nn_count(dist, nn_dist)

    print(base)
    n = -np.log(p) / np.log(base)
    n += 2
    print(n)
    n *= nn_dist
    dist[dist > n] = n


def get_nn_count(dist: np.ndarray, r: float, r_min: float = 0.5, r_max: float = 1.5) -> float:
    """Return the number of elements in **dist** whose value are within the range :math:`[r*r_{min}, r*r_{max}]`."""  # noqa
    valid = (dist > r_min * r) & (dist < r_max * r)
    valid.shape = valid.size
    n = np.bincount(valid)[1]
    n /= len(dist)
    return n
