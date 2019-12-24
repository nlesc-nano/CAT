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
from typing import Generator, Optional, Iterable, FrozenSet, Any, Union
from itertools import islice

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import Molecule

__all__ = ['distribute_idx']

#: A set of allowed values for the **mode** parameter in :func:`get_distribution`.
MODE_SET: FrozenSet[str] = frozenset({'uniform', 'random', 'cluster'})


def distribute_idx(core: Union[Molecule, np.ndarray], idx: Union[int, Iterable[int]], p: float,
                   mode: str = 'uniform', **kwargs: Any) -> np.ndarray:
    r"""Create a new distribution of atomic indices from **idx** of length :code:`p * len(idx)`.

    Parameters
    ----------
    core : array-like [:class:`int`]
        A 2D array-like object (such as a :class:`Molecule` instance) consisting
        of Cartesian coordinates.

    idx : array-like [:class:`int`]
        An integer or iterable of integers representing the 0-based indices of
        all anchor atoms in **core**.

    p : :class:`float`
        A float obeying the following condition: :math:`0.0 < p <= 1.0`.
        Represents the fraction of **idx** that will be returned.

    mode : :class:`str`
        How the subset of to-be returned indices will be generated.
        Accepts one of the following values:

        * ``"random"``: A random distribution.
        * ``"uniform"``: A uniform distribution; the distance between each successive point and
          all previous points is maximized.
        * ``"cluster"``: A clustered distribution; the distance between each successive point and
          all previous points is minmized.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for the **mode**-specific functions.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`int`]
        An array of atomic indices.
        If **idx** has :math:`i` elements,
        then the length of the returned list is equal to :math:`\max(1, p*i)`.

    See Also
    --------
    :func:`uniform_idx`
        Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    :func:`random_idx`
        Yield random elements from an **iterable**.

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
    if mode in ('uniform', 'cluster'):
        xyz = np.array(core, dtype=float, ndmin=2, copy=False)[idx_ar]
        dist = cdist(xyz, xyz)
        operation = 'max' if mode == 'cluster' else 'min'
        iterable = (idx_ar[i] for i in uniform_idx(dist, operation, kwargs.get('start', None)))
    elif mode == 'random':
        iterable = np.random.permutation(idx_ar)

    # Return a list of `p * len(idx)` atomic indices
    stop = max(1, round(p * len(idx_ar)))
    return np.fromiter(islice(iterable, stop), count=stop, dtype=int)


def uniform_idx(dist: np.ndarray, operation: str = 'max',
                start: Optional[int] = None) -> Generator[int, None, None]:
    r"""Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    Given the symmetric distance matrix :math:`D` and
    a vector :math:`\boldsymbol{d}` (representing a set of indices in :math:`D`),
    then the :math:`i`'th element :math:`\boldsymbol{d}_{i}` is defined as following:

    .. math::

        \DeclareMathOperator*{\argmax}{\arg\!\max}

        \boldsymbol{d}_{i} = \argmax_{k} {\boldsymbol{d}_{k}}^{*}
        \quad \text{with} \quad
        \boldsymbol{d}^{*} = \sqrt{{D_{i,:}}^{2} + \sum_{0<j<i} {D_{j, \boldsymbol{d}_{j}}}^{2}}

    The row in :math:`D` corresponding to :math:`\boldsymbol{d}_{i=0}`
    can be specified by **start**.

    The :math:`\text{argmax}` operation can be exchanged for :math:`\text{argmin}` by settings
    **operation** to ``"min"``, thus yielding a clustered- rather than uniform-distribution.

    Parameters
    ----------
    dist : :math:`(m, m)` :class:`numpy.ndarray`
        A symmetric 2D NumPy array (:code:`(dist == dist.T).all()`) representing the
        distance matrix :math:`D`.

    operation : :class:`str`
        Whether to minimize or maximize the distance between points.
        Accepted values are ``"min"`` and ``"max"``.

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
    dist_sqr = np.asarray(dist)**2  # Squared distance matrix
    norm_sqr = 0.0  # Squared norm

    # Use either argmin or argmax
    if operation == 'min':
        fill_value = np.inf
        np.fill_diagonal(dist_sqr, np.inf)
        arg_func = np.ndarray.argmin
    else:
        fill_value = 0.0
        arg_func = np.ndarray.argmax
    shift = start if start is not None else np.unravel_index(arg_func(dist), dist.shape)[0]
    dist_sqr[:] = np.roll(dist_sqr, shift, axis=0)  # Shift the first row to **shift**

    # Yield indices
    for i, ar in enumerate(dist_sqr):
        ar_sqrt = np.sqrt(ar + norm_sqr)

        j = arg_func(ar_sqrt)
        norm_sqr += ar[j]

        dist_sqr[:, j] = fill_value
        yield j
