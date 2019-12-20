"""
CAT.attachment.distribution
===========================

Functions for creating distributions of atomic indices (*i.e.* core anchor atoms).

"""

import random
import reprlib
from typing import Generator, Optional, Iterable, TypeVar, List, Sequence, FrozenSet, Any
from itertools import islice

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import Molecule

__all__ = ['get_distribution']


#: A set of allowed values for the **mode** parameter in :func:`get_distribution`.
MODE_SET: FrozenSet[str] = frozenset({'uniform', 'random'})


def get_distribution(core: Molecule, idx: Sequence[int], mode: str = 'uniform',
                     p: float = 0.5, **kwargs: Any) -> Sequence[int]:
    # Validate the input
    if mode not in MODE_SET:
        raise ValueError(f"Invalid value for 'mode' ({reprlib.repr(mode)}); "
                         f"accepted values: {tuple(MODE_SET)}")
    elif not (0.0 < p <= 1.0):
        raise ValueError("'p' should be larger than 0.0 and smaller than or equal to 1.0; "
                         f"observed value: {repr(p)}")
    elif p == 1.0:
        return idx

    # Create an iterator of atomic indices
    if mode == 'uniform':
        xyz = np.array(core)[idx]
        dist = cdist(xyz, xyz)
        generator = (idx[i] for i in uniform(dist, kwargs.get('start', None)))
    elif mode == 'random':
        generator = random_element(idx)

    # Return a list of `p * len(idx)` atomic indices
    stop = max(1, round(p * len(idx)))
    return [i for i in islice(generator, 0, stop)]


def uniform(dist: np.ndarray, start: Optional[int] = None) -> Generator[int, None, None]:
    r"""Return the column-indices of **dist**, :math:`d`, which yield a uniform distribution.

    Given a square distance matrix :math:`D` and
    a vector :math:`\boldsymbol{d}` (representing a set of indices in :math:`D`),
    then the :math:`i`'th element :math:`\boldsymbol{d}_{i}` is defined as following:

    .. math::

        \DeclareMathOperator*{\argmin}{\arg\!\min}

        \boldsymbol{d}_{i} = \argmin_{k} {\boldsymbol{d}_{k}}^{*}
        \quad \text{with} \quad
        \boldsymbol{d}^{*} = \sqrt{{D_{i,:}}^{2} + \sum_{0<j<i} {D_{j, \boldsymbol{d}_{j}}}^{2}}

    The row in :math:`D` corresponding to :math:`\boldsymbol{d}_{i=0}`
    can be specified by **start**.

    Parameters
    ----------
    dist : :math:`(m, m)` :class:`numpy.ndarray`
        A square 2D NumPy array representing the distance matrix :math:`D`.

    start : :class:`int`, optional
        The index of the starting row.
        Will be randomized if ``None``.

    Yields
    ------
    :class:`int`
        Column-indices specified in :math:`\boldsymbol{d}`.

    """
    shift = start if start is not None else random.randint(0, len(dist)-1)
    dist_sqr = np.asarray(dist)**2  # Squared distance matrix
    dist_sqr[:] = np.roll(dist_sqr, shift, axis=0)  # Shift the first row to **shift**
    norm_sqr: float = 0.0  # Squared norm

    for i, ar in enumerate(dist_sqr):
        ar_sqrt = np.sqrt(ar + norm_sqr)

        j = ar_sqrt.argmax()
        norm_sqr += ar[j]

        dist_sqr[:, j] = 0.0
        yield j


T = TypeVar('T')


def random_element(iterable: Iterable[T]) -> Generator[T, None, None]:
    r"""Yield random elements from an iterable.

    No single element is returned more than once.

    Parameters
    ----------
    iterable : :class:`Iterable<collections.abc.Iterable>`
        An iterable.

    Yields
    ------
    :class:`object`
        Random elements from **iterable**.

    """
    try:
        k = len(iterable)
    except TypeError:  # It's an iterator
        iterable = list(iterable)
        k = len(iterable)

    rand_element = random.sample(iterable, k=k)
    for i in rand_element:
        yield i
