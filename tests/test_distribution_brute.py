"""Tests for :mod:`CAT.attachment.distribution_brute`."""

from pathlib import Path
from itertools import cycle
from typing import Iterable, Generator

import numpy as np

from scm.plams import Molecule
from assertionlib import assertion

from CAT.attachment import brute_uniform_idx

PATH = Path('tests') / 'test_files'
MOL = Molecule(PATH / 'core' / 'Cd68Se55.xyz')
IDX = np.fromiter([i for i, at in enumerate(MOL) if at.symbol == 'Cd'], dtype=int)
IDX.shape = len(IDX) // 4, 4
IDX.setflags(write=False)


def cycle_slice(iterable: Iterable[int], start: int = 0) -> Generator[slice, None, None]:
    """``cycle_slice(range(3))`` --> ``slice(0, 0)``, ``slice(0, 1)``, ``slice(1, 3)``, ``slice(3, 3)``, ..."""  # noqa
    for i in cycle(iterable):
        yield slice(start, start+i)
        start += i


def test_distribute_idx() -> None:
    """Test :func:`CAT.attachment.distribution_brute.brute_uniform_idx`."""
    ar_list = [
        brute_uniform_idx(MOL, IDX, n=1, operation='min'),
        brute_uniform_idx(MOL, IDX, n=2, operation='min'),
        brute_uniform_idx(MOL, IDX, n=3, operation='min'),
        brute_uniform_idx(MOL, IDX, n=4, operation='min'),

        brute_uniform_idx(MOL, IDX, n=1, operation='max'),
        brute_uniform_idx(MOL, IDX, n=2, operation='max'),
        brute_uniform_idx(MOL, IDX, n=3, operation='max'),
        brute_uniform_idx(MOL, IDX, n=4, operation='max')
    ]

    ref_tot = np.load(PATH / 'brute_uniform_idx.npy')
    for i, ar in zip(cycle_slice(range(1, 5)), ar_list):
        ref = ref_tot[:, i]
        np.testing.assert_array_equal(ar, ref)

    assertion.assert_(brute_uniform_idx, MOL, IDX, operation='min', exception=ValueError)
    assertion.assert_(brute_uniform_idx, MOL, IDX, n=0, exception=ValueError)
    assertion.assert_(brute_uniform_idx, MOL, IDX, n=99, exception=ValueError)
