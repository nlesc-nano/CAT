"""Tests for :mod:`CAT.attachment.distribution`."""

from pathlib import Path

import numpy as np

from scm.plams import Molecule
from assertionlib import assertion

from CAT.attachment.distribution import distribute_idx

PATH = Path('tests') / 'test_files'
MOL = Molecule(str(PATH / 'core' / 'Cd68Se55.xyz'))
IDX = np.array([i for i, at in enumerate(MOL) if at.symbol == 'Cl'])
IDX.setflags(write=False)


def test_distribute_idx() -> None:
    """Test :func:`CAT.attachment.distribution.random_idx`."""
    out1 = distribute_idx(MOL, IDX, p=0.5, mode='uniform')
    out2 = distribute_idx(MOL, IDX, p=0.5, mode='cluster')
    out3 = distribute_idx(MOL, IDX, p=0.5, mode='uniform', start=-1)
    out4 = distribute_idx(MOL, IDX, p=0.5, mode='cluster', start=-1)
    out5 = distribute_idx(MOL, IDX, mode='random', p=0.5)

    np.testing.assert_array_equal(out1, [129, 139, 141, 144, 133, 137, 123, 148, 132, 131, 127, 145, 142])  # noqa
    np.testing.assert_array_equal(out2, [134, 146, 144, 147, 128, 138, 143, 136, 145, 124, 129, 126, 137])  # noqa
    np.testing.assert_array_equal(out3, [142, 135, 138, 144, 139, 123, 134, 140, 148, 131, 126, 132, 143])  # noqa
    np.testing.assert_array_equal(out4, [132, 138, 145, 131, 143, 141, 144, 125, 135, 140, 133, 142, 147])  # noqa
    assertion.len_eq(out5, round(0.5 * len(IDX)))
    assertion.len_eq(np.intersect1d(out5, IDX), len(out5))

    assertion.assert_(distribute_idx, MOL, IDX, p=2, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=-1, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, mode='bob', exception=ValueError)
