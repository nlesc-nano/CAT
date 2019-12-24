"""Tests for :mod:`CAT.attachment.distribution`."""

from pathlib import Path

import numpy as np

from scm.plams import Molecule
from assertionlib import assertion

from CAT.attachment.distribution import distribute_idx

PATH = Path('tests') / 'test_files'
PATH = Path(r'D:\hardd\Documents\GitHub\CAT\tests\test_files')
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

    np.testing.assert_array_equal(out1, [132, 141, 134, 137, 143, 138, 148, 133, 139, 147, 130, 125, 131])  # noqa
    np.testing.assert_array_equal(out2, [130, 142, 140, 143, 124, 134, 139, 132, 141, 146, 125, 148, 133])  # noqa
    np.testing.assert_array_equal(out3, [141, 134, 137, 143, 138, 148, 133, 139, 147, 130, 125, 131, 142])  # noqa
    np.testing.assert_array_equal(out4, [131, 137, 144, 130, 142, 140, 143, 124, 134, 139, 132, 141, 146])  # noqa
    assertion.len_eq(out5, round(0.5 * len(IDX)))
    assertion.len_eq(np.intersect1d(out5, IDX), len(out5))

    assertion.assert_(distribute_idx, MOL, IDX, p=2, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=-1, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, mode='bob', exception=ValueError)
