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
    out6 = distribute_idx(MOL, IDX, p=0.5, mode='uniform', follow_edge=True)
    out7 = distribute_idx(MOL, IDX, p=0.5, mode='cluster', follow_edge=True)

    np.testing.assert_array_equal(out1, [127, 123, 124, 126, 128, 129, 131, 144, 140, 142, 137, 139, 145])  # noqa
    np.testing.assert_array_equal(out2, [136, 142, 129, 148, 135, 127, 143, 125, 131, 147, 123, 134, 126])  # noqa
    np.testing.assert_array_equal(out3, [148, 123, 124, 125, 126, 128, 131, 127, 140, 145, 136, 137, 139])  # noqa
    np.testing.assert_array_equal(out4, [148, 129, 142, 136, 133, 123, 130, 147, 124, 125, 132, 127, 144])  # noqa
    np.testing.assert_array_equal(out6, [128, 123, 124, 125, 126, 129, 131, 143, 140, 145, 142, 132, 141])  # noqa
    np.testing.assert_array_equal(out7, [135, 131, 143, 139, 127, 136, 126, 129, 142, 134, 147, 125, 148])  # noqa
    assertion.len_eq(out5, round(0.5 * len(IDX)))
    assertion.len_eq(np.intersect1d(out5, IDX), len(out5))

    assertion.assert_(distribute_idx, MOL, IDX, p=2, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=-1, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, mode='bob', exception=ValueError)
