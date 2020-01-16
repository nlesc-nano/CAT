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
    out5 = distribute_idx(MOL, IDX, p=0.5, mode='random')
    out6 = distribute_idx(MOL, IDX, p=0.5, mode='uniform', follow_edge=True)
    out7 = distribute_idx(MOL, IDX, p=0.5, mode='cluster', follow_edge=True)
    out8 = distribute_idx(MOL, IDX, p=0.5, mode='uniform', cluster_size=2)
    out9 = distribute_idx(MOL, IDX, p=0.5, mode='uniform', cluster_size=range(1, 4))

    np.testing.assert_array_equal(out1, [127, 130, 138, 147, 144, 148, 131, 146, 137, 136, 124, 134, 123])  # noqa
    np.testing.assert_array_equal(out2, [132, 147, 126, 137, 123, 130, 133, 141, 124, 129, 148, 136, 142])  # noqa
    np.testing.assert_array_equal(out3, [148, 138, 127, 141, 126, 144, 123, 131, 146, 136, 134, 133, 132])  # noqa
    np.testing.assert_array_equal(out4, [148, 129, 136, 142, 123, 132, 147, 126, 137, 130, 133, 141, 124])  # noqa
    np.testing.assert_array_equal(out6, [127, 130, 128, 129, 144, 126, 139, 140, 142, 132, 124, 134, 135])  # noqa
    np.testing.assert_array_equal(out7, [135, 131, 139, 143, 127, 134, 125, 144, 145, 136, 142, 148, 129])  # noqa
    np.testing.assert_array_equal(out8, [127, 123, 146, 128, 124, 141, 126, 137, 148, 129, 145, 144, 139])  # noqa
    np.testing.assert_array_equal(out9, [127, 130, 141, 126, 137, 147, 140, 148, 129, 144, 145, 125, 139])  # noqa
    assertion.len_eq(out5, round(0.5 * len(IDX)))
    assertion.len_eq(np.intersect1d(out5, IDX), len(out5))

    assertion.assert_(distribute_idx, MOL, IDX, p=2, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=-1, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, mode='bob', exception=ValueError)
    assertion.assert_(distribute_idx, MOL, ['bob'], p=0.5, exception=TypeError)
    assertion.assert_(distribute_idx, MOL, ['bob'], p=0.5, exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, cluster_size=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, cluster_size='bob', exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, cluster_size=['bob'], exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, start='bob', exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, p=0.5, start=9999, exception=ValueError)
