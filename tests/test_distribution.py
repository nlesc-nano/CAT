"""Tests for :mod:`CAT.attachment.distribution`."""

from pathlib import Path

import yaml
import pytest
import numpy as np
from scm.plams import Molecule, Settings, readpdb
from assertionlib import assertion
from nanoutils import delete_finally

from CAT.test_utils import assert_mol_allclose
from CAT.base import prep
from CAT.workflows import MOL as MOL_KEY
from CAT.attachment.distribution import distribute_idx

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'

MOL = Molecule(PATH / 'core' / 'Cd68Se55.xyz')

IDX = np.array([i for i, at in enumerate(MOL) if at.symbol == 'Cl'])
IDX.setflags(write=False)


def test_distribute_idx() -> None:
    """Test :func:`CAT.attachment.distribution.random_idx`."""
    out1 = distribute_idx(MOL, IDX, f=0.5, mode='uniform')
    out2 = distribute_idx(MOL, IDX, f=0.5, mode='cluster')
    out3 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', start=-1)
    out4 = distribute_idx(MOL, IDX, f=0.5, mode='cluster', start=-1)
    out5 = distribute_idx(MOL, IDX, f=0.5, mode='random')
    out6 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', follow_edge=True)
    out7 = distribute_idx(MOL, IDX, f=0.5, mode='cluster', follow_edge=True)
    out8 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', cluster_size=2)
    out9 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', cluster_size=range(1, 4))
    out10 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', randomness=0.5)

    out11 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', weight=lambda x: x**-2)
    out12 = distribute_idx(MOL, IDX, f=0.5, mode='cluster', weight=lambda x: x**-2)
    out13 = distribute_idx(MOL, IDX, f=0.5, mode='uniform', weight=lambda x: np.log(x)**-1)
    out14 = distribute_idx(MOL, IDX, f=0.5, mode='cluster', weight=lambda x: np.log(x)**-1)

    np.testing.assert_array_equal(out1, [127, 130, 138, 126, 148, 144, 131, 123, 146, 136, 124, 134, 143])  # noqa
    np.testing.assert_array_equal(out2, [130, 133, 124, 141, 123, 132, 147, 126, 137, 129, 148, 142, 136])  # noqa
    np.testing.assert_array_equal(out3, [148, 138, 135, 141, 134, 137, 123, 146, 144, 136, 139, 127, 133])  # noqa
    np.testing.assert_array_equal(out4, [148, 129, 142, 136, 123, 132, 147, 126, 137, 130, 133, 124, 141])  # noqa
    np.testing.assert_array_equal(out6, [127, 130, 138, 126, 148, 144, 131, 132, 146, 136, 124, 134, 143])  # noqa
    np.testing.assert_array_equal(out7, [130, 133, 124, 141, 123, 132, 147, 137, 126, 129, 148, 136, 142])  # noqa
    np.testing.assert_array_equal(out8, [127, 123, 146, 128, 124, 141, 126, 137, 145, 144, 148, 129, 139])  # noqa
    np.testing.assert_array_equal(out9, [127, 130, 141, 126, 147, 137, 140, 148, 129, 144, 145, 125, 139])  # noqa

    np.testing.assert_array_equal(out11, [127, 130, 138, 147, 144, 148, 131, 146, 137, 136, 124, 134, 123])  # noqa
    np.testing.assert_array_equal(out12, [132, 147, 126, 137, 123, 130, 133, 141, 124, 129, 148, 136, 142])  # noqa
    np.testing.assert_array_equal(out13, [127, 130, 138, 147, 144, 148, 131, 146, 137, 124, 136, 134, 123])  # noqa
    np.testing.assert_array_equal(out14, [132, 147, 137, 126, 123, 130, 133, 141, 124, 129, 148, 136, 142])  # noqa

    assertion.len_eq(out5, round(0.5 * len(IDX)))
    assertion.len_eq(np.intersect1d(out5, IDX), len(out5))
    assertion.len_eq(out10, round(0.5 * len(IDX)))
    assertion.len_eq(np.intersect1d(out10, IDX), len(out10))

    assertion.assert_(distribute_idx, MOL, IDX, f=2, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=-1, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, mode='bob', exception=ValueError)
    assertion.assert_(distribute_idx, MOL, ['bob'], f=0.5, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, cluster_size=0, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, cluster_size='bob', exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, cluster_size=['bob'], exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, start='bob', exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, start=9999, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, randomness='bob', exception=TypeError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, randomness=-10, exception=ValueError)
    assertion.assert_(distribute_idx, MOL, IDX, f=0.5, randomness=1.5, exception=ValueError)


@pytest.mark.slow
@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_cat() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT_subset.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

    arg.path = PATH
    qd_df, _, _ = prep(arg)

    assertion.len_eq(qd_df, 1)
    qd = qd_df[MOL_KEY].iloc[0]
    qd_ref = readpdb(PATH / "qd_test_distribution.pdb")
    assert_mol_allclose(qd, qd_ref)
