"""Tests for :mod:`CAT.attachment.edge_distance`."""

from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

from scm.plams import Molecule
from assertionlib import assertion

from CAT.attachment.edge_distance import array_combinations, to_convex, edge_dist

PATH = Path('tests') / 'test_files'
MOL = Molecule(str(PATH / 'core' / 'Cd68Se55.xyz'))
XYZ = MOL.as_array(atom_subset=[at for at in MOL if at.symbol == 'Cl'])
XYZ.setflags(write=False)


def test_array_combinations() -> None:
    """Test :func:`CAT.attachment.edge_distance.array_combinations`."""
    ar = np.arange(32).reshape(8, 4)

    for i in range(1, 5):
        ref = np.load(PATH / f'combinations_r{i}.npy')
        comb = array_combinations(ar, r=i)
        np.testing.assert_array_equal(comb, ref)
    assertion.assert_(array_combinations, ar, r=-1, exception=ValueError)
    assertion.assert_(array_combinations, ar, r=5, exception=ValueError)
    assertion.assert_(array_combinations, np.random.rand(10, 10, 10), r=2, exception=ValueError)


def test_to_convex() -> None:
    """Test :func:`CAT.attachment.edge_distance.to_convex`."""
    ar = np.array([to_convex(XYZ, n=i) for i in np.arange(0.1, 1.1, 0.1)])
    ref = np.load(PATH / 'to_convex.npy')
    np.testing.assert_allclose(ar, ref)

    assertion.assert_(to_convex, XYZ, n=-1, exception=ValueError)
    assertion.assert_(to_convex, XYZ, n=0, exception=ValueError)
    assertion.assert_(to_convex, XYZ, n=1.5, exception=ValueError)


def test_edge_dist() -> None:
    """Test :func:`CAT.attachment.edge_distance.edge_dist`."""
    dist1 = edge_dist(XYZ, n=0)

    hull = ConvexHull(XYZ)
    edges = array_combinations(hull.simplices, r=2).reshape((-1, 2))
    dist2 = edge_dist(XYZ, edges=edges)

    ref1 = np.load(PATH / 'edge_dist1.npy')
    ref2 = np.load(PATH / 'edge_dist2.npy')
    np.testing.assert_allclose(dist1, ref1)
    np.testing.assert_allclose(dist2, ref2)
