"""Tests for :mod:`CAT.attachment.perp_surface`."""

from pathlib import Path

import numpy as np

from scm.plams import Molecule

from CAT.attachment import get_surface_vec

PATH = Path('tests') / 'test_files'
MOL = Molecule(str(PATH / 'core' / 'Cd68Se55.xyz'))
ANCHOR = MOL.as_array(atom_subset=[at for at in MOL if at.symbol == 'Cl'])
ANCHOR.setflags(write=False)


def test_get_surface_vec() -> None:
    """Test :func:`CAT.attachment.perp_surface.get_surface_vec`."""
    vec1 = get_surface_vec(MOL)
    vec2 = get_surface_vec(MOL, ANCHOR)

    ref1 = np.load(PATH / 'get_surface_vec1.npy')
    ref2 = np.load(PATH / 'get_surface_vec2.npy')

    np.testing.assert_allclose(vec1, ref1)
    np.testing.assert_allclose(vec2, ref2)
