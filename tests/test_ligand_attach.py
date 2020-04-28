"""Tests for :mod:`CAT.attachment.ligand_attach`."""

from os.path import join

import numpy as np

from assertionlib import assertion

from CAT.attachment.ligand_attach import (_get_rotmat1, _get_rotmat2)

PATH = join('tests', 'test_files')


def test_get_rotmat1() -> None:
    """Test :func:`CAT.attachment.ligand_attach._get_rotmat1`."""
    vec1 = np.array([1, 0, 0], dtype=float)
    vec2 = np.array([0, 1, 0], dtype=float)
    vec3 = np.array([[0, 1, 0], [0, 0, 1]], dtype=float)
    vec4 = np.zeros(3, dtype=float)

    ref1 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], ndmin=3, dtype=float)
    _ref2 = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]], ndmin=3, dtype=float)
    ref2 = np.vstack([ref1, _ref2])
    ref3 = np.zeros((2, 3, 3))
    ref3[:] = np.identity(3)
    ref4 = np.array([[[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                     [[0, 0, -1], [0, 1, 0], [1, 0, 0]]], dtype=float)
    ref5 = np.identity(3)[None, ...]

    rotmat1 = _get_rotmat1(vec1, vec2)
    np.testing.assert_allclose(rotmat1, ref1)

    rotmat2 = _get_rotmat1(vec1, vec3)
    np.testing.assert_allclose(rotmat2, ref2)

    rotmat3 = _get_rotmat1(vec3, vec3)
    np.testing.assert_allclose(rotmat3, ref3)

    rotmat4 = _get_rotmat1(vec3, vec1)
    np.testing.assert_allclose(rotmat4, ref4)

    rotmat5 = _get_rotmat1(vec4, vec1)
    np.testing.assert_allclose(rotmat5, ref5)

    assertion.assert_(_get_rotmat1, vec3[..., None], vec3, exception=ValueError)


def test_get_rotmat2() -> None:
    """Test :func:`CAT.attachment.ligand_attach._get_rotmat2`."""
    vec1 = np.array([1, 0, 0], dtype=float)
    vec2 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)

    ref1 = np.load(join(PATH, 'rotmat2_1.npy'))
    ref2 = np.load(join(PATH, 'rotmat2_2.npy'))

    rotmat1 = _get_rotmat2(vec1)
    np.testing.assert_allclose(rotmat1, ref1)

    rotmat2 = _get_rotmat2(vec2)
    np.testing.assert_allclose(rotmat2, ref2)
