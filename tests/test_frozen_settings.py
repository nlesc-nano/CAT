"""Tests for :mod:`CAT.frozen_settings`."""

from assertionlib import assertion

from CAT.frozen_settings import FrozenSettings

SETTINGS = FrozenSettings({'a': True, 'b': False, 'c': [1, 2, 3, 4]})


def test_missing() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.__missing__`."""
    item = SETTINGS.d
    assertion.eq(item, FrozenSettings())


def test_delitem() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.__delitem__`."""
    args = ('a')
    assertion.assert_(SETTINGS.__delitem__, args, exception=TypeError)


def test_setitem() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.__setitem__`."""
    args = ('d', True)
    assertion.assert_(SETTINGS.__setitem__, args, exception=TypeError)


def test_copy() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.copy`."""
    settings = SETTINGS.copy()
    assertion.eq(settings, SETTINGS)
    assertion.is_not(settings, SETTINGS)
