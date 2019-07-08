"""Tests for :mod:`CAT.frozen_settings`."""

from CAT.frozen_settings import FrozenSettings
from CAT.assertion_functions import (assert_value, assert_exception, assert_id, Invert)

SETTINGS = FrozenSettings({'a': True, 'b': False, 'c': [1, 2, 3, 4]})


def test_missing() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.__missing__`."""
    item = SETTINGS.d
    assert_value(item, FrozenSettings())


def test_delitem() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.__delitem__`."""
    args = ('a')
    assert_exception(TypeError, SETTINGS.__delitem__, args)


def test_setitem() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.__setitem__`."""
    args = ('d', True)
    assert_exception(TypeError, SETTINGS.__setitem__, args)


def test_copy() -> None:
    """Tests for :meth:`CAT.frozen_settings.FrozenSettings.copy`."""
    settings = SETTINGS.copy()
    assert_value(settings, SETTINGS)
    with Invert(assert_id) as func:
        func(settings, SETTINGS)
