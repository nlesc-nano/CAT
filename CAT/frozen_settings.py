"""
CAT.frozen_settings
===================

A module which adds the :class:`.FrozenSettings` class, an immutable counterpart to plams.Settings_.

.. _plams.Settings: https://www.scm.com/doc/plams/components/settings.html

Index
-----
.. currentmodule:: CAT.frozen_settings
.. autosummary::
    FrozenSettings

API
---
.. autoclass:: FrozenSettings
    :members:
    :private-members:
    :special-members:

"""

from __future__ import annotations

from typing import (Any, Union, Iterable)

from scm.plams import Settings

__all__ = ['FrozenSettings']

# Various immutable objects suited as dictionary keys
Immutable = Union[tuple, int, float, str, frozenset]


class FrozenSettings(Settings):
    """An inmutable subclass of plams.Settings_.

    .. _plams.Settings: https://www.scm.com/doc/plams/components/settings.html
    """

    def __init__(self, *args: Iterable, **kwargs: dict) -> None:
        """Initialize the construction of a :class:`FrozenSettings` instance."""
        dict.__init__(self, *args, **kwargs)

        # Fill the FrozenSettings instance by means of the dict.__setitem__ method
        for key, value in self.items():
            if isinstance(value, dict):
                Settings.__setitem__(self, key, FrozenSettings(value))
            elif isinstance(value, list):
                value = [FrozenSettings(i) if isinstance(i, dict) else i for i in value]
                Settings.__setitem__(self, key, value)

    def __missing__(self, name: Immutable) -> FrozenSettings:
        """Return a new (empty) :class:`FrozenSettings` instance."""
        return FrozenSettings()

    def __delitem__(self, name: Immutable) -> None:
        """Raise a :exc:`TypeError`; :class:`FrozenSettings` instances are immutable."""
        raise TypeError("'FrozenSettings' object does not support item deletion")

    def __setitem__(self, name: Immutable, value: Any) -> None:
        """Raise a :exc:`TypeError`; :class:`FrozenSettings` instances are immutable."""
        raise TypeError("'FrozenSettings' object does not support item assignment")

    def copy(self) -> FrozenSettings:
        """Create a copy of this instance."""
        ret = FrozenSettings()
        for key, value in self.items():
            try:
                Settings.__setitem__(ret, key, value.copy())
            except AttributeError:
                Settings.__setitem__(ret, key, value)
        return ret

    def __copy__(self) -> FrozenSettings:
        """Create a copy of this instance by calling :meth:`FrozenSettings.copy`."""
        return self.copy()
