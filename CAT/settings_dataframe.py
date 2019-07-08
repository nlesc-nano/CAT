"""
CAT.settings_dataframe
======================

A module for holding the :class:`.SettingsDataFrame` and :class:`.SettingsSeries` classes.

Index
-----
.. currentmodule:: CAT.settings_dataframe
.. autosummary::
    SettingsSeries
    SettingsDataFrame

API
---
.. autoclass:: CAT.settings_dataframe.SettingsSeries
.. autoclass:: CAT.settings_dataframe.SettingsDataFrame

"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from scm.plams import Settings

__all__ = ['SettingsSeries', 'SettingsDataFrame']


class SettingsSeries(pd.Series):
    """A subclass of the Pandas Series with an additional :attr:`.settings` attribute.

    Parameters
    ----------
    settings : dict
        Optional: A dictionary with additional user-defined (meta-)data.
        See :attr:`.SettingsSeries.settings`.

    Attributes
    ----------
    settings : |plams.Settings|_
        A dictionary with additional user-defined (meta-)data.

    """

    _metadata = ['settings']

    def __init__(self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False,
                 settings: Optional[dict] = None) -> None:
        """Initialize the :class:`.SettingsSeries` construction."""
        self.settings = self._sanitize_settings(settings)
        super().__init__(data, index, dtype, name, copy, fastpath)

    @property
    def _constructor(self) -> SettingsSeries:
        """Construct a :class:`.SettingsSeries` instance."""
        def _series(*args, **kwargs) -> SettingsSeries:
            return SettingsSeries(*args, **kwargs).__finalize__(self)
        return _series

    @property
    def _constructor_expanddim(self) -> SettingsDataFrame:
        """Construct a :class:`.SettingsDataFrame` instance."""
        def _df(*args, **kwargs) -> SettingsDataFrame:
            return SettingsDataFrame(*args, **kwargs).__finalize__(self)
        return _df

    @staticmethod
    def _sanitize_settings(settings: Optional[dict]) -> Settings:
        """Sanitize the **settings** parameter for :attr:`.settings`."""
        if settings is None:
            return Settings()
        elif isinstance(settings, dict):
            return Settings(settings)
        else:
            err = "The settings argument expects an instance of 'dict', not '{}'"
            raise TypeError(err.format(settings.__class__.__name__))


class SettingsDataFrame(pd.DataFrame):
    """A subclass of the Pandas DataFrame with an additional :attr:`.settings` attribute.

    Parameters
    ----------
    settings : dict
        Optional: A dictionary with additional user-defined (meta-)data.
        See :attr:`.SettingsDataFrame.settings`.

    Attributes
    ----------
    settings : |plams.Settings|_
        A dictionary with additional user-defined (meta-)data.

    """

    _metadata = ['settings']

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False,
                 settings: Optional[dict] = None) -> None:
        """Initialize the :class:`.SettingsDataFrame` construction."""
        self.settings = self._sanitize_settings(settings)
        super().__init__(data, index, columns, dtype, copy)

    @property
    def _constructor(self) -> SettingsDataFrame:
        """Construct a :class:`.SettingsDataFrame` instance."""
        def _df(*args, **kwargs) -> SettingsDataFrame:
            return SettingsDataFrame(*args, **kwargs).__finalize__(self)
        return _df

    @property
    def _constructor_sliced(self) -> SettingsSeries:
        """Construct a :class:`.SettingsSeries` instance."""
        def _series(*args, **kwargs) -> SettingsSeries:
            return SettingsSeries(*args, **kwargs).__finalize__(self)
        return _series

    @staticmethod
    def _sanitize_settings(settings: Optional[dict]) -> Settings:
        """Sanitize the **settings** parameter for :attr:`.settings`."""
        if settings is None:
            return Settings()
        elif isinstance(settings, dict):
            return Settings(settings)
        else:
            err = "The settings argument expects an instance of 'dict', not '{}'"
            raise TypeError(err.format(settings.__class__.__name__))
