"""A module for holding the :class:`.PropertiesDataFrame` and :class:`.PropertiesSeries` classes."""

from __future__ import annotations

from typing import Optional

import pandas as pd

from scm.plams import Settings

__all__ = ['PropertiesSeries', 'PropertiesDataFrame']


class PropertiesSeries(pd.Series):
    """A subclass of the Pandas Series with an additional :attr:`.properties` attribute.

    Parameters
    ----------
    properties : dict
        Optional: A dictionary with additional user-defined (meta-)data.
        See :attr:`.PropertiesSeries.properties`.

    Attributes
    ----------
    properties : |plams.Settings|_
        A dictionary with additional user-defined (meta-)data.

    """

    _metadata = ['properties']

    def __init__(self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False,
                 properties: Optional[dict] = None) -> None:
        """Initialize the :class:`.PropertiesSeries` construction."""
        self.properties = self._sanitize_properties(properties)
        super().__init__(data, index, dtype, name, copy, fastpath)

    @property
    def _constructor(self) -> PropertiesSeries:
        """Construct a :class:`.PropertiesSeries` instance."""
        def _series(*args, **kwargs) -> PropertiesSeries:
            return PropertiesSeries(*args, **kwargs).__finalize__(self)
        return _series

    @property
    def _constructor_expanddim(self) -> PropertiesDataFrame:
        """Construct a :class:`.PropertiesDataFrame` instance."""
        def _df(*args, **kwargs) -> PropertiesDataFrame:
            return PropertiesDataFrame(*args, **kwargs).__finalize__(self)
        return _df

    @staticmethod
    def _sanitize_properties(properties: Optional[dict]) -> Settings:
        """Sanitize the **properties** parameter for :attr:`.properties`."""
        if properties is None:
            return Settings()
        elif isinstance(properties, dict):
            return Settings(properties)
        else:
            err = "The properties argument expects an instance of 'dict', not '{}'"
            raise TypeError(err.format(properties.__class__.__name__))


class PropertiesDataFrame(pd.DataFrame):
    """A subclass of the Pandas DataFrame with an additional :attr:`.properties` attribute.

    Parameters
    ----------
    properties : dict
        Optional: A dictionary with additional user-defined (meta-)data.
        See :attr:`.PropertiesDataFrame.properties`.

    Attributes
    ----------
    properties : |plams.Settings|_
        A dictionary with additional user-defined (meta-)data.

    """

    _metadata = ['properties']

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False,
                 properties: Optional[dict] = None) -> None:
        """Initialize the :class:`.PropertiesDataFrame` construction."""
        self.properties = self._sanitize_properties(properties)
        super().__init__(data, index, columns, dtype, copy)

    @property
    def _constructor(self) -> PropertiesDataFrame:
        """Construct a :class:`.PropertiesDataFrame` instance."""
        def _df(*args, **kwargs) -> PropertiesDataFrame:
            return PropertiesDataFrame(*args, **kwargs).__finalize__(self)
        return _df

    @property
    def _constructor_sliced(self) -> PropertiesSeries:
        """Construct a :class:`.PropertiesSeries` instance."""
        def _series(*args, **kwargs) -> PropertiesSeries:
            return PropertiesSeries(*args, **kwargs).__finalize__(self)
        return _series

    @staticmethod
    def _sanitize_properties(properties: Optional[dict]) -> Settings:
        """Sanitize the **properties** parameter for :attr:`.properties`."""
        if properties is None:
            return Settings()
        elif isinstance(properties, dict):
            return Settings(properties)
        else:
            err = "The properties argument expects an instance of 'dict', not '{}'"
            raise TypeError(err.format(properties.__class__.__name__))
