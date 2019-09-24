"""
CAT.nd_repr
===========

A module for holding the :class:`NDRepr` class.

Index
-----
.. currentmodule:: CAT.assertion_functions
.. autosummary::
    NDRepr

API
---
.. autoclass:: NDRepr
    :members:

"""

import reprlib
from typing import Any, Dict, Union
from textwrap import indent
from itertools import chain

try:
    import numpy as np
    ndarray: Union[type, str] = np.ndarray
except ImportError:
    ndarray: Union[type, str] = 'numpy.ndarray'

try:
    import pandas as pd
    DataFrame: Union[type, str] = pd.DataFrame
    Series: Union[type, str] = pd.Series
except ImportError:
    DataFrame: Union[type, str] = 'pandas.DataFrame'
    Series: Union[type, str] = 'pandas.Series'


class NDRepr(reprlib.Repr):
    """A subclass of :class:`reprlib.Repr`.

    Has additional methods for handling instances of:

    * :class:`float` : :meth:`NDRepr.repr_float`
    * :class:`numpy.ndarray` : :meth:`NDRepr.repr_ndarray`
    * :class:`pandas.Series` : :meth:`NDRepr.repr_Series`
    * :class:`pandas.DataFrame` : :meth:`NDRepr.repr_DataFrame`

    Parameters
    ----------
    **kwargs : object
        User-specified values for one or more :class:`NDRepr` instance attributes.
        An :exc:`AttributeError` is raised upon encountering unrecognized keys.

    Attributes
    ----------
    maxfloat : :class:`int`
        The maximum number of to-be displayed :class:`float` decimals.

    maxndarray : :class:`int`
        The maximum number of items in a :class:`numpy.ndarray` row.
        Passed as argument to the :func:`numpy.printoptions` function:

        * :code:`threshold = self.maxndarray`
        * :code:`edgeitems = self.maxndarray // 2`

    maxseries : :class:`int`
        The maximum number of rows per :class:`pandas.Series` instance.
        Passed as value to :attr:`pandas.options.display`.

        * :code:`pandas.options.display.max_rows = self.series`

    maxdataframe : :class:`int`
        The maximum number of rows per :class:`pandas.DataFrame` instance.
        Passed as values to :attr:`pandas.options.display`:

        * :code:`pandas.options.display.max_rows = self.maxdataframe`
        * :code:`pandas.options.display.max_columns = self.maxdataframe // 2`

    np_printoptions : :class:`dict`
        Additional keyword arguments for :func:`numpy.printoptions`.

        .. note::
            Arguments provided herein will take priority over those specified internally
            in :meth:`NDRepr.repr_ndarray`.

    pd_printoptions : :class:`dict`
        Additional "keyword arguments" for :attr:`pandas.options`.

        .. note::
            Arguments provided herein will take priority over those specified internally
            in :meth:`NDRepr.repr_DataFrame` and :meth:`NDRepr.repr_Series`.

    """

    def __init__(self, **kwargs) -> None:
        """Initialize a :class:`NDRepr` instance."""
        super().__init__()
        self.maxfloat: int = 4
        self.maxndarray: int = 6
        self.maxseries: int = 12
        self.maxdataframe: int = 12
        self.np_printoptions: Dict[str, Any] = {}
        self.pd_printoptions: Dict[str, Any] = {}

        # Update attributes based on **kwargs; raise an error if a key is unrecognized
        for k, v in kwargs.items():
            if not hasattr(self, k):
                err = f'{repr(self.__class__.__name__)} object has no attribute {repr(k)}'
                raise AttributeError(err)
            setattr(self, k, v)

    """################################ Non-magic methods ################################"""

    def repr_float(self, obj: float, level: int) -> str:
        """Create a :class:`str` represenation of a :class:`float` with 4 decimals."""
        i = self.maxfloat
        return f'{obj:{i}.{i}f}'

    def repr_ndarray(self, obj: ndarray, level: int) -> str:
        """Create a :class:`str` represenation of a :class:`numpy.ndarray` instance."""
        if level <= 0:
            return f'{obj.__class__.__name__}(...)'

        kwargs = {'threshold': self.maxndarray,
                  'edgeitems': self.maxndarray // 2,
                  'formatter': self._get_ndformatter(obj)}
        kwargs.update(self.np_printoptions)

        with np.printoptions(**kwargs):
            return repr(obj)

    def repr_DataFrame(self, obj: DataFrame, level: int) -> str:
        """Create a :class:`str` represenation of a :class:`pandas.DataFrame` instance."""
        if level <= 0:
            return f'{obj.__class__.__name__}(...)'

        kwargs = {'display.max_rows': self.maxdataframe,
                  'display.max_columns': self.maxdataframe // 2}
        kwargs.update(self.pd_printoptions)
        args = chain.from_iterable(kwargs.items())

        with pd.option_context(*args):
            return repr(obj)

    def repr_Series(self, obj: Series, level: int) -> str:
        """Create a :class:`str` represenation of a :class:`pandas.Series` instance."""
        if level <= 0:
            return f'{obj.__class__.__name__}(...)'

        kwargs = {'display.max_rows': self.maxseries}
        kwargs.update(self.pd_printoptions)
        args = chain.from_iterable(kwargs.items())

        with pd.option_context(*args):
            return repr(obj)

    def _get_ndformatter(self, obj: ndarray) -> dict:
        """Return a value for the **formatter** argument in :func:`numpy.printoptions`."""
        if obj.dtype not in (np.dtype(float), np.dtype(int)):
            return {}

        max_len = len(str(int(obj.max())))
        min_len = len(str(int(obj.min())))
        width = max(max_len, min_len)

        if obj.dtype == np.dtype(float):
            width += 5
            value = '{' + f':{width}.{self.maxfloat}f' + '}'
            return {'float': value.format}
        else:  # obj.dtype == np.dtype(int)
            value = '{' + f':{width}d' + '}'
            return {'int': value.format}


#: An instance of :class:`NDRepr`.
aNDRepr: NDRepr = NDRepr()
