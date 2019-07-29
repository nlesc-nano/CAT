"""
CAT.assertion_functions
=======================

Various generic assertion functions for the testing of CAT.

Index
-----
.. currentmodule:: CAT.assertion_functions
.. autosummary::
    Invert
    assert_hasattr
    assert_isfile
    assert_isdir
    assert_len
    assert_eq
    assert_id
    assert_subclass
    assert_instance
    assert_exception
    assert_isin
    assert_lt
    assert_le
    assert_gt
    assert_ge
    _err_msg
    _exc_msg

API
---
.. autoclass:: Invert
.. autofunction:: assert_hasattr
.. autofunction:: assert_isfile
.. autofunction:: assert_isdir
.. autofunction:: assert_len
.. autofunction:: assert_eq
.. autofunction:: assert_id
.. autofunction:: assert_subclass
.. autofunction:: assert_instance
.. autofunction:: assert_exception
.. autofunction:: assert_isin
.. autofunction:: assert_lt
.. autofunction:: assert_le
.. autofunction:: assert_gt
.. autofunction:: assert_ge
.. autofunction:: _err_msg
.. autofunction:: _exc_msg

"""

from typing import (Any, Callable, Tuple, Container, Sized)
from reprlib import Repr
from os.path import (isfile, isdir)
from textwrap import TextWrapper
from functools import wraps


class Invert():
    """Context manager for inverting assertion result.

    Instances of :exc:`AssertionError` raised by the passed callables are supressed and
    *vice versa*.

    Examples
    --------
    .. code:: python

        >>> def assert_true(value):
        >>>     assert value is True, repr("value is not 'True'")

        >>> assert_true(False)
        AssertionError: "value is not 'True'"

        >>> with Invert(assert_true) as func:
        >>>     func(False)  # Raises no exception

    Parameters
    ----------
    func : |Callable|_
        A callable that can raise an :exc:`AssertionError`.

    Attributes
    ----------
    func : |Callable|_
        A callable constructed from the **func** parameter.
        Operations that previously did *not* raise an :exc:`AssertionError` now do and *vice versa*.

    """

    def __init__(self, func: Callable) -> None:
        """Initialize the :class:`.Invert` context manager."""
        self.func = self.invert(func)

    def __enter__(self) -> Callable:
        """Return the inverted assertion function."""
        return self.func

    def __exit__(self, *args) -> None:
        """Close the :class:`.Invert` context manager."""
        return

    @staticmethod
    def get_err_msg(func: Callable,
                    args: Tuple[str, Any, Any]) -> str:
        """Create an error message for :meth:`Invert.invert`."""
        if not args:
            return ''
        elif func is assert_exception:
            return _exc_msg(*args)
        else:
            return _err_msg(*args)

    def invert(self, func: Callable) -> Callable:
        """Invert a function that may or may not raise an :exc:`AssertionError`."""
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                tup = func(*args, **kwargs)
            except AssertionError:
                pass
            else:
                raise AssertionError(self.get_err_msg(func, tup))
        return wrapper


def assert_hasattr(value: str, ref: Any) -> Tuple[str, str, str]:
    """Assert :code:`hasattr(ref, value); returns arguments for :func:`._err_msg`."""
    assertion = 'assert hasattr(ref, value)'
    if callable(ref):
        _ref = str(ref)
    else:
        _ref = f'<{ref.__class__.__module__}.{ref.__class__.__name__} at {hex(id(ref))}>'
    assert hasattr(ref, value), _err_msg(assertion, value, _ref)

    _assertion = 'assert not hasattr(ref, value)'
    return _assertion, value, _ref


def assert_isfile(value: str) -> Tuple[str, str, None]:
    """Assert :code:`os.path.isfile(value)`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert os.path.isfile(value)'
    assert isfile(value), _err_msg(assertion, value, None)

    _assertion = 'assert not os.path.isfile(value)'
    return _assertion, value, None


def assert_isdir(value: str) -> Tuple[str, str, None]:
    """Assert :code:`os.path.isdir(value)`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert os.path.isdir(value)'
    assert isdir(value), _err_msg(assertion, value, None)

    _assertion = 'assert not os.path.isdir(value)'
    return _assertion, value, None


def assert_len(value: Sized,
               ref: int) -> Tuple[str, Any, Any]:
    """Assert :code:`len(value) == ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert len(value) == reference'
    assert len(value) == ref, _err_msg(assertion, value, ref)

    _assertion = 'assert len(value) != reference'
    return _assertion, value, ref


def assert_le(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value <= ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value <= reference'
    assert value <= ref, _err_msg(assertion, value, ref)

    _assertion = 'assert value > reference'
    return _assertion, value, ref


def assert_ge(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value >= ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value => reference'
    assert value >= ref, _err_msg(assertion, value, ref)

    _assertion = 'assert value < reference'
    return _assertion, value, ref


def assert_lt(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value < ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value < reference'
    assert value < ref, _err_msg(assertion, value, ref)

    _assertion = 'assert value >= reference'
    return _assertion, value, ref


def assert_gt(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value > ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value > reference'
    assert value > ref, _err_msg(assertion, value, ref)

    _assertion = 'assert value <= reference'
    return _assertion, value, ref


def assert_isin(value: Any,
                ref: Container) -> Tuple[str, Any, Container]:
    """Assert :code:`value in ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value in reference'
    assert value in ref, _err_msg(assertion, value, ref)

    _assertion = 'assert value not in reference'
    return _assertion, value, ref


def assert_instance(value: Any,
                    ref: type) -> Tuple[str, str, str]:
    """Assert :code:`isinstance(value, ref)`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert isinstance(value, reference)'
    ref_name = ref.__name__
    value_name = value.__class__.__name__
    assert isinstance(value, ref), _err_msg(assertion, ref_name, value_name)

    _assertion = 'assert not isinstance(value, reference)'
    return _assertion, ref_name, value_name


def assert_subclass(value: type,
                    ref: type) -> Tuple[str, str, str]:
    """Assert :code:`issubclass(value, ref)`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert issubclass(value, reference)'
    ref_name = ref.__name__
    value_name = value.__name__
    assert issubclass(value, ref), _err_msg(assertion, ref_name, value_name)

    _assertion = 'assert not issubclass(value, reference)'
    return _assertion, ref_name, value_name


def assert_eq(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value == ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value == reference'
    assert value == ref, _err_msg(assertion, value, ref)

    _assertion = 'assert value != reference'
    return _assertion, value, ref


def _str(value: Any) -> str:
    args = value.__class__.__module__, value.__class__.__name__, hex(id(value))
    return '<{}.{} at {}>'.format(*args)


def assert_id(value: Any,
              ref: Any) -> Tuple[str, str, str]:
    """Assert :code:`value is ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'assert value is reference'
    value_id = f'{_str(value)}'
    ref_id = f'{_str(ref)}'
    assert ref is value, _err_msg(assertion, value_id, ref_id)

    _assertion = 'assert value is not reference'
    return _assertion, value_id, ref_id


def assert_exception(exc: Callable[..., Exception],
                     func: Callable[..., Any],
                     *args: Any,
                     **kwargs: Any) -> Tuple[str, str, str]:
    """Assert **exc** is raised by :code:`func(*args, **kwargs)`."""
    err_ref = exc.__name__
    err = 'None'

    # Prepare arguments (as str) for **func**
    arguments = ''
    if args is not None:
        arguments += ', '.join(repr(i) for i in args)
    if kwargs is not None:
        arguments += ', '.join(f'{k}={v}' for k, v in kwargs.items())

    # Construct the assertion as str
    assertion = f'assert {func.__qualname__}({arguments}) -> {exc.__name__}'

    # Conduct the assertion
    try:
        func(*args, **kwargs)
    except exc:  # The desired exception is raised
        pass
    except Exception as ex:  # An undesired exception is raised
        err = ex.__class__.__name__
        raise AssertionError(_exc_msg(assertion, err, err_ref))
    else:  # No exception is raised; this is undesired
        raise AssertionError(_exc_msg(assertion, err, err_ref))

    _assertion = f'assert {func.__qualname__}({arguments}) -/> {exc.__name__}'
    return _assertion, err, err_ref


class FloatRepr(Repr):
    """A subclass of :class:`reprlib.Repr` with an additional method for handling floats."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.maxstring = 65

    def repr_float(self, x: float, level: int) -> str:
        """Create a :class:`str`-represenation of a :class:`float` with 4 decimals."""
        try:
            return f'{x:4.4f}'
        except ValueError:  # x is somehow not a float; convert into a string and try again.
            return self.repr(str(x))


WRAPPER: TextWrapper = TextWrapper(width=70, initial_indent='    ', subsequent_indent='     ')
FLOATREPR: FloatRepr = FloatRepr()


def _err_msg(assertion: Any,
             value: Any,
             ref: Any) -> str:
    """Return a formatted error message."""
    args = repr(assertion), WRAPPER.fill(FLOATREPR.repr(value)), WRAPPER.fill(FLOATREPR.repr(ref))
    return '{}\nSupplied value:\n{}\n\nSupplied reference:\n{}'.format(*args)


def _exc_msg(assertion: Any,
             value: Any,
             ref: Any) -> str:
    """Return a formatted error message for :func:`.assert_exception`."""
    args = repr(assertion), WRAPPER.fill(FLOATREPR.repr(value)), WRAPPER.fill(FLOATREPR.repr(ref))
    return '{}\nObserved exception:\n{}\n\nExpected exception:\n{}'.format(*args)
