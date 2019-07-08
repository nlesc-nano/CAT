"""
CAT.assertion_functions
=======================

Various generic assertion functions for the testing of CAT.

Index
-----
.. currentmodule:: CAT.assertion_functions
.. autosummary::
    invert
    assert_value
    assert_id
    assert_instance
    assert_exception
    assert_isin
    assert_lt
    assert_gt
    _err_msg
    _exc_msg

API
---
.. autoclass:: CAT.assertion_functions.Invert
.. autofunction:: CAT.assertion_functions.assert_value
.. autofunction:: CAT.assertion_functions.assert_id
.. autofunction:: CAT.assertion_functions.assert_instance
.. autofunction:: CAT.assertion_functions.assert_exception
.. autofunction:: CAT.assertion_functions.assert_isin
.. autofunction:: CAT.assertion_functions.assert_lt
.. autofunction:: CAT.assertion_functions.assert_gt
.. autofunction:: CAT.assertion_functions._err_msg
.. autofunction:: CAT.assertion_functions._exc_msg

"""

from functools import wraps
from typing import (Any, Callable, Tuple, Sequence, Container)


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
    func: |Callable|_
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

    def get_err_msg(self, func: Callable,
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
                raise AssertionError(self.get_err_msg(func, *tup))
        return wrapper


def assert_lt(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value < ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'value < reference'
    assert value < ref, _err_msg(assertion, value, ref)
    return assertion, value, ref


def assert_gt(value: Any,
              ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value > ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'value > reference'
    assert value > ref, _err_msg(assertion, value, ref)
    return assertion, value, ref


def assert_isin(value: Any,
                ref: Container) -> Tuple[str, Any, Container]:
    """Assert :code:`value in ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'value is in reference'
    assert value in ref, _err_msg(assertion, value, ref)
    return assertion, value, ref


def assert_instance(value: Any,
                    ref: type) -> Tuple[str, str, str]:
    """Assert :code:`isinstance(value, ref)`; returns arguments for :func:`._err_msg`."""
    assertion = 'isinstance(value, reference)'
    ref_name = ref.__name__
    value_name = value.__class__.__name__
    assert isinstance(value, ref), _err_msg(assertion, ref_name, value_name)
    return assertion, ref_name, value_name


def assert_value(value: Any,
                 ref: Any) -> Tuple[str, Any, Any]:
    """Assert :code:`value == ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'value == reference'
    assert value == ref, _err_msg(assertion, value, ref)
    return assertion, value, ref


def assert_id(value: Any,
              ref: Any) -> Tuple[str, int, int]:
    """Assert :code:`value is ref`; returns arguments for :func:`._err_msg`."""
    assertion = 'value is reference'
    assert_value(ref, value)
    value_id = f'{id(value)} = id({value})'
    ref_id = f'{id(ref)} = id({ref})'
    assert ref is value, _err_msg(assertion, value_id, ref_id)
    return assertion, value_id, ref_id


def assert_exception(exc: Exception,
                     func: Callable,
                     *args: Sequence,
                     **kwargs: dict) -> Tuple[str, str, str]:
    """Assert that :code:`func(*args, **kwargs)` raises **exc**."""
    err_ref = exc.__name__
    err = None.__class__.__name__
    assertion = f'{func}(*{args}, **{kwargs}) -> {exc.__name__}'

    try:
        func(*args, **kwargs)
    except exc:  # The desired exception is raised
        pass
    except Exception as ex:  # An undesired exception is raised
        err = repr(ex.__class__.__name__)
        raise AssertionError(_exc_msg(assertion, err, err_ref))
    else:  # No exception is raised; this is undesired
        raise AssertionError(_exc_msg(assertion, err, err_ref))

    _assertion = f'{func}(*{args}, **{kwargs}) -/> {exc.__name__}'
    return _assertion, err, err_ref


def _err_msg(assertion: Any,
             value: Any,
             ref: Any) -> str:
    """Return a formatted error message."""
    i, j, k = repr(assertion), repr(value), repr(ref)
    return f'{i}\nSupplied value:\n\t{j}\n\nSupplied reference:\n\t{k}'


def _exc_msg(assertion: Any,
             value: Any,
             ref: Any) -> str:
    """Return a formatted error message for :func:`.assert_exception`."""
    i, j, k = repr(assertion), repr(value), repr(ref)
    return f'{i}\nSupplied exception:\n\t{j}\n\nReference exception:\n\t{k}'
