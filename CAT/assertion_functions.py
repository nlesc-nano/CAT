"""
CAT.assertion_functions
=======================

Various generic assertion functions for the testing of CAT.

Index
-----
.. currentmodule:: CAT.assertion_functions
.. autosummary::
    assert_value
    assert_id
    assert_instance
    assert_exception
    assert_isin
    assert_lt
    assert_gt


API
---
.. autofunction:: CAT.assertion_functions.assert_value
.. autofunction:: CAT.assertion_functions.assert_id
.. autofunction:: CAT.assertion_functions.assert_instance
.. autofunction:: CAT.assertion_functions.assert_exception
.. autofunction:: CAT.assertion_functions.assert_isin
.. autofunction:: CAT.assertion_functions.assert_lt
.. autofunction:: CAT.assertion_functions.assert_gt


"""

from functools import wraps
from typing import (Hashable, Any, Callable, Tuple, Sequence, Container)


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

    err_dict : |dict|_ [|Callable|_, |str|_]
        A dictionary with error messages for specific functions:
            * :func:`.assert_value`
            * :func:`.assert_id`
            * :func:`.assert_instance`
            * :func:`.assert_exception`
            * :func:`.assert_isin`
            * :func:`.assert_lt`
            * :func:`.assert_gt`

    """

    def __init__(self, func: Callable) -> None:
        """Initialize the :class:`.Invert` context manager."""
        self.err_dict = {
            assert_value: '\ninvalid value observed:\t{}',
            assert_id: '\ninvalid id observed:\t{}',
            assert_instance: '\ninvalid instance observed:\t{}',
            assert_exception: '\ninvalid exception observed:\t{}',
            assert_isin: '\nvalue:\t{}\n\nis in:\t{}',
            assert_lt: '\nvalue:\t{}\n\nis smaller than:\t{}',
            assert_gt: '\nvalue\t{}\n\nis larger than:\t{}',
        }
        self.func = self.invert(func)

    def __enter__(self) -> Callable:
        """Return the inverted assertion function."""
        return self.func

    def __exit__(self, *args) -> None:
        """Close the :class:`.Invert` context manager."""
        pass

    def invert(self, func: Callable) -> Callable:
        """Invert a function that may or may not raise an :exc:`AssertionError`."""
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                i = func(*args, **kwargs)
            except AssertionError:
                pass
            else:
                raise AssertionError(self.err_dict.setdefault(func, '').format(*i))
        return wrapper


def assert_lt(value: Any, ref) -> Tuple[str, str]:
    """Assert :code:`value < ref`; returns a tuple of **value** and **ref**."""
    assert value < ref, f'\nvalue:\t{repr(value)}\n\nis larger than:\t{repr(value)}'
    return repr(value), repr(ref)


def assert_gt(value: Any, ref) -> Tuple[str, str]:
    """Assert :code:`value > ref`; returns a tuple of **value** and **ref**."""
    assert value > ref, f'\nvalue\t{repr(value)}\n\nis smaller than:\t{repr(value)}'
    return repr(value), repr(ref)


def assert_isin(value: Any, ref: Container) -> str:
    """Assert :code:`value in ref`; returns **value** as :class:`str`."""
    assert value in ref, f'\nvalue\t{repr(value)}\n\nis not in:\t{repr(ref)}'
    return repr(value)


def assert_instance(value: Any, ref: Callable) -> str:
    """Assert :code:`isinstance(value, ref)`; returns the class name of **value**."""
    ref_name = repr(ref.__name__)
    value_name = repr(value.__class__.__name__)
    err = f'\ninstance expected:\t{ref_name}; \n\ninstance observed:\t{value_name}'
    assert isinstance(value, ref), err
    return value_name


def assert_value(value: Hashable,
                 ref: Hashable) -> str:
    """Assert :code:`value == ref`; returns **value** as :class:`str`."""
    assert value == ref, f'\nvalue expected:\t{repr(ref)}; \n\nvalue observed:\t{repr(value)}'
    return repr(value)


def assert_id(value: Hashable,
              ref: Hashable) -> str:
    """Assert :code:`value is ref`; returns the :class:`id` of **value**."""
    assert_value(ref, value)
    assert ref is value, f'\nid expected:\t{repr(id(ref))}; \n\nid observed:\t{repr(id(value))}'
    return repr(id(value))


def assert_exception(exc: Exception,
                     func: Callable,
                     *args: Sequence,
                     **kwargs: dict) -> str:
    """Assert that :code:`func(*args, **kwargs)` raises **exc**; returns the name of **exc**."""
    err_ref = repr(exc.__name__)
    err = repr(None.__class__.__name__)

    try:
        func(*args, **kwargs)
    except exc:  # The desired exception is raised
        pass
    except Exception as ex:  # An undesired exception is raised
        err = repr(ex.__class__.__name__)
        raise AssertionError(f'\nexception expected:\t{err_ref}; \n\nexception observed:\t{err}')
    else:  # No exception is raised; this is undesired
        raise AssertionError(f'\nexception expected:\t{err_ref}; \n\nexception observed:\t{err}')
    return err_ref
