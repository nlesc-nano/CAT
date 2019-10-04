"""
CAT.assertion.assertion_manager
===============================

Various assertion-related classes for the testing of CAT.

Index
-----
.. currentmodule:: CAT.assertion.assertion_manager
.. autosummary::
    AssertionManager
    _MetaAM
    _get_sphinx_domain
    _wrap_docstring
    _bind_callable

API
---
.. autoclass:: AssertionManager
    :members:
    :special-members:
    :private-members:

.. autoclass:: _MetaAM
    :members:
    :special-members:
    :private-members:

.. autofunction:: _get_sphinx_domain
.. autofunction:: _wrap_docstring
.. autofunction:: _bind_callable

"""

import os
import types
import inspect
import reprlib
import textwrap
import operator
from string import ascii_lowercase
from typing import Callable, Any, Type, FrozenSet, Optional, Union, Sized

from .ndrepr import aNDRepr
from ..abc.dataclass import AbstractDataClass

__all__ = ['AssertionManager', 'assertion']


def _get_sphinx_domain(func: Callable) -> str:
    """Create a Sphinx domain from func.

    Accepts function, classes and methods; raises a :exc:`TypeError` otherwise.

    """
    name = func.__qualname__ if hasattr(func, '__qualname__') else func.__name__
    if inspect.isbuiltin(func) or inspect.isfunction(func):
        return f':func:`{name}<{func.__module__}.{name}>`'
    elif inspect.ismethod(func):
        return f':meth:`{name}<{func.__module__}.{name}>`'
    elif inspect.isclass(func):
        return f':class:`{name}<{func.__module__}.{name}>`'
    raise TypeError(f"'{name}' is neither a (builtin) function, method nor class")


def _wrap_docstring(func: Callable) -> str:
    """Create a new NumPy style assertion docstring from the docstring of **func**."""
    domain = _get_sphinx_domain(func)

    # Extract the first line from the func docstring
    try:
        func_summary = func.__doc__.split('\n', 1)[0]
    except AttributeError:
        func_summary = 'No description.'

    # Return a new docstring
    name = func.__qualname__ if hasattr(func, '__qualname__') else func.__name__
    return ('Perform the following assertion: '
            f':code:`assert {name}(*args, **kwargs)`.\n\n    See also\n    --------\n'
            f'    {domain}:\n'
            f'        {func_summary}\n\n    ')


def _bind_callable(class_type: Union[type, Any], func: Callable,
                   name: Optional[str] = None) -> None:
    """Take a callable and use it to create a new assertion method for **class_type**.

    Parameter
    ---------
    class_type : :class:`type` or :class:`Any<typing.Any>`
        A class (*i.e.* a :class:`type` instance) or class instance.

    func : :class:`Callable<typing.Callable>`
        A callable object whose output will be asserted by the created method.

    name : :class:`str`, optional
        The name of the name of the new method.
        If ``None``, use the name of **func**.

    """
    func_name = name if name is not None else func.__name__

    # Create the to-be added method

    def method(self, *args, invert=False, **kwargs):
        __tracebackhide__ = True
        self.assert_(func, *args, invert=invert, **kwargs)

    # Update docstrings and annotations
    method.__doc__ = _wrap_docstring(func)
    method.__annotations__ = {'args': Any, 'kwargs': Any, 'invert': bool, 'return': None}

    # Set the new method
    if isinstance(class_type, type):  # A class
        setattr(class_type, func_name, method)
    else:  # A class instance
        _method = types.MethodType(method, class_type)
        setattr(class_type, func_name, _method)


def len_eq(a: Sized, b: int) -> bool:
    """Check if the length of **a** is equivalent to **b**."""
    return len(a) == b


def allclose(a: float, b: float, rtol: float = 1e-07) -> bool:
    """CHeck if the absolute differnce between **a** and **b** is smaller than **rtol**."""
    return (abs(a) - abs(b)) < rtol


class _MetaAM(type):
    """The meta-class of :class:`AssertionManager`.

    The :meth:`_MetaAM.__new__` method iterates over all (almost) functions in the :mod:`operator`
    module and creates a matching a assertion method for the :class:`AssertionManager` class.

    """
    #: A :class:`frozenset` of to-be ignored functions in :mod:`operator`.
    EXCLUDE: FrozenSet[str] = frozenset({'setitem', 'delitem'})

    #: A :class:`frozenset` of callables which need an assertion function.
    INCLUDE: FrozenSet[Callable] = frozenset({
        os.path.isfile, os.path.isdir, isinstance, issubclass, callable, hasattr, len_eq, allclose
    })

    def __new__(cls, name, bases, namespace, **kwargs) -> type:
        sub_cls = super().__new__(cls, name, bases, namespace, **kwargs)

        # Iterature over the __all__ attribute of the operator builtin module
        for name in operator.__all__:
            if name[1:] in operator.__all__ or name in _MetaAM.EXCLUDE:
                pass  # Exclude inplace operations

            func = getattr(operator, name)
            _bind_callable(sub_cls, func, name)

        for func in _MetaAM.INCLUDE:
            name = func.__name__
            _bind_callable(sub_cls, func, name)

        return sub_cls


class AssertionManager(AbstractDataClass, metaclass=_MetaAM):
    """An assertion manager.

    Paramaters
    ----------
    repr_instance : :class:`reprlib.Repr`
        An instance of :class:`reprlib.Repr` for formatting Exception messages.
        The passed instance should have access to a bound callable by the name of ``repr``,
        which in turn should produce a string representation of any passed objects.
        See also :attr:`AssertionManager.repr_instance`

    Attributes
    ----------
    repr_instance : :class:`reprlib.Repr`
        An instance of :class:`reprlib.Repr` for formatting Exception messages.
        The passed instance should have access to a bound callable by the name of ``repr``,
        which in turn should produce a string representation of passed objects.

    """

    def __init__(self, repr_instance: Type[reprlib.Repr] = aNDRepr) -> None:
        """Initialize an :class:`AssertionManager` instance."""
        self.repr_instance = repr_instance

    @property
    def repr(self) -> Callable[[Any], str]:
        """Return the :meth:`repr<reprlib.Repr.repr>` method of :attr:`AssertionManager.repr_instance`."""  # noqa
        return self.repr_instance.repr

    def _get_exc_message(self, ex: Exception, func: Callable, *args: Any,
                         invert: bool = False, **kwargs: Any) -> str:
        """Return a formatted exception message for failed assertions.

        Examples
        --------
        .. code:: python

            >>> import operator

            >>> ex = TypeError('Fancy custom exception')
            >>> func = operator.contains
            >>> a = [1, 2, 3, 4]
            >>> b = 5

            >>> assertion = AssertionManager()
            >>> msg = assertion._get_exc_message(ex, func, a, b)
            >>> raise AssertionError(msg)
            AssertionError: contains(a, b)
            Exception: TypeError('Fancy custom exception')

            Value a:
                [1, 2, 3, 4]

            Value b:
                5

        Parameters
        ----------
        ex : :class:`Exception`
            The exception raised by :meth:`AssertionManager.assert_`.

        func : :class:`Callable<typing.Callable>`
            The callable whose output has been evaluated.

        \*args : :class:`Any<typing.Any>`
            Positional arguments supplied to **func**.

        invert :class:`bool`
            If ``True``, invert the output of the assertion: :code:`not func(a, b, **kwargs)`.

        \**kwargs : :class:`Any<typing.Any>`, optional
            Further optional keyword arguments supplied to **func**.

        Returns
        -------
        :class:`str`
            A newly-formatted exception message to-be raised by :meth:`AssertionManager.assert_`.

        """  # noqa
        __tracebackhide__ = True

        indent = 4 * ' '

        name = func.__qualname__ if hasattr(func, '__qualname__') else func.__name__
        ret = f'{name}('
        for i, (_, j) in enumerate(zip(args, ascii_lowercase)):
            ret += f', {j}' if i >= 1 else j
        if kwargs:
            ret += ', **kwargs'
        ret += ')'

        if invert:
            ret = 'not ' + ret

        ret += f'\nException: {repr(ex)}'
        for i, j in zip(args, ascii_lowercase):
            ret += '\n\nValue {}:\n{}'.format(j, textwrap.indent(self.repr(i), indent))

        return ret

    def assert_(self, func: Callable, *args: Any, invert: bool = False, **kwargs: Any) -> None:
        """Perform the :func:`assert` operation on the output of :code:`func(a, b, **kwargs)`.

        Examples
        --------
        For example :code:`assert 5 == 5` is equivalent to
        :code:`AssertionManager().assert_(operator.eq, 5, 5)`.

        Parameters
        ----------
        func : :class:`Callable<typing.Callable>`
            The callable whose output will be evaluated.

        \*args : :class:`Any<typing.Any>`
            Positional arguments for **func**.

        invert :class:`bool`
            If ``True``, invert the output of the assertion: :code:`not func(a, b, **kwargs)`.

        \**kwargs : :class:`Any<typing.Any>`, optional
            Keyword arguments for **func**.

        """  # noqa
        __tracebackhide__ = True

        try:
            if invert:
                assert not func(*args, **kwargs)
            else:
                assert func(*args, **kwargs)
        except Exception as ex:
            err = self._get_exc_message(ex, func, *args, invert=invert, **kwargs)
            ex_new = AssertionError(err)
            ex_new.__traceback__ = ex.__traceback__
            raise ex_new

    def add_to_instance(self, func: Callable, override_attr: bool = False,
                        name: Optional[str] = None) -> None:
        """Add a new custom assertion method to this instance.

        Parameters
        ----------
        func : :class:`Callable<typing.Callable>`
            The callable whose output will be asserted in the to-be created method.

        override_attr : :class:`bool`
            If ``False``, raise an :exc:`AttributeError` if a method with the same name already
            exists in this instance.

        name : :class:`str`, optional
            The name of the name of the new method.
            If ``None``, use the same name as **func**.

        Exception
        ---------
        AttributeError
            Raised if ``override_attr=False`` and a method with the same name already
            exists in this instance.

        """
        meth_name = name if name is not None else func.__name__
        if not override_attr and hasattr(self, meth_name):
            raise AttributeError(f"'{self.__class__.__name__}' instance already has an attribute "
                                 f"by the name of '{meth_name}'")
        _bind_callable(self, func, name)

    def exception(self, exception: Type[Exception], func: Callable,
                  *args: Any, **kwargs: Any) -> None:
        """Assert that **exception** is raised by :code:`func(*args, **kwargs)`.

        Parameters
        ----------
        exception : class`type` [:exc:`Exception`]
            An exception that should be raised by :code:`func(*args, **kwargs)`.
            Note: :exc:`AssertionError` is dissallowed as value.

        func : :class:`Callable<typing.Callable>`
            The callable whose output will be evaluated.

        \*args : :class:`Any<typing.Any>`
            Positional arguments for **func**.

        \**kwargs : :class:`Any<typing.Any>`, optional
            Keyword arguments for **func**.

        See also
        --------
        :exc:`Exception`
            Common base class for all non-exit exceptions.

        """  # noqa
        __tracebackhide__ = True

        if exception is AssertionError:  # AssertionError
            raise TypeError("'AssertionError' is a disallowed value for the 'exception' parameter")

        try:
            func(*args, **kwargs)
            raise AssertionError(f"Failed to raise '{exception.__name__}'")
        except exception:
            pass  # This is the expected exception
        except Exception as ex:  # This is an unexpected exception
            err = self._get_exc_message(ex, func, *args, **kwargs)
            raise AssertionError(err)


#: An instance of :class:`AssertionManager`.
assertion: AssertionManager = AssertionManager()
