"""A module for containing the :class:`SetAttr` class."""

import reprlib
from threading import RLock
from typing import Generic, TypeVar, NoReturn, Dict, Any, Optional

__all__ = ['SetAttr']

T1 = TypeVar('T1')
T2 = TypeVar('T2')
ST = TypeVar('ST', bound='SetAttr')


class SetAttr(Generic[T1, T2]):
    """A context manager for temporarily changing an attribute's value.

    The :class:`SetAttr` context manager is thread-safe, reusable and reentrant.

    Wanrings
    --------
    Note that while :meth:`SetAttr.__enter__` and :meth:`SetAttr.__exit__` are thread-safe,
    the same does *not* hold for :meth:`SetAttr.__init__`.

    Examples
    --------
    .. code:: python

        >>> from CAT.utils import SetAttr

        >>> class Test:
        ...     a = False

        >>> print(Test.a)
        False

        >>> set_attr = SetAttr(Test, 'a', True)
        >>> with set_attr:
        ...     print(Test.a)
        True

    """

    __slots__ = ('__weakref__', '_obj', '_name', '_value', '_value_old', '_lock', '_hash')

    @property
    def obj(self) -> T1:
        """The to-be modified object."""
        return self._obj

    @property
    def name(self) -> str:
        """The name of the to-be modified attribute."""
        return self._name

    @property
    def value(self) -> T2:
        """The value to-be assigned to the :attr:`name` attribute of :attr:`obj`."""
        return self._value

    @property
    def attr(self) -> T2:
        """Get or set the :attr:`~SetAttr.name` attribute of :attr:`SetAttr.obj`."""
        return getattr(self.obj, self.name)

    @attr.setter
    def attr(self, value: T2) -> None:
        with self._lock:
            setattr(self.obj, self.name, value)

    def __init__(self, obj: T1, name: str, value: T2) -> None:
        """Initialize the :class:`SetAttr` context manager.

        Parameters
        ----------
        obj : :data:`~typing.Any`
            The to-be modified object.
        name : :class:`str`
            The name of the to-be modified attribute.
        value : :data:`~typing.Any`
            The value to-be assigned to the **name** attribute of **obj**.

        """
        self._obj = obj
        self._name = name
        self._value = value

        self._value_old = self.attr
        self._lock = RLock()

    def __repr__(self) -> str:
        """Implement :func:`str(self)<str>` and :func:`repr(self)<repr>`."""
        obj = object.__repr__(self.obj)
        value = reprlib.repr(self.value)
        return f'{self.__class__.__name__}(obj={obj}, name={self.name!r}, value={value})'

    def __eq__(self, value: Any) -> bool:
        """Implement :func:`self == value<object.__eq__>`."""
        if type(self) is not type(value):
            return False
        return self.obj is value.obj and self.name == value.name and self.value == value.value

    def __reduce__(self) -> NoReturn:
        """Unsupported operation, raise a :exc:`TypeError`."""
        raise TypeError(f"can't pickle {self.__class__.__name__} objects")

    def __copy__(self: ST) -> ST:
        """Implement :func:`copy.copy(self)<copy.copy>`."""
        return self

    def __deepcopy__(self: ST, memo: Optional[Dict[int, Any]] = None) -> ST:
        """Implement :func:`copy.deepcopy(self, memo=...)<copy.deepcopy>`."""
        return self

    def __hash__(self) -> int:
        """Implement :func:`hash(self)<hash>`.

        Warnings
        --------
        A :exc:`TypeError` will still be raised if :attr:`SetAttr.value` is not hashable.

        """
        try:
            return self._hash
        except AttributeError:
            args = (type(self), (id(self.obj), self.name, self.value))
            self._hash: int = hash(args)
            return self._hash

    def __enter__(self) -> None:
        """Enter the context manager, modify :attr:`SetAttr.obj`."""
        self.attr = self.value

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager, restore :attr:`SetAttr.obj`."""
        self.attr = self._value_old
