import warnings
from types import MappingProxyType
from typing import Mapping, Callable, NoReturn, Union

__all__ = ['WARN_MAP']


def _warn(exc: Exception) -> None:
    """Perform a warning using **exc**."""
    warnings.warn(str(exc), category=RuntimeWarning, stacklevel=2)


def _raise(exc: Exception) -> NoReturn:
    """Raise **exc**."""
    raise exc


def _ignore(exc: Exception) -> None:
    """Do nothing."""
    return None


WARN_MAP: Mapping[str, Callable[[Exception], Union[NoReturn, None]]] = MappingProxyType({
    'raise': _raise,
    'warn': _warn,
    'ignore': _ignore
})
