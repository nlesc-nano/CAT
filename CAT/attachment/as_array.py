import warnings
import reprlib
import types
from typing import Any, NoReturn, Iterable, Callable
from collections.abc import Sequence
from contextlib import AbstractContextManager

import numpy as np

from scm.plams import MoleculeError, Atom, Molecule


class AsArray(AbstractContextManager):
    r"""A context manager for temporary calling the :meth:`Molecule.as_array` method.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule

        # Create a H2 example molecule
        >>> h1 = Atom(symbol='H', coords=(0.0, 0.0, 0.0))
        >>> h2 = Atom(symbol='H', coords=(1.0, 0.0, 0.0))
        >>> mol = Molecule()
        >>> mol.add_atom(h1)
        >>> mol.add_atom(h2)

        >>> print(mol)
          Atoms:
            1         H      0.000000      0.000000      0.000000
            2         H      1.000000      0.000000      0.000000

        # Example: Translate the molecule along the Cartesian Z-axis by 5 Angstroem
        >>> with AsArray(mol) as xyz:
        ...     xyz[:, 2] += 5

        >>> print(mol)
          Atoms:
            1         H      0.000000      0.000000      5.000000
            2         H      1.000000      0.000000      5.000000

    Parameters
    ----------
    mol : |plams.Molecule|
        An iterable consisting of PLAMS atoms.
        See :attr:`AsArray.mol`.

    delete_atom : :class:`str`
        The action which is to be taken when calling the :meth:`Molecule.delete_atom` method
        of **mol**:

        * ``"raise"``: Raise a :exc:`MoleculeError`.
        * ``"warn"``: Issue a warning before calling the method.
        * ``"pass"``: Just call the method.

        See :attr:`AsArray.delete_atom`.

    \**kwargs
        Keyword arguments for the :func:`numpy.array` function.
        See :attr:`AsArray.kwargs`.

    """

    _MOl: Molecule = Molecule()

    @classmethod
    def from_array(cls) -> Callable[..., None]: return cls._MOl.from_array

    def __init__(self, mol: Iterable[Atom], delete_atom: str = 'raise', **kwargs: Any) -> None:
        """Initialize a :class:`AsArray` instance."""
        self.mol = mol if isinstance(mol, (Sequence, Molecule)) else tuple(mol)
        self.kwargs = kwargs
        self.delete_atom = delete_atom.lower()

        self._xyz = None

        if self.delete_atom == 'pass':
            return None

        try:
            func_new = getattr(self, f'_{self.delete_atom}')
        except AttributeError as ex:
            raise ValueError("An invalid value was passed to the 'delete_atom' parameter: "
                             f"{reprlib.repr(delete_atom)}; accepted values are 'raise', 'warn' "
                             "and 'pass'").with_traceback(ex.__traceback__)
        method = types.MethodType(func_new, mol)
        setattr(mol, 'delete_atom', method)
        return None

    def __enter__(self) -> np.ndarray:
        """Enter the context manager; return an array of Cartesian coordinates."""
        self._xyz = np.array(self.mol, **self.kwargs)
        return self._xyz

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; update the Cartesian coordinates of :attr:`AsArray.mol`."""
        mol = self.mol
        xyz = self._xyz
        self.from_array(xyz, atom_subset=mol)

        # Deleting removes the 'delete_atom' method from the molecules' instance variables;
        # thus reverting back to the method originally defined in the class itself
        if self.delete_atom in ('raise', 'warn'):
            delattr(mol, 'delete_atom')
        self.__dict__ = {}

    """######################### Warning- and exception-related methods #########################"""

    @staticmethod
    def _raise(*args, **kwargs) -> NoReturn:
        raise MoleculeError("Atoms should not be deleted while a molecule is opened with 'AsArray'")

    @staticmethod
    def _warn(*args, **kwargs) -> None:
        warnings.warn("Atoms should not be deleted while a molecule is opened with 'AsArray'")
        return Molecule.delete_atom(*args, **kwargs)

    @staticmethod
    def _pass(*args, **kwargs) -> None:
        return Molecule.delete_atom(*args, **kwargs)

    _raise.__doc__ = _warn.__doc__ = _pass.__doc__ = Molecule.delete_atom.__doc__
