"""
CAT.attachment.remove_atoms_cm
==============================

A context manager for temporary removing a set of atoms from a molecule.

Index
-----
.. currentmodule:: CAT.attachment.remove_atoms_cm
.. autosummary::
    RemoveAtoms

API
---
.. autoclass:: RemoveAtoms
    :members:
    :private-members:
    :special-members:

"""

from typing import Iterable, Union, Sequence
from contextlib import AbstractContextManager
from collections import OrderedDict, abc

from scm.plams import Molecule, Atom, MoleculeError

__all__ = ['RemoveAtoms']


class RemoveAtoms(AbstractContextManager):
    """A context manager for temporary removing a set of atoms from a molecule.

    The *relative* ordering of the to-be removed atoms (and matching bonds),
    as specified in **atoms**, is preserved during the removal and reattachment process.
    Note that reattaching will (re-)append the removed atoms/bonds,
    a process which is thus likelly to affect the *absolute* ordering of
    atoms/bonds within the entire molecule.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule, Atom

        >>> mol: Molecule = ...  # A random molecule
        >>> atom1: Atom = mol[1]
        >>> atom2: Atom = mol[2]

        >>> atom_set = {atom1, atom2}
        >>> with RemoveAtoms(mol, atom_set):
        ...     print(atom1 in mol, atom2 in mol)
        False, False

        >>> print(atom1 in mol, atom2 in mol)
        True, True

    Parameters
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.
        See :attr:`RemoveAtoms.mol`.

    atoms : |plams.Atom| or |Iterable| [|plams.Atom|]
        A PLAMS atom or an iterable consisting of unique PLAMS atoms.
        All supplied atoms should belong to **mol**.
        See :attr:`RemoveAtoms.atoms`.

    Attributes
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.

    atoms : |Sequence| [|plams.Atom|]
        A sequence of PLAMS atoms belonging to :attr:`RemoveAtoms.mol`.
        Setting a value will convert it into a sequence of atoms.

    bonds : |OrderedDict| [|plams.Bond|, ``None``], optional
        A ordered dictionary of PLAMS bonds connected to one or more atoms in
        :attr:`RemoveAtoms.atoms`.
        All values are ``None``, the dictionary serving as an improvised ``OrderedSet``.
        Set to ``None`` until :meth:`RemoveAtoms.__enter__` is called.

    """

    @property
    def atoms(self) -> Sequence[Atom]:
        """Get or set :attr:`RemoveAtoms.atoms`. Setting will convert the supplied value into a sequence of atoms."""  # noqa
        return self._atoms

    @atoms.setter
    def atoms(self, value: Union[Atom, Iterable[Atom]]) -> None:
        if isinstance(value, Atom):
            self._atoms = (value,)
        elif not isinstance(value, abc.Sequence):
            self._atoms = tuple(value)
        else:
            self._atoms = value

    def __init__(self, mol: Molecule, atoms: Union[Atom, Iterable[Atom]]) -> None:
        """Initialize a :class:`RemoveAtoms` instance."""
        self.mol = mol
        self.atoms = atoms
        self.bonds = None

    def __enter__(self) -> None:
        """Enter the context manager; delete all atoms in :class:`RemoveAtoms.atoms`."""
        mol = self.mol
        self.bonds = bonds_set = OrderedDict()  # An improvised "OrderedSet"
        for atom in self.atoms:
            for bond in atom.bonds:
                bonds_set[bond] = None
            mol.delete_atom(atom)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; reassign all atoms in :class:`RemoveAtoms.atoms`."""
        mol = self.mol
        for atom in reversed(self.atoms):
            mol.add_atom(atom)
        for bond in reversed(self.bonds):
            try:
                mol.add_bond(bond)
            except MoleculeError:
                pass  # One of the bonded atoms has been manually deleted by the user: skip it
        self.bonds = None
