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

from typing import Iterable
from contextlib import AbstractContextManager

from scm.plams import Molecule, Atom, MoleculeError


class RemoveAtoms(AbstractContextManager):
    """A context manager for temporary removing a set of atoms from a molecule.

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

    atoms : |Iterable| [|plams.Atom|]
        An iterable with PLAMS atoms beloning to **mol**.
        See :attr:`RemoveAtoms.atoms`.

    Attributes
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.

    atoms : :class:`tuple` [|plams.Atom|]
        A tuple of PLAMS atoms belonging to :attr:`RemoveAtoms.mol`.

    bonds : :class:`set` [|plams.Bond|], optional
        A set of PLAMS bonds connected to one or more atoms in :attr:`RemoveAtoms.atoms`.
        Set to ``None`` until :meth:`RemoveAtoms.__enter__` is called.

    """

    def __init__(self, mol: Molecule, atoms: Iterable[Atom]) -> None:
        """Initialize a :class:`RemoveAtoms` instance."""
        self.mol = mol
        self.atoms = tuple(atoms)  # This should be a collection to preserve atom ordering
        self.bonds = None

    def __enter__(self) -> None:
        """Enter the :class:`RemoveAtoms` context manager."""
        mol = self.mol
        self.bonds = bonds_set = set()
        for atom in reversed(self.atoms):
            for bond in atom.bonds:
                bonds_set.add(bond)
            mol.delete_atom(atom)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the :class:`RemoveAtoms` context manager."""
        mol = self.mol
        for atom in self.atoms:
            mol.add_atom(atom)
        for bond in self.bonds:
            try:
                mol.add_bond(bond)
            except MoleculeError:
                pass  # One of the bonded atoms has been manually deleted by the user: skip it
        self.bonds = None
