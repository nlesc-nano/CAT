"""
CAT.attachment.mol_split_cm
===========================

A context manager for temporary splitting (and reassembling) molecules.

Index
-----
.. currentmodule:: CAT.attachment.mol_split_cm
.. autosummary::
    SplitMol

API
---
.. autoclass:: SplitMol
    :members:
    :private-members:
    :special-members:

"""

import reprlib
from typing import Iterable, Union, Dict, Tuple
from contextlib import AbstractContextManager

from scm.plams import Molecule, Atom, PT, Bond, MoleculeError, PTError, rotation_matrix

from ..mol_utils import separate_mod


class SplitMol(AbstractContextManager):
    """A context manager for temporary splitting a single molecule into multiple components.

    The context manager splits the provided molecule into multiple components,
    capping all broken bonds in the process.
    The exact amount of fragments depends on the number of specified bonds.

    These moleculair fragments are returned upon opening the context manager and merged back
    into the initial molecule once the context manager is closed. While opened, the initial molecule
    is cleared of all atoms and bonds, while the same hapens to the moleculair fragments upon
    closing.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule, Bond

        >>> mol: Molecule = ...  # Placeholder for a random molecule
        >>> bond: Bond = mol[1, 2]

        # A backup of all bonds and atoms
        >>> bonds_backup = mol.bonds.copy()
        >>> atoms_backup = mol.atoms.copy()

        >>> with SplitMol(mol, bond) as mol_list:
        ...     mol1, mol2 = mol_list
        ...     fancy_operation1(mol1)
        ...     fancy_operation2(mol2)
        ...     print(
        ...         mol.bonds == bonds_backup,
        ...         mol.atoms == atoms_backup,
        ...         bond in mol.bonds
        ...     )
        False, False, False

        # The context manager is closed; all atoms and bonds have been restored
        >>> print(
        ...     mol.bonds == bonds_backup,
        ...     mol.atoms == atoms_backup,
        ...     bond in mol.bonds
        ... )
        True, True, True

    Parameters
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.
        See :attr:`SplitMol.mol`.

    bond_list : |plams.Bond| or |Iterable| [|plams.Bond|]
        An iterable consisting of PLAMS bonds.
        All bonds must be part of **mol**.
        See :attr:`SplitMol.bonds`.

    cap_type : :class:`str` or :class:`int`
        An atomic number or symbol of the atom type used for capping the to-be split molecule.
        See :attr:`SplitMol.cap_type`.

    Attributes
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.

    bonds : :class:`dict` [|plams.Bond|, :class:`int`]
        An dictionary with PLAMS bonds as key and their (0-based) index as value.
        All bonds must be part of **mol**.

    cap_type : :class:`str`
        An atomic symbol of the atom type used for capping the to-be split molecule.

    _at_pairs : :class:`list` [:class:`dict` [|plams.Atom|, |plams.Atom|]], optional
        A list of dictionaries.
        Each dictionary contains an atom (see :attr:`SplitMol.bond_list`) and its respective
        capping atom as value.
        Used for reassembling :attr:`SplitMol.mol` once the context manager is closed.
        Set internally by :meth:`SplitMol.__enter__`.

    _vars_backup : :class:`dict` [:class:`str`, |Any|], optional
        A backup of all instance variables of :attr:`SplitMol.mol`.
        Set internally by :meth:`SplitMol.__enter__`.

    _tmp_mol_list : :class:`tuple` [|plams.Molecule|], optional
        A list of PLAMS molecules obtained by splitting :attr:`SplitMol.mol`.
        Set internally by :meth:`SplitMol.__enter__`.

    """

    def __init__(self, mol: Molecule,
                 bond_list: Union[Bond, Iterable[Bond]],
                 cap_type: Union[str, int, Atom] = 'H') -> None:
        """Initialize a :class:`SplitMol` instance."""
        self.mol = mol
        self.bonds = bond_list
        self.cap_type = cap_type

        # Private variables which are modified in __enter__ and __exit__
        self._at_pairs = None
        self._vars_backup = None
        self._tmp_mol_list = None

    """####################################### Properties #######################################"""

    @property
    def bonds(self) -> Dict[Bond, int]:
        """Getter: Return :attr:`SplitMol.bonds`. Setter: Convert the provided value into a dictionary with bonds as keys and their matching index as value."""  # noqa
        return self._bonds

    @bonds.setter
    def bonds(self, value: Union[Bond, Iterable[Bond]]) -> None:
        bonds = (value,) if isinstance(value, Bond) else value
        try:
            self._bonds = {bond: self.mol.bonds.index(bond) for bond in bonds}
        except ValueError:
            raise MoleculeError('Passed bonds should belong to mol')

    @property
    def cap_type(self) -> str:
        """Getter: Return :attr:`SplitMol.cap_type`. Setter: Convert the provided value into a valid atom type (:class:`str`)."""  # noqa
        return self._cap_type

    @cap_type.setter
    def cap_type(self, value: Union[str, int, Atom]) -> None:
        if isinstance(value, int):
            self._cap_type = PT.get_symbol(value)
        elif isinstance(value, Atom):
            self._cap_type = value.symbol
        elif isinstance(value, str):
            if value.capitalize() not in PT.symtonum:
                raise PTError(f"Trying to convert an invalid atomic symbol: {reprlib.repr(value)}")
            self._cap_type = value
        else:
            raise TypeError(f"Invalid atom type: {repr(type(value))}")

    """Methods for manipulating the molecule."""

    def split_bond(self, bond: Bond) -> Dict[Atom, Atom]:
        """Delete a bond from :attr:`SplitMol.mol` and cap the resulting fragments with :attr:`SplitMol.cap_type`.

        Parameters
        ----------
        bond : |plams.Bond|
            A PLAMS bond.

        Returns
        -------
        :class:`dict` [|plams.Atom|, |plams.Atom|]
            A dictionary with the old atoms in **bond** as keys and their new capping atoms
            as values.

        """  # noqa
        mol = self.mol
        symbol = self.cap_type

        # Construct the capping atoms
        at_1a, at_2a = bond.atom1, bond.atom2
        at_1b = Atom(symbol=symbol, coords=at_2a.coords)
        at_2b = Atom(symbol=symbol, coords=at_1a.coords)
        mol.add_atom(at_1b, adjacent=[at_1a])
        mol.add_atom(at_2b, adjacent=[at_2a])

        # Resize the capping atom bond
        length1 = at_1a.radius + at_1b.radius
        length2 = at_2a.radius + at_2b.radius
        mol.bonds[-2].resize(at_1b, length1)
        mol.bonds[-1].resize(at_2b, length2)

        # Delete the old bond and return a dictionary containg marking all new bonds
        mol.delete_bond(bond)
        return {at_1a: at_1b, at_2a: at_2b}

    def reassemble(self) -> None:
        """Reassemble :attr:`SplitMol.mol` from its consitituent components in :attr:`SplitMol._tmp_mol_list`.

        All capping atoms (and their respective bonds) specified in :attr:`SplitMol._at_pairs` are
        removed and the previously broken bonds (stored in :attr:`SplitMol.bonds`) are restored.

        """  # noqa
        mol = self.mol
        mark = self._at_pairs
        bond_iterator = sorted(self.bonds.items(), key=lambda x: x[-1])

        for atom_dict, (bond, idx) in zip(mark, bond_iterator):
            # Extract atoms
            iterator = iter(atom_dict.items())
            at_1a, at_1b = next(iterator)
            at_2a, at_2b = next(iterator)

            # Allign the molecules
            vec1 = at_2a.vector_to(at_2b)
            vec2 = at_1a.vector_to(at_1b)
            rotmat = rotation_matrix(vec1, vec2)
            at_2a.mol.rotate(rotmat)

            # Delete the capping atoms and bonds
            self._uncap_mol(at_1b, at_2b)

            # Add a new bond
            bond.atom1.mol = bond.atom2.mol = mol
            mol.add_bond(bond)
            mol.bonds.insert(idx, mol.bonds.pop())

            # Resize the bond
            length = at_1a.radius + at_2a.radius
            bond.resize(at_1a, length)

        # Ensure all atoms and bonds belong to mol
        for at in mol.atoms:
            at.mol = mol
        for bond in mol.bonds:
            bond.mol = mol

    def _uncap_mol(self, *cap_atoms: Atom) -> None:
        """Remove the specified capping atoms (and their respective bonds) from :attr:`SplitMol.mol`.

        Parameters
        ----------
        \*cap_atoms : |plams.Atom|
            The capping atoms which are to-be removed.

        """  # noqa
        mol = self.mol
        for atom in cap_atoms:
            bond = atom.bonds[0]
            other_end = bond.other_end(atom)

            mol.atoms.remove(atom)
            mol.bonds.remove(bond)
            other_end.bonds.remove(bond)

    def reset_vars(self) -> None:
        """Reset :attr:`Molecule.__dict__` and assign a backup to :attr:`SplitMol._vars_backup`."""
        mol = self.mol

        self._vars_backup = vars(mol)
        mol.__dict__ = vars(mol).copy()
        mol.bonds = []
        mol.atoms = []

    """############################# Context manager magic methods #############################"""

    def __enter__(self) -> Tuple[Molecule]:
        """Enter the :class:`SplitMol` context manager; return a list of molecules."""
        # Break bonds and add capping atoms to mol
        # Split mol into multiple seperate molecules without copying atoms
        mol = self.mol
        self._at_pairs = [self.split_bond(bond) for bond in self.bonds]
        self._tmp_mol_list = mol.separate_mod()

        # Pop all instance variables of mol and return a list of temporary molecules
        self.reset_vars()
        return self._tmp_mol_list

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the :class:`SplitMol` context manager."""
        # Restore all instance variables of mol
        mol = self.mol
        vars(mol).update(self._vars_backup)

        # Delete all capping atoms recreate the previously broken bond
        self.reassemble()

        # Delete all instance variables of the temporary molecules
        for mol_tmp in self._tmp_mol_list:
            mol_tmp.__dict__.update({'atoms': [], 'bonds': []})

        # Reset all private instance variables of the context manager
        self._mark = self._vars_backup = self._tmp_mol_list = None
