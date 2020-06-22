.. _PDBContainer:

The PDBContainer Class
======================
.. currentmodule:: dataCAT

An immutable class for holding array-like representions of a set of .pdb files.

Index
-----

.. autosummary::

    PDBContainer.atoms
    PDBContainer.bonds
    PDBContainer.atom_count
    PDBContainer.bond_count

    PDBContainer.__init__
    PDBContainer.__getitem__
    PDBContainer.__len__
    PDBContainer.keys
    PDBContainer.values
    PDBContainer.items

    PDBContainer.from_molecules
    PDBContainer.to_molecules
    PDBContainer.create_hdf5_group
    PDBContainer.validate_hdf5
    PDBContainer.from_hdf5
    PDBContainer.to_hdf5

    DTYPE_ATOM
    DTYPE_BOND


API
---
.. autoclass:: PDBContainer
    :members: atoms, bonds, atom_count, bond_count, __init__

.. automethod:: PDBContainer.__getitem__
.. automethod:: PDBContainer.__len__
.. automethod:: PDBContainer.keys
.. automethod:: PDBContainer.values
.. automethod:: PDBContainer.items
.. automethod:: PDBContainer.from_molecules
.. automethod:: PDBContainer.to_molecules
.. automethod:: PDBContainer.create_hdf5_group
.. automethod:: PDBContainer.validate_hdf5
.. automethod:: PDBContainer.from_hdf5
.. automethod:: PDBContainer.to_hdf5

.. data:: DTYPE_ATOM
    :type: Mapping[str, np.dtype]
    :value: ...

    A mapping representing the dtype of :attr:`PDBContainer.atoms`.

    Most field names are based on to their, identically named, counterpart as produced by
    :func:`readpdb()<scm.plams.interfaces.molecule.rdkit.readpdb>`,
    the data in question being stored in the
    :class:`Atom.properties.pdb_info<scm.plams.mol.atom.Atom>` block.

    There are six exception to this general rule:

    * ``x``, ``y`` & ``z``: Based on :class:`Atom.x<scm.plams.mol.atom.Atom>`,
      :class:`Atom.y<scm.plams.mol.atom.Atom>` and :class:`Atom.z<scm.plams.mol.atom.Atom>`.
    * ``symbol``: Based on :class:`Atom.symbol<scm.plams.mol.atom.Atom>`.
    * ``charge``: Based on :class:`Atom.properties.charge<scm.plams.mol.atom.Atom>`.
    * ``charge_float``: Based on :class:`Atom.properties.charge_float<scm.plams.mol.atom.Atom>`.

    .. code:: python

        mappingproxy({
            'IsHeteroAtom':  dtype('bool'),
            'SerialNumber':  dtype('int16'),
            'Name':          dtype('S4'),
            'ResidueName':   dtype('S3'),
            'ChainId':       dtype('S1'),
            'ResidueNumber': dtype('int16'),
            'x':             dtype('float32'),
            'y':             dtype('float32'),
            'z':             dtype('float32'),
            'Occupancy':     dtype('float32'),
            'TempFactor':    dtype('float32'),
            'symbol':        dtype('S4'),
            'charge':        dtype('int8'),
            'charge_float':  dtype('float64')
        })


.. data:: DTYPE_BOND
    :type: Mapping[str, np.dtype]
    :value: ...

    A mapping representing the dtype of :attr:`PDBContainer.bonds`.

    Field names are based on to their, identically named,
    counterpart in :class:`~scm.plams.mol.bond.Bond`.

    .. code:: python

        mappingproxy({
            'atom1': dtype('int32'),
            'atom2': dtype('int32'),
            'order': dtype('int8')
        })
