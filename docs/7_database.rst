.. _Database:

The Database Class
==================
.. currentmodule:: dataCAT

A Class designed for the storing, retrieval and updating of results.

.. image:: _images/Database.png
    :scale: 35 %
    :align: center

The methods of the Database class can be divided into three categories
accoring to their functionality:

-   Opening & closing the database - these methods serve as context managers
    for loading and unloading parts of the database from the harddrive.

    The context managers can be accessed by calling
    either :attr:`.Database.csv_lig`, :attr:`.Database.csv_qd`,
    :attr:`.Database.yaml` or :attr:`.Database.hdf5`, with the option
    of passing additional positional or keyword arguments.

    .. code:: python

        >>> from dataCAT import Database

        >>> database = Database()
        >>> with database.csv_lig(write=False) as db:
        >>>     print(repr(db))
        DFProxy(ndframe=<pandas.core.frame.DataFrame at 0x7ff8e958ce80>)

        >>> with database.yaml() as db:
        >>>     print(type(db))
        <class 'scm.plams.core.settings.Settings'>

        >>> with database.hdf5('r') as db:
        >>>     print(type(db))
        <class 'h5py._hl.files.File'>

-   Importing to the database - these methods handle the importing of new data
    from python objects to the Database class:

    ============================  =============================  =============================  ================================
    :meth:`~Database.update_csv`  :meth:`~Database.update_yaml`  :meth:`~Database.update_hdf5`  :meth:`~Database.update_mongodb`
    ============================  =============================  =============================  ================================

-   Exporting from the database - these methods handle the exporting of data
    from the Database class to other python objects or remote locations:

    ==========================  ===========================
    :meth:`~Database.from_csv`  :meth:`~Database.from_hdf5`
    ==========================  ===========================


Index
-----

.. autosummary::

    Database.dirname
    Database.csv_lig
    Database.csv_qd
    Database.hdf5
    Database.yaml
    Database.mongodb

    Database.update_mongodb
    Database.update_csv
    Database.update_yaml
    Database.update_hdf5
    Database.from_csv
    Database.from_hdf5

    DFProxy
    OpenLig
    OpenQD
    OpenYaml

    PDBContainer
    DTYPE_ATOM
    DTYPE_BOND


API
---
.. autoclass:: Database
    :members:

.. autoclass:: DFProxy

.. autoclass:: OpenLig
    :members:
    :inherited-members:

.. autoclass:: OpenQD
    :members:
    :inherited-members:

.. autoclass:: OpenYaml
    :members:
    :inherited-members:

.. autoclass:: PDBContainer
    :members:

.. data:: DTYPE_ATOM
    :type: Mapping[str, numpy.dtype]
    :value: ...

    A mapping representing the dtype of :attr:`PDBContainer.atoms`.

    Most field names are based on to their, identically named, counterpart as produced by
    :func:`readpdb()<scm.plams.interfaces.molecule.rdkit.readpdb>`,
    the data in question being stored in the
    :class:`Atom.properties.pdb_info<scm.plams.mol.atom.Atom>` block.

    There are four exception to this general rule:

    * ``x``, ``y`` & ``z``: Based on :class:`Atom.coords<scm.plams.mol.atom.Atom>`.
    * ``symbol``: Based on :class:`Atom.symbol<scm.plams.mol.atom.Atom>`.
    * ``charge``: Based on :class:`Atom.properties.charge<scm.plams.mol.atom.Atom>`.
    * ``charge_float``: Based on :class:`Atom.properties.charge_float<scm.plams.mol.atom.Atom>`.

    .. code:: python

        >>> from dataCAT import DTYPE_ATOM
        >>> print(DTYPE_ATOM)
        mappingproxy({
            'IsHeteroAtom': dtype('bool'),
            'SerialNumber': dtype('int16'),
            'Name': dtype('S4'),
            'ResidueName': dtype('S3'),
            'ChainId': dtype('S1'),
            'ResidueNumber': dtype('int16'),
            'x': dtype('float32'),
            'y': dtype('float32'),
            'z': dtype('float32'),
            'Occupancy': dtype('float32'),
            'TempFactor': dtype('float32'),
            'symbol': dtype('S4'),
            'charge': dtype('int8'),
            'charge_float': dtype('float64')
        })


.. data:: DTYPE_BOND
    :type: Mapping[str, numpy.dtype]
    :value: ...

    A mapping representing the dtype of :attr:`PDBContainer.bonds`.

    Field names are based on to their, identically named,
    counterpart in :class:`~scm.plams.mol.bond.Bond`.

    .. code:: python

        >>> from dataCAT import DTYPE_BOND
        >>> print(DTYPE_BOND)
        mappingproxy({
            'atom1': dtype('int32'),
            'atom2': dtype('int32'),
            'order': dtype('int8')
        })
