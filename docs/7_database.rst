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

    PDBTuple
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

.. autoclass:: PDBTuple
    :members:

.. data:: DTYPE_ATOM
    :value: pandas.DataFrame(...)

    A :class:`pandas.DataFrame` representing the dtype of :attr:`PDBTuple.atoms`.

    Most field names are based on to their, identically named, counterpart as produced by
    :func:`readpdb()<scm.plams.interfaces.molecule.rdkit.readpdb>`,
    the data in question being stored in the
    :class:`Atom.properties.pdb_info<scm.plams.mol.atom.Atom>` block.

    There are four exception to this general rule:

    * ``coords``: Based on :class:`Atom.coords<scm.plams.mol.atom.Atom>`.
    * ``symbol``: Based on :class:`Atom.symbol<scm.plams.mol.atom.Atom>`.
    * ``charge``: Based on :class:`Atom.properties.charge<scm.plams.mol.atom.Atom>`.
    * ``charge_float``: Based on :class:`Atom.properties.charge_float<scm.plams.mol.atom.Atom>`.

    .. code:: python

        >>> from dataCAT import DTYPE_ATOM
        >>> print(DTYPE_ATOM)  # doctest: +SKIP
                         dtype  shape
        name
        IsHeteroAtom      bool      0
        SerialNumber     int16      0
        Name               |S4      0
        ResidueName        |S3      0
        ChainId            |S1      0
        ResidueNumber    int16      0
        coords         float32      3
        Occupancy      float32      0
        TempFactor     float32      0
        symbol             |S4      0
        charge            int8      0
        charge_float   float64      0


.. data:: DTYPE_BOND
    :value: pandas.DataFrame(...)

    A :class:`pandas.DataFrame` representing the dtype of :attr:`PDBTuple.bonds`.

    .. code:: python

        >>> from dataCAT import DTYPE_BOND
        >>> print(DTYPE_BOND)
               dtype  shape
        name
        atoms  int32      2
        order   int8      0
