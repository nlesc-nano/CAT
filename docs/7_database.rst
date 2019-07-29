.. _Database:

The Database Class
==================

A Class designed for the storing, retrieval and updating of results.


.. image:: _images/Database.png
    :scale: 35 %
    :align: center

The methods of the Database class can be divided into three categories
accoring to their functionality:

-   Opening & closing the database - these methods serve as context managers
    for loading and unloading parts of the database from the harddrive.

    The context managers can be accessed via the :meth:`.MetaManager.open`
    method of :attr:`.Database.csv_lig`, :attr:`.Database.csv_qd`,
    :attr:`.Database.yaml` or :attr:`.Database.hdf5`, with the option
    of passing additional positional or keyword arguments.

    .. code:: python

        >>> import CAT

        >>> database = CAT.Database()
        >>> with database.csv_lig.open(write=False) as db:
        >>>     print(repr(db))
        DFCollection(df=<pandas.core.frame.DataFrame at 0x7ff8e958ce80>)

        >>> with database.yaml.open() as db:
        >>>     print(type(db))
        <class 'scm.plams.core.settings.Settings'>

        >>> with database.hdf5.open('r') as db:
        >>>     print(type(db))
        <class 'h5py._hl.files.File'>

-   Importing to the database - these methods handle the importing of new data
    from python objects to the Database class:

    ===================  ====================  ====================  =======================
    :meth:`.update_csv`  :meth:`.update_yaml`  :meth:`.update_hdf5`  :meth:`.update_mongodb`
    ===================  ====================  ====================  =======================

-   Exporting from the database - these methods handle the exporting of data
    from the Database class to other python objects or remote locations:

    =================  ==================
    :meth:`.from_csv`  :meth:`.from_hdf5`
    =================  ==================


Index
-----

.. currentmodule:: dataCAT.database.Database
.. autosummary::

    dirname
    csv_lig
    csv_qd
    hdf5
    yaml
    mongodb

    update_mongodb
    update_csv
    update_yaml
    update_hdf5
    from_csv
    from_hdf5

.. currentmodule:: dataCAT
.. autosummary::

    df_collection.get_df_collection
    database_functions.as_pdb_array
    database_functions.from_pdb_array
    database_functions.sanitize_yaml_settings


Class API
---------

Database
~~~~~~~~

.. autoclass:: dataCAT.database.Database
    :members:


DFCollection
~~~~~~~~~~~~

.. autoclass:: dataCAT.df_collection._DFCollection
    :members:


MetaManager
~~~~~~~~~~~

.. autoclass:: dataCAT.context_managers.MetaManager
    :members:


OpenLig
~~~~~~~

.. autoclass:: dataCAT.context_managers.OpenLig
    :members:


OpenQD
~~~~~~

.. autoclass:: dataCAT.context_managers.OpenQD
    :members:


OpenYaml
~~~~~~~~

.. autoclass:: dataCAT.context_managers.OpenYaml
    :members:


Function API
------------

.. autofunction:: dataCAT.df_collection.get_df_collection

.. autofunction:: dataCAT.database_functions.as_pdb_array

.. autofunction:: dataCAT.database_functions.from_pdb_array

.. autofunction:: dataCAT.database_functions.sanitize_yaml_settings


.. _dataCAT.DFCollection: 7_database.html#dfcollection
.. _dataCAT.MetaManager: 7_database.html#metamanager
.. _rdkit.Chem.Mol: http://rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol
.. _h5py.File: http://docs.h5py.org/en/stable/high/file.html
.. _plams.Settings: https://www.scm.com/doc/plams/components/settings.html
.. _plams.Molecule: https://www.scm.com/doc/plams/components/molecule.html#id1
.. _np.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
.. _np.float64: https://docs.scipy.org/doc/numpy/user/basics.types.html#array-types-and-conversions-between-types
.. _np.int64: https://docs.scipy.org/doc/numpy/user/basics.types.html#array-types-and-conversions-between-types
.. _np.str_: https://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html#arrays-dtypes
.. _np.bytes: https://docs.scipy.org/doc/numpy-1.15.1/reference/arrays.dtypes.html
.. _pd.Series: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.html
.. _pd.DataFrame: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
.. _dict: https://docs.python.org/3/library/stdtypes.html#dict
.. _list: https://docs.python.org/3/library/stdtypes.html#list
.. _tuple: https://docs.python.org/3/library/stdtypes.html#tuple
.. _str: https://docs.python.org/3/library/stdtypes.html#str
.. _int: https://docs.python.org/3/library/functions.html#int
.. _None: https://docs.python.org/3.7/library/constants.html#None
.. _bool: https://docs.python.org/3/library/functions.html?highlight=bool#bool
.. _with: https://docs.python.org/3/reference/compound_stmts.html#with
.. _type: https://docs.python.org/3/library/functions.html#type
.. _Sequence: https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence
.. _AbstractContextManager: https://docs.python.org/3/library/contextlib.html#contextlib.AbstractContextManager

.. |dataCAT.DFCollection| replace:: *dataCAT.DFCollection*
.. |dataCAT.MetaManager| replace:: *dataCAT.MetaManager*
.. |rdkit.Chem.Mol| replace:: *rdkit.Chem.Mol*
.. |h5py.File| replace:: *h5py.File*
.. |plams.Molecule| replace:: *plams.Molecule*
.. |plams.Settings| replace:: *plams.Settings*
.. |pd.Series| replace:: *pd.Series*
.. |pd.DataFrame| replace:: *pd.DataFrame*
.. |np.ndarray| replace:: *np.ndarray*
.. |np.float64| replace:: *np.float64*
.. |np.int64| replace:: *np.int64*
.. |np.str_| replace:: *np.str_*
.. |np.bytes| replace:: *np.bytes*
.. |dict| replace:: *dict*
.. |list| replace:: *list*
.. |tuple| replace:: *tuple*
.. |str| replace:: *str*
.. |int| replace:: *int*
.. |None| replace:: *None*
.. |bool| replace:: *bool*
.. |with| replace:: ``with``
.. |type| replace:: *type*
.. |Sequence| replace:: *Sequence*
.. |AbstractContextManager| replace:: *AbstractContextManager*
