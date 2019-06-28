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
    These methods should be used in conjunction with |with|_ statements:

    ::

        import CAT

        database = CAT.Database()
        with database.open_csv_lig(db.csv_lig) as db:
            print('my ligand database')
        with database.open_yaml(db.yaml) as db:
            print('my job settings database')
        with h5py.File(db.hdf5) as db:
            print('my structure database')

    ======================  =====================  ===================  ==================
    :class:`.open_csv_lig`  :class:`.open_csv_qd`  :class:`.open_yaml`  :class:`h5py.File`
    ======================  =====================  ===================  ==================

-   Importing to the database - these methods handle the importing of new data
    from python objects to the Database class:

    ===================  ====================  ====================
    :meth:`.update_csv`  :meth:`.update_yaml`  :meth:`.update_hdf5`
    ===================  ====================  ====================

-   Exporting from the database - these methods handle the exporting of data
    from the Database class to other python objects or remote locations:

    =================  ==================
    :meth:`.from_csv`  :meth:`.from_hdf5`
    =================  ==================

Index
~~~~~

.. currentmodule:: data_CAT.data_handling.database.Database
.. autosummary::

    open_yaml
    open_csv_lig
    open_csv_qd
    DF
    update_csv
    update_yaml
    update_hdf5
    from_csv
    from_hdf5


.. currentmodule:: data_CAT.data_handling.database_functions
.. autosummary::

    mol_to_file
    as_pdb_array
    from_pdb_array
    sanitize_yaml_settings


Class API
~~~~~~~~~

.. autoclass:: data_CAT.data_handling.database.Database
    :members:

Function API
~~~~~~~~~~~~

.. autofunction:: data_CAT.data_handling.database_functions.mol_to_file

.. autofunction:: data_CAT.data_handling.database_functions.as_pdb_array

.. autofunction:: data_CAT.data_handling.database_functions.from_pdb_array

.. autofunction:: data_CAT.data_handling.database_functions.sanitize_yaml_settings

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
.. _with: https://docs.python.org/3/reference/compound_stmts.html#with

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
.. |with| replace:: ``with``
