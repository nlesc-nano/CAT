""" A module which holds the Database class. """

__all__ = ['Database', 'mol_to_file']

from os import getcwd
from os.path import (join, isfile, isdir)

import yaml
import h5py
import numpy as np
import pandas as pd

from scm.plams import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..mol_utils import from_rdmol
from ..utils import (get_time, get_template)


def mol_to_file(mol_list, path=None, overwrite=False, mol_format=['xyz', 'pdb']):
    """ Export all molecules in **mol_list** to .pdb and/or .xyz files.

    :parameter mol_list: A list of PLAMS molecules.
    :type mol_list: |list|_ [|plams.Molecule|_]
    :parameter path: The path to the directory where the molecules will be stored. Defaults
        to the current working directory if *None*.
    :type path: |None|_ or |str|_
    :parameter bool overwrite: If previously generated structures can be overwritten or not.
    :parameter mol_format: A list of strings with the to-be exported file types. Accepted values
        are *xyz* and/or *pdb*.
    :type mol_format: |list|_ [|str|_]
    """
    # Set the export path
    path = path or getcwd()
    assert isdir(path)

    if overwrite:  # Export molecules while allowing for file overriding
        for mol in mol_list:
            mol_path = join(path, mol.properties.name)
            if 'pdb' in mol_format:
                molkit.writepdb(mol, mol_path + '.pdb')
            if 'xyz' in mol_format:
                mol.write(mol_path + '.xyz')

    else:  # Export molecules without allowing for file overriding
        for mol in mol_list:
            mol_path = join(path, mol.properties.name)
            if 'pdb' in mol_format and not isfile(mol_path + '.pdb'):
                molkit.writepdb(mol, mol_path + '.pdb')
            if 'xyz' in mol_format and not isfile(path + '.xyz'):
                mol.write(path + '.xyz')


def get_nan_row(df):
    """ Return a list of None-esque objects for each column in **df**.
    The object in question depends on the data type of the column.
    Will default to *None* if a specific data type is not recognized

        * |np.int64|_: *-1*

        * |np.float64|_: *np.nan*

        * |object|_: *None*

    :parameter df: A dataframe
    :type df: |pd.DataFrame|_
    :return: A list of non-esque objects, one for each column in **df**.
    :rtype: |list|_ [|int|_, |float|_ and/or |None|_]
    """
    dtype_dict = {np.dtype('int64'): -1, np.dtype('float64'): np.nan, np.dtype('O'): None}
    if not isinstance(df.index, pd.MultiIndex):
        return [dtype_dict[df[i].dtype] for i in df]
    else:
        ret = []
        for i in df:
            try:
                j = dtype_dict[df[i].dtype]
            except KeyError:  # dtype is neither int, float nor object
                j = None
            ret.append(j)
        return ret


def as_pdb_array(mol_list, min_size=0):
    """ Converts a list of PLAMS molecule into an array of strings representing (partially)
    de-serialized .pdb files.

    :parameter mol_list: A list of PLAMS molecules.
    :type mol_list: |list|_ [|plams.Molecule|_]
    :parameter int min_size: The minimumum length of the pdb_array. The array is padded with empty
        strings if required.
    :return: An array with *m* partially deserialized .pdb files with up to *n* lines each.
    :rtype: *m*n* |np.ndarray|_ [|np.bytes|_ *|S80*]
    """
    pdb_list = []
    shape = min_size
    for mol in mol_list:
        pdb_block = Chem.MolToPDBBlock(molkit.to_rdmol(mol)).splitlines()
        pdb_list.append(pdb_block)
        shape = max(shape, len(pdb_block))

    # Construct, fill and return the pdb array
    shape = len(mol_list), shape
    ret = np.zeros(shape, dtype='S80')
    for i, item in enumerate(pdb_list):
        ret[i][:len(item)] = item

    return ret


def from_pdb_array(array, rdmol=True):
    """ Converts an array with a (partially) de-serialized .pdb file into an
    RDKit or PLAMS molecule.

    :parameter array: A (partially) de-serialized .pdb file with *n* lines.
    :type array: *n* |np.ndarray|_ [|np.bytes|_ / S80]
    :parameter bool rdmol: If *True*, return an RDKit molecule instead of a PLAMS molecule.
    :return: A PLAMS or RDKit molecule build from **array**.
    :rtype: |plams.Molecule|_ or |rdkit.Chem.Mol|_
    """
    pdb_str = ''.join([item.decode() + '\n' for item in array if item])
    ret = Chem.MolFromPDBBlock(pdb_str, removeHs=False, proximityBonding=False)
    if not rdmol:
        return molkit.from_rdmol(ret)
    return ret


def sanitize_yaml_settings(settings, job_type):
    """ Remove a predetermined set of unwanted keys and values from a settings object.

    :param settings: A settings object with, potentially, undesired keys and values.
    :type settings: |plams.Settings|_ (superclass: |dict|_)
    :return: A (nested) dictionary with unwanted keys and values removed.
    :rtype: |dict|_
    """
    def recursive_del(s, s_del):
        for key in s:
            if key in s_del:
                if isinstance(s_del[key], dict):
                    recursive_del(s[key], s_del[key])
                else:
                    del s[key]
            if not s[key]:
                del s[key]

    # Prepare a blacklist of specific keys
    blacklist = get_template('settings_blacklist.yaml')
    settings_del = blacklist['generic']
    settings_del.update(blacklist[job_type])

    # Recursivelly delete all keys from **s** if aforementioned keys are present in the s_del
    recursive_del(settings, settings_del)
    return settings


def _create_csv(path, database='ligand'):
    """ Create a ligand or QD database (csv format) and, if it does not exist, and return
    its absolute path.

    :param str path: The path to the database.
    :param str database: The type of database, accepted values are *ligand* and *qd*.
    :return: The absolute path to the ligand or QD database.
    :rtype: |str|_
    """
    path = join(path, database + '_database.csv')

    # Check if the database exists and has the proper keys; create it if it does not
    if not isfile(path):
        print(get_time() + database + '_database.csv not found in ' +
              path + ', creating ' + database + ' database')
        if database == 'ligand':
            _create_csv_lig(path)
        elif database == 'QD':
            _create_csv_qd(path)
        else:
            raise TypeError(str(database) + " is not an accepated value for the 'database' \
                            argument")
    return path


def _create_csv_lig(path):
    """ Create a ligand database and and return its absolute path.

    :param str path: The path to the database.
    """
    idx = pd.MultiIndex.from_tuples([('-', '-')], names=['smiles', 'anchor'])
    columns_tups = [('hdf5 index', ''), ('formula', ''), ('settings', 1)]
    columns = pd.MultiIndex.from_tuples(columns_tups, names=['index', 'sub index'])
    df = pd.DataFrame(None, index=idx, columns=columns)
    df['hdf5 index'] = -1
    df['formula'] = 'str'
    df['settings'] = 'str'
    df.to_csv(path)


def _create_csv_qd(path):
    """ Create a QD database and and return its absolute path.

    :param str path: The path to the database.
    """
    idx = pd.MultiIndex.from_tuples(
            [('-', '-', '-', '-')],
            names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    )
    columns_tups = [('hdf5 index', ''), ('ligand count', ''), ('settings', 1), ('settings', 2)]
    columns = pd.MultiIndex.from_tuples(columns_tups, names=['index', 'sub index'])
    df = pd.DataFrame(None, index=idx, columns=columns)
    df['hdf5 index'] = -1
    df['ligand count'] = -1
    df['settings'] = 'str'
    df.to_csv(path)


def _create_hdf5(path, name='structures.hdf5'):
    """ Create a pdb structure database (hdf5 format), populate it with the *core*, *ligand*
    and *QD* datasets and finally return its absolute path.

    :param str path: The path to the database.
    :param str name: The filename of the database (excluding its path)
    :return: The absolute path to the pdb structure database.
    :rtype: |str|_
    """
    # Define arguments
    path = join(path, name)
    dataset_names = 'core', 'ligand', 'QD'
    kwarg = {'chunks': True, 'maxshape': (None, None), 'compression': 'gzip'}

    # Create new datasets
    with h5py.File(path, 'a') as f:
        for name in dataset_names:
            if name not in f:
                f.create_dataset(name=name, data=np.empty((0, 1), dtype='S80'), **kwarg)

    return path


def _create_yaml(path, name='job_settings.yaml'):
    """ Create a job settings database (.yaml

    :param str path: The path to the database.
    :param str name: The filename of the database (excluding its path)
    :return: The absolute path to the pdb structure database.
    :rtype: |str|_
    """
    # Define arguments
    path = join(path, name)

    # Create a new .yaml file if it does not yet exist
    if not isfile(path):
        with open(path, 'w') as f:
            f.write(yaml.dump({None: [None]}, default_flow_style=False, indent=4))
    return path


class Database():
    """ The Database class.

    :Atributes:     * **csv_lig** (|str|_) – Path and filename of the .csv file containing all \
                    ligand related results.

                    * **csv_qd** (|str|_) – Path and filename of the .csv file containing all \
                    quantum dot related results.

                    * **yaml** (|str|_) – Path and filename of the .yaml file containing all \
                    job settings.

                    * **hdf5** (|str|_) – Path and filename of the .hdf5 file containing all \
                    structures (as partiallize de-serialized .pdb files).

                    * **mongodb** (|None|_) – *None*.
    """

    def __init__(self, path=None):
        path = path or getcwd()

        # Attributes which hold the absolute paths to various components of the database
        self.csv_lig = _create_csv(path, database='ligand')
        self.csv_qd = _create_csv(path, database='QD')
        self.yaml = _create_yaml(path)
        self.hdf5 = _create_hdf5(path)
        self.mongodb = None  # Placeholder

    def __str__(self):
        ret = Settings()
        attr_dict = vars(self)
        for key in attr_dict:
            ret[key] = type(attr_dict[key])
        return str(ret)

    """ ###########################  Opening and closing the database ######################### """

    class open_yaml():
        """ Context manager for opening and closing the job settings database.

        :param str path: The path+filename to the database component.
        :param bool write: Whether or not the database file should be updated after
            closing **self**.
        """
        def __init__(self, path=None, write=True):
            self.path = path or getcwd()
            self.write = write
            self.settings = None

        def __enter__(self):
            with open(self.path, 'r') as f:
                self.settings = Settings(yaml.load(f, Loader=yaml.FullLoader))
                return self.settings

        def __exit__(self, type, value, traceback):
            if self.write:
                yml_dict = self.settings.as_dict()

                # A fix for Settings.as_dict() not functioning when containg a lists of Settings
                for key in yml_dict:
                    for i, value in enumerate(yml_dict[key]):
                        if isinstance(value, Settings):
                            yml_dict[key][i] = value.as_dict()

                # Write to the .yaml file
                with open(self.path, 'w') as f:
                    f.write(yaml.dump(yml_dict, default_flow_style=False, indent=4))
            self.settings = False

    class open_csv_lig():
        """ Context manager for opening and closing the ligand database.

        :param str path: The path+filename to the database component.
        :param bool write: Whether or not the database file should be updated after
            closing **self**.
        """
        def __init__(self, path=None, write=True):
            self.path = path or getcwd()
            self.write = write
            self.df = None

        def __enter__(self):
            # Open the .csv file
            dtype = {'hdf5 index': int, 'formula': str, 'settings': str}
            self.df = pd.read_csv(self.path, index_col=[0, 1], header=[0, 1], dtype=dtype)

            # Fix the columns
            idx_tups = [(i, '') if 'Unnamed' in j else (i, j) for i, j in self.df.columns]
            columns = pd.MultiIndex.from_tuples(idx_tups, names=self.df.columns.names)
            self.df.columns = columns
            return self.df

        def __exit__(self, type, value, traceback):
            if self.write:
                self.df.to_csv(self.path)
            self.df = None

    class open_csv_qd():
        """ Context manager for opening and closing the quantum dot database.

        :param str path: The path+filename to the database component.
        :param bool write: Whether or not the database file should be updated after
            closing **self**.
        """
        def __init__(self, path=None, write=True):
            self.path = path or getcwd()
            self.write = write
            self.df = None

        def __enter__(self):
            # Open the .csv file
            dtype = {'hdf5 index': int, 'ligand count': np.int64, 'settings': str}
            self.df = pd.read_csv(self.path, index_col=[0, 1, 2, 3], header=[0, 1], dtype=dtype)

            # Fix the columns
            idx_tups = [(i, '') if 'Unnamed' in j else (i, j) for i, j in self.df.columns]
            columns = pd.MultiIndex.from_tuples(idx_tups, names=self.df.columns.names)
            self.df.columns = columns
            return self.df

        def __exit__(self, type, value, traceback):
            if self.write:
                self.df.to_csv(self.path)
            self.df = None

    """ #################################  Updating the database ############################## """

    def update_csv(self, df, database='ligand', columns=None, overwrite=False, job_recipe=None):
        """ Update **self.csv_lig** or **self.csv_qd** with
        (potentially) new user provided settings.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter columns: A list of column keys in **df** which
            (potentially) are to be added to **self**. If *None*: Add all columns.
        :type columns: |None|_ or |list|_ [|tuple|_ [|str|_]]
        :parameter bool overwrite: Whether or not previous entries can be overwritten or not.
        :parameter job_recipe: A Settings object with settings specific to a job.
        :type job_recipe: |None|_ or |plams.Settings|_ (superclass: |dict|_)
        """
        # Operate on either the ligand or quantum dot database
        if database == 'ligand':
            path = self.csv_lig
            open_csv = self.open_csv_lig
        elif database == 'QD':
            path = self.csv_qd
            open_csv = self.open_csv_qd

        # Update **self.yaml**
        if job_recipe is not None:
            job_settings = self.update_yaml(job_recipe)
            for key in job_settings:
                df[('settings', key)] = job_settings[key]

        with open_csv(path, write=True) as db:
            # Update **db.index**
            nan_row = get_nan_row(db)
            for i in df.index:
                if i not in db.index:
                    db.at[i, :] = nan_row
            db['hdf5 index'] = db['hdf5 index'].astype(int, copy=False)  # Fix the data type

            # Filter columns
            if not columns:
                df_columns = df.columns
            else:
                df_columns = columns + [i for i in df.columns if i[0] == 'settings']

            # Update **db.columns**
            for i in df_columns:
                if i not in db.columns:
                    try:
                        db[i] = np.array((None), dtype=df[i].dtype)
                    except TypeError:  # e.g. if csv[i] consists of the datatype np.int64
                        db[i] = -1

            # Update **self.hdf5**; returns a new series of indices
            hdf5_series = self.update_hdf5(df, database=database, overwrite=overwrite)

            # Update **db.values**
            db.update(df[columns], overwrite=overwrite)
            db.update(hdf5_series, overwrite=True)

    def update_yaml(self, job_recipe):
        """ Update **self.yaml** with (potentially) new user provided settings.

        :parameter job_recipe: A settings object with one or more settings specific to a job.
        :type job_recipe: |plams.Settings|_ (superclass: |dict|_)
        :return: A dictionary with the column names as keys and the key for **self.yaml** as
            matching values.
        :rtype: |dict|_ (keys: |str|_, values: |str|_)
        """
        ret = {}
        with self.open_yaml(self.yaml) as db:
            for item in job_recipe:
                # Unpack and sanitize keys
                key = job_recipe[item].key
                if isinstance(key, type):
                    key = str(key).rsplit("'", 1)[0].rsplit('.', 1)[-1]

                # Unpack and sanitize values
                value = job_recipe[item].value
                if isinstance(value, dict):
                    value = sanitize_yaml_settings(value, key)

                # Check if the appropiate key is available in **self.yaml**
                if key not in db:
                    db[key] = []

                # Check if the appropiate value is available in **self.yaml**
                if value in db[key]:
                    ret[item] = key + ' ' + str(db[key].index(value))
                else:
                    db[key].append(value)
                    ret[item] = key + ' ' + str(len(db[key]) - 1)
        return ret

    def update_hdf5(self, df, database='ligand', overwrite=False):
        """ Export molecules (see the *mol* column in **df**) to the structure database.
        Returns a series with the **self.hdf5** indices of all new entries.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool overwrite: Whether or not previous entries can be overwritten or not.
        :return: A series with the index of all new molecules in **self.hdf5**
        :rtype: |pd.Series|_ (index: |str|_, values: |np.int64|_)
        """
        # Identify new and preexisting entries
        new = df['hdf5 index'][df['hdf5 index'] == -1]
        old = df['hdf5 index'][df['hdf5 index'] != -1]

        # Add new entries to the database
        with h5py.File(self.hdf5, 'r+') as f:
            i, j = f[database].shape

            if new.any():
                pdb_array = as_pdb_array(df['mol'][new.index], min_size=j)

                # Reshape and update **self.hdf5**
                k = i + pdb_array.shape[0]
                f[database].shape = k, pdb_array.shape[1]
                f[database][i:k] = pdb_array
                ret = pd.Series(np.arange(i, k), index=new.index, name=('hdf5 index', ''))
            else:
                ret = pd.Series(name=('hdf5 index', ''), dtype=int)

            # If **overwrite** is *True*
            if overwrite and old.any():
                f[database][old] = as_pdb_array(df['mol'][old.index], min_size=j)

        return ret

    """ ########################  Pulling results from the database ########################### """

    def from_csv(self, df, database='ligand', get_mol=True, inplace=True):
        """ Pull results from **self.csv_lig** or **self.csv_qd**.
        Performs in inplace update of **df** if **inplace** = *True*, returing *None*.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter columns: A list of to be updated columns in **df**.
        :parameter bool get_mol: Attempt to pull preexisting molecules from the database.
            See **inplace** for more details.
        :parameter bool inplace: If *True* perform an inplace update of the *mol* column in **df**.
            Otherwise Return a new series of PLAMS molecules.
        :return: If **inplace** = *False*: return a new series of PLAMS molecules pulled
            from **self**, else return |None|_
        :rtype: |None|_ or |pd.Series|_ (index: |str|_, values: |plams.Molecule|_)
        """
        # Operate on either the ligand or quantum dot database
        if database == 'ligand':
            path = self.csv_lig
            open_csv = self.open_csv_lig
        elif database == 'QD':
            path = self.csv_qd
            open_csv = self.open_csv_qd

        # Update the *hdf5 index* column in **df**
        with open_csv(path, write=False) as db:
            df.update(db, overwrite=True)
            df['hdf5 index'] = df['hdf5 index'].astype(int, copy=False)

        # **df** has been updated and **get_mol** = *False*
        if get_mol:
            ret = self._get_csv_mol(df, database, inplace)
        else:
            ret = None

        # Return a new series if **inplace** = *False*; return *None* otherwise
        return ret

    def _get_csv_mol(self, df, database='ligand', inplace=True):
        """ A method which handles the retrieval and subsequent formatting of molecules.
        Called internally by :meth:`Database.from_csv`.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool inplace: If *True* perform an inplace update of the *mol* column in **df**.
            Otherwise Return a new series of PLAMS molecules.
        :parameter bool close: If the database component should be closed afterwards.
        :return: If **inplace** = *False*: return a new series of PLAMS molecules pulled
            from **self**, else return |None|_
        :rtype: |None|_ or |pd.Series|_ (index: |str|_, values: |plams.Molecule|_)
        """
        # Sort and find all valid HDF5 indices
        df.sort_values(by=['hdf5 index'], inplace=True)
        df_slice = df['hdf5 index'] >= 0
        idx = df['hdf5 index'][df_slice].values

        # If no HDF5 indices are availble in **df** then abort the function
        if not df_slice.any():
            if inplace:
                return None
            return pd.Series(None, name=('mol', ''), dtype=object)

        # Update **df** with preexisting molecules from **self**, returning *None*
        if inplace:
            mol_list = self.from_hdf5(idx, database=database)
            for i, rdmol in zip(df_slice.index, mol_list):
                df.loc[i, ('mol', '')].from_rdmol(rdmol)
            ret = None

        # Create and return a new series of PLAMS molecules
        else:
            mol_list = self.from_hdf5(idx, database=database, rdmol=False)
            ret = pd.Series(mol_list, index=df[df_slice].index, name=('mol', ''))

        return ret

    def from_hdf5(self, index, database='ligand', rdmol=True, close=True):
        """ Import structures from the hdf5 database as RDKit or PLAMS molecules.

        :parameter index: The indices of the to be retrieved structures.
        :type index: |list|_ [|int|_]
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool rdmol: If *True*, return an RDKit molecule instead of a PLAMS molecule.
        :parameter bool close: If the database component should be closed afterwards.
        :return: A list of PLAMS or RDKit molecules.
        :rtype: |list|_ [|plams.Molecule|_ or |rdkit.Chem.Mol|_]
        """
        # Convert **index** to an array if it is a series or dataframe
        if isinstance(index, (pd.Series, pd.DataFrame)):
            index = index.values.tolist()
        elif isinstance(index, np.ndarray):
            index = index.tolist()

        # Open the database and pull entries
        with h5py.File(self.hdf5, 'r') as f:
            pdb_array = f[database][index]

        # Return a list of RDKit or PLAMS molecules
        return [from_pdb_array(mol, rdmol=rdmol) for mol in pdb_array]
