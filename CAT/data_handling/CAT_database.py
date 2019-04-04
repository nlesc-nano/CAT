""" A module which holds the Database class. """

__all__ = ['Database', 'mol_to_file']

from os import getcwd
from os.path import (join, isfile, isdir)

import yaml
import h5py
import numpy as np
import pandas as pd

from scm.plams.core.settings import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from CAT import utils as CAT
from CAT.mol_utils import from_rdmol
from CAT.utils import (get_time, type_to_string)


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

    # Raise an error if mol_format is empty
    if not mol_format:
        return

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
    try:  # Plan A
        return [dtype_dict[df[i].dtype] for i in df]
    except KeyError:  # Plan B (for when dealing with multiindices)
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


def sanitize_yaml_settings(s, job_type):
    """ Remove a predetermined set of unwanted keys and values from a settings object.

    :param s: A settings object with, potentially, undesired keys and values.
    :type s: |plams.Settings|_ (superclass: |dict|_)
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
    blacklist = CAT.get_template('settings_blacklist.yaml')
    s_del = blacklist['generic']
    s_del.update(blacklist[job_type])

    # Recursivelly delete all keys from **s** if aforementioned keys are present in the s_del
    recursive_del(s, s_del)
    return s


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
    path = join(path, name)
    hdf5 = h5py.File(path, 'a')

    # Check if the dataset is already available in structures.hdf5
    dataset_names = 'core', 'ligand', 'QD'
    for name in dataset_names:
        if name not in hdf5:
            hdf5.create_dataset(name=name, data=np.empty((0, 1), dtype='S80'),
                                chunks=True, maxshape=(None, None), compression='gzip')
    hdf5.close()
    return path


class Database():
    """ The Database class.

    :Atributes:     * **path** (|plams.Settings|_) – A settings object with the absolute paths to \
                    the various database components (see below).

                    * **csv_lig** (|None|_ or |pd.DataFrame|_) – A dataframe with all ligand \
                    related results.

                    * **csv_qd** (|None|_ or |pd.DataFrame|_) – A dataframe with all quantum dot \
                    related results.

                    * **yaml** (|None|_ or |plams.Settings|_) – A settings object with all job \
                    settings.

                    * **hdf5** (|None|_ or |h5py.File|_) – A file object with all \
                    (partially) deserialized .pdb files.

                    * **mongodb** (|None|_) – *None*.
    """

    def __init__(self, path=None):
        path = path or getcwd()

        # Attributes which hold the absolute paths to various components of the database
        self.path = Settings()
        self.path.csv_lig = _create_csv(path, database='ligand')
        self.path.csv_qd = _create_csv(path, database='QD')
        self.path.yaml = join(path, 'job_settings.yaml')
        self.path.hdf5 = _create_hdf5(path)
        self.path.mongodb = None  # Placeholder

        # Attributes which hold the actual components of the database (when opened)
        self.csv_lig = None
        self.csv_qd = None
        self.yaml = None
        self.hdf5 = None
        self.mongodb = None  # Placeholder

    def __str__(self):
        ret = Settings()
        attr_dict = vars(self)
        for key in attr_dict:
            ret[key] = type(attr_dict[key])
        return str(ret)

    """ #################################  Opening the database ############################### """

    def open_all(self, database=None):
        """ Open all components of the database, populating their respective attributes:

            * :func:`Database.open_csv`

            * :func:`Database.open_yaml`

            * :func:`Database.open_hdf5`

        :parameter database: If not *None*, open all components of a specific database.
            Accepted values (besides *None*) are *ligand* and *QD*.
        :type database: |None|_ or |str|_
        """
        self.open_yaml()
        self.open_hdf5()
        if database is None:
            self.open_csv(database='ligand')
            self.open_csv(database='QD')
        else:
            self.open_csv(database=database)

    def open_csv(self, database='ligand'):
        """ Open the ligand or quantum dot database, populating **self.csv_lig** or **self.csv_qd**
        with a |pd.DataFrame|_ object.

        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        """
        if database == 'ligand':
            self._open_csv_lig()
        elif database == 'QD':
            self._open_csv_qd()
        else:
            raise ValueError('open_csv: {} is not an accepted value'.format(str(database)))

    def _open_csv_lig(self):
        """ Open the ligand database, populating **self.csv_lig** with a |pd.DataFrame|_ object. """
        if self.csv_lig is None:
            dtype = {'hdf5 index': int, 'formula': str, 'settings': str}
            self.csv_lig = pd.read_csv(self.path.csv_lig, index_col=[0, 1],
                                       header=[0, 1], dtype=dtype)
            idx_tups = [(i, '') if 'Unnamed' in j else (i, j) for i, j in self.csv_lig.columns]
            columns = pd.MultiIndex.from_tuples(idx_tups, names=self.csv_lig.columns.names)
            self.csv_lig.columns = columns

    def _open_csv_qd(self):
        """ Open the quantum dot database, populating **self.csv_qd**
        with a |pd.DataFrame|_ object.
        """
        if self.csv_qd is None:
            dtype = {'hdf5 index': int, 'ligand count': np.int64, 'settings': str}
            self.csv_qd = pd.read_csv(self.path.csv_qd, index_col=[0, 1, 2, 3],
                                      header=[0, 1], dtype=dtype)
            idx_tups = [(i, '') if 'Unnamed' in j else (i, j) for i, j in self.csv_qd.columns]
            columns = pd.MultiIndex.from_tuples(idx_tups, names=self.csv_qd.columns.names)
            self.csv_qd.columns = columns

    def open_yaml(self):
        """ Open the job settings database, populating **self.yaml** with a
        |plams.Settings|_ object.
        """
        if self.yaml is None:
            if isfile(self.path.yaml):
                with open(self.path.yaml, 'r') as file:
                    self.yaml = Settings(yaml.load(file, Loader=yaml.FullLoader))
            else:
                self.yaml = Settings()

    def open_hdf5(self):
        """ Open the .pdb structure database, populating **self.hdf5**
        with a |h5py.File|_ object.
        """
        if self.hdf5 is None:
            self.hdf5 = h5py.File(self.path.hdf5, 'a')

    """ #################################  Closing the database ############################### """

    def close_all(self, write=True):
        """ Close all components of the database, depopulating their respective attributes and
        exporting its content:

            * :func:`Database.close_csv`

            * :func:`Database.close_yaml`

            * :func:`Database.close_hdf5`

        :parameter bool write: Whether or not the database content should be exported.
        """
        self.close_csv(database='ligand', write=write)
        self.close_csv(database='QD', write=write)
        self.close_yaml(write=write)
        self.close_hdf5()

    def close_csv(self, database='ligand', write=True):
        """ Close the ligand or quantum dot database, depopulating **self.csv_lig** or
        **self.csv_qd** and exporting its content to a .csv file.

        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool write: Whether or not **self.csv_lig** or **self.csv_qd** should be
            exported to a .csv file.
        """
        if database == 'ligand':
            self._close_csv_lig(write)
        elif database == 'QD':
            self._close_csv_qd(write)
        else:
            raise ValueError('close_csv: {} is not an accepted value'.format(str(database)))

    def _close_csv_lig(self, write=True):
        """ Close the ligand database, depopulating **self.csv_lig** and exporting its content to
        a .csv file.

        :parameter bool write: Whether or not **self.csv_lig** should be exported to a .csv file.
        """
        if self.csv_lig is not None:
            if write:
                self.csv_lig.to_csv(self.path.csv_lig)
            self.csv_lig = None

    def _close_csv_qd(self, write=True):
        """ Close the quantum dot database, depopulating **self.csv_qd** and exporting its content
        to a .csv file.

        :parameter bool write: Whether or not **self.csv_qd** should be exported to a .csv file.
        """
        if self.csv_qd is not None:
            if write:
                self.csv_qd.to_csv(self.path.csv_qd)
            self.csv_qd = None

    def close_yaml(self, write=True):
        """ Close the job settings database, depopulating **self.yaml** and exporting
        its content to a .yaml file.

        :parameter bool write: Whether or not **self.yaml** should be exported to a .yaml file.
        """
        if self.yaml is not None:
            if write:
                yml_dict = self.yaml.as_dict()
                for key in yml_dict:
                    for i, value in enumerate(yml_dict[key]):
                        if isinstance(value, Settings):
                            yml_dict[key][i] = value.as_dict()

                with open(self.path.yaml, 'w') as file:
                    file.write(yaml.dump(yml_dict, default_flow_style=False, indent=4))
            self.yaml = None

    def close_hdf5(self):
        """ Close the .pdb structure database, depopulating **self.hdf5** and closing
        the .hdf5 file.
        """
        if self.hdf5 is not None:
            self.hdf5.close()
            self.hdf5 = None

    """ #################################  Updating the database ############################## """

    def update_csv(self, df, database='ligand', columns=None, overwrite=False, job_recipe=None,
                   close=True):
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
        :parameter bool close: If the database component should be closed afterwards.
        """
        # Open the database
        self.open_csv(database)

        # Operate on either *self.csv_lig* or *self.csv_qd*
        if database == 'ligand':
            csv = self.csv_lig
        elif database == 'QD':
            csv = self.csv_qd

        # Update **self.yaml**
        if job_recipe is not None:
            job_settings = self.update_yaml(job_recipe)
            for key in job_settings:
                df[('settings', key)] = job_settings[key]

        # Update **csv.index**
        nan_row = get_nan_row(csv)
        for i in df.index:
            if i not in csv.index:
                csv.at[i, :] = nan_row
        csv['hdf5 index'] = csv['hdf5 index'].astype(int, copy=False)  # Fix the data type

        # Filter columns
        if not columns:
            df_columns = df.columns
        else:
            df_columns = columns + [i for i in df.columns if i[0] == 'settings']

        # Update **csv.columns**
        for i in df_columns:
            if i not in csv.columns:
                try:
                    csv[i] = np.array((None), dtype=df[i].dtype)
                except TypeError:  # i.e. if csv[i] consists of np.int64
                    csv[i] = -1

        # Update **self.hdf5**; returns a new series of indices
        hdf5_series = self.update_hdf5(df, database=database, overwrite=overwrite)

        # Update **csv.values**
        csv.update(df[columns], overwrite=overwrite)
        csv.update(hdf5_series, overwrite=True)

        # Close the database
        if close:
            self.close_csv(database)

    def update_yaml(self, job_recipe, close=True):
        """ Update **self.yaml** with (potentially) new user provided settings.

        :parameter job_recipe: A settings object with one or more settings specific to a job.
        :type job_recipe: |plams.Settings|_ (superclass: |dict|_)
        :parameter bool close: If the database component should be closed afterwards.
        :return: A dictionary with the column names as keys and the key for **self.yaml** as
            matching values.
        :rtype: |dict|_ (keys: |str|_, values: |str|_)
        """
        if self.yaml is None:
            self.open_yaml()

        ret = {}
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
            if key not in self.yaml:
                self.yaml[key] = []

            # Check if the appropiate value is available in **self.yaml**
            if value in self.yaml[key]:
                ret[item] = key + ' ' + str(self.yaml[key].index(value))
            else:
                self.yaml[key].append(value)
                ret[item] = key + ' ' + str(len(self.yaml[key]) - 1)

        if close:
            self.close_yaml()
        return ret

    def update_hdf5(self, df, database='ligand', overwrite=False, close=True):
        """ Export molecules (see the *mol* column in **df**) to the structure database.
        Returns a series with the **self.hdf5** indices of all new entries.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool overwrite: Whether or not previous entries can be overwritten or not.
        :parameter bool close: If the database component should be closed afterwards.
        :return: A series with the index of all new molecules in **self.hdf5**
        :rtype: |pd.Series|_ (index: |str|_, values: |np.int64|_)
        """
        self.open_hdf5()

        # Identify new and preexisting entries
        new = df['hdf5 index'][df['hdf5 index'] == -1]
        old = df['hdf5 index'][df['hdf5 index'] != -1]

        i, j = self.hdf5[database].shape

        # Add new entries to the database
        if new.any():
            pdb_array = as_pdb_array(df['mol'][new.index], min_size=j)

            # Reshape and update **self.hdf5**
            k = i + pdb_array.shape[0]
            self.hdf5[database].shape = k, pdb_array.shape[1]
            self.hdf5[database][i:k] = pdb_array
            ret = pd.Series(np.arange(i, k), index=new.index, name=('hdf5 index', ''))
        else:
            ret = pd.Series(name=('hdf5 index', ''), dtype=int)

        # If **overwrite** is *True*
        if overwrite and old.any():
            self.hdf5[database][old] = as_pdb_array(df['mol'][old.index], min_size=j)

        if close:
            self.close_hdf5()

        return ret

    """ ########################  Pulling results from the database ########################### """

    def from_csv(self, df, database='ligand', get_mol=True, inplace=True, close=True):
        """ Pull results from **self.csv_lig** or **self.csv_qd**.
        Performs in inplace update of **df** if **inplace** = *True*.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter columns: A list of to be updated columns in **df**.
        :parameter bool get_mol: Attempt to pull preexisting molecules from the database.
            See **inplace** for more details.
        :parameter bool inplace: If *True* perform an inplace update of the *mol* column in **df**.
            Otherwise Return a new series of PLAMS molecules.
        :parameter bool close: If the database component should be closed afterwards.
        :return: If **inplace** = *False*: return a new series of PLAMS molecules pulled
            from **self**.
        :rtype: |pd.Series|_ (index: |str|_, values: |plams.Molecule|_)
        """
        # Open the database
        self.open_csv(database)

        # Operate on either **self.csv_lig** or **self.csv_qd**
        if database == 'ligand':
            csv = self.csv_lig
        elif database == 'QD':
            csv = self.csv_qd

        # Update the *hdf5 index* column in **df**
        df.update(csv, overwrite=True)
        df['hdf5 index'] = df['hdf5 index'].astype(int, copy=False)

        # **df** has been updated and **get_mol** = *False*
        if get_mol:
            ret = self._from_csv_mol()

        if close:
            self.close_csv(database, write=False)

        # hdf5 fancy indexing only supports use of sorted lists of indices
        df.sort_values(by=['hdf5 index'], inplace=True)

        # Close the database
        if close:
            self.close_csv(database, write=False)

        # Return a new series if **inplace** = *False*
        if not inplace:
            if df_slice.any():
                return ret
            return pd.Series(None, name=('mol', ''), dtype=object)

    def _from_csv_mol(self):
        """ """
        # Update the *mol* column in **df** or return a new series
        df_slice = df['hdf5 index'] >= 0
        if df_slice.any():
            idx = df['hdf5 index'][df_slice].values
            if inplace:  # Update **df** with preexisting molecules from **self**
                mol_list = self.from_hdf5(idx, database=database)
                for i, rdmol in zip(df_slice.index, mol_list):
                    df.loc[i, ('mol', '')].from_rdmol(rdmol)
            else:  # Create and return a new series of PLAMS molecules
                mol_list = self.from_hdf5(idx, database=database, rdmol=False)
                ret = pd.Series(mol_list, index=df[df_slice].index, name=('mol', ''))

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
        # Open the database
        self.open_hdf5()

        # Convert **index** to an array if it is a series or dataframe
        if isinstance(index, (pd.Series, pd.DataFrame)):
            index = index.values.tolist()
        elif isinstance(index, np.ndarray):
            index = index.tolist()

        # Pull entries from the database as a list of RDKit or PLAMS molecules
        pdb_array = self.hdf5[database][index]
        ret = [from_pdb_array(mol, rdmol=rdmol) for mol in pdb_array]

        # Close the database
        if close:
            self.close_hdf5()
        return ret

    """ #####################################  Utilities  #################################### """
