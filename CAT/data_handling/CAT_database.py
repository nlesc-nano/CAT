""" A module which holds the Database class. """

__all__ = ['Database']

from os import getcwd
from os.path import (join, isfile)

import yaml
import h5py
import numpy as np
import pandas as pd

from scm.plams.core.settings import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from qmflows.templates import templates as qmflows

from rdkit import Chem

from CAT.mol_utils import from_rdmol
from CAT.utils import (get_time, type_to_string)


def as_pdb_array(mol_list, min_size=0):
    """ Converts a PLAMS molecule into an array of strings, the array consisting of a
    (partially) de-serialized .pdb file.

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
    """ Converts a an array with a (partially) de-serialized .pdb file into an
    RDKit or PLAMS molecule.

    :parameter array: A (partially) de-serialized .pdb file.
    :type array: |np.ndarray|_ [|np.bytes|_ / *|S80*]
    :parameter bool rdmol: If *True*, return an RDKit molecule instead of a PLAMS molecule.
    :return: A PLAMS or RDKit molecule build from **array**.
    :rtype: |plams.Molecule|_ or |rdkit.Chem.Mol|_
    """
    pdb_str = ''.join([item.decode() + '\n' for item in array if item])
    rdmol = Chem.MolFromPDBBlock(pdb_str, removeHs=False, proximityBonding=False)
    if not rdmol:
        return molkit.from_rdmol(rdmol)
    return rdmol


def _sanitize_yaml_settings(s):
    """ Remove a predetermined set of unwanted keys and values from a settings object.

    :param s: A settings object with, potentially, undesired keys and values.
    :type s: |plams.Settings|_ (superclass: |dict|_)
    :return: A (nested) dictionary with unwanted keys and values removed.
    :rtype: |dict|_.
    """
    s.description = s.ignore_molecule = s.input.ams = s.input.Compound = s.input.compound._h = None
    del s.description
    del s.ignore_molecule
    del s.input.ams
    del s.input.Compound
    del s.input.compound._h
    if not s.input.compound:
        del s.input.compound
    return s.as_dict()


def _create_csv(path, database='ligand'):
    """ Create a ligand or QD database (csv format) and, if it does not exist, and return
    its absolute path.

    :param str path: The path to the database.
    :param str database: The type of database, accepted values are *ligand* and *qd*.
    :return: The absolute path to the ligand or QD database.
    :rtype: |str|_.
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
            raise TypeError()
    return path


def _create_csv_lig(path):
    """ Create a ligand database and and return its absolute path.

    :param str path: The path to the database.
    """
    idx = sorted(['hdf5 index', 'formula'])
    idx = pd.MultiIndex.from_tuples([(i, '') for i in idx], names=['index', 'sub index'])
    columns = pd.MultiIndex.from_tuples([(None, None)], names=['smiles', 'anchor'])
    df = pd.DataFrame(None, index=idx, columns=columns)
    df.to_csv(path)


def _create_csv_qd(path):
    """ Create a QD database and and return its absolute path.

    :param str path: The path to the database.
    """
    idx = sorted(['hdf5 index', 'ligand count'])
    idx = pd.MultiIndex.from_tuples([(i, '') for i in idx], names=['index', 'sub index'])
    columns = pd.MultiIndex.from_tuples(
            [(None, None, None, None)],
            names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    )
    df = pd.DataFrame(None, index=idx, columns=columns)
    df.to_csv(path)


def _create_hdf5(path, name='structures.hdf5'):
    """ Create a pdb structure database (hdf5 format), populate it with the *core*, *ligand*
    and *QD* datasets and finally return its absolute path.

    :param str path: The path to the database.
    :param str name: The filename of the database (excluding its path)
    :return: The absolute path to the pdb structure database.
    :rtype: |str|_.
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

                    * **csv_lig** (|None|_ or |pd.DataFrame|_) – A dataframe with all results \
                    related to the ligands.

                    * **csv_qd** (|None|_ or |pd.DataFrame|_) – A dataframe with all results \
                    related to the quantum dots.

                    * **yaml** (|None|_ or |plams.Settings|_) – A settings object with all job \
                    settings.

                    * **hdf5** (|None|_ or |h5py.File|_) – A file object with all \
                    (partially) deserialized .pdb files.

                    * **mongodb** (|None|_) – *None*.
    """

    def __init__(self, path=None):
        """ """
        path = path or getcwd()

        # Attributes which hold the absolute paths to various components of the database
        self.path = Settings()
        self.path.csv_lig = _create_csv(path, database='ligand')
        self.path.csv_qd = _create_csv(path, database='QD')
        self.path.yaml = join(path, 'job_settings.yaml')
        self.path.hdf5 = _create_hdf5(path)
        self.path.mongodb = None

        # Attributes which hold the actual components of the database (when opened)
        self.csv_lig = None
        self.csv_qd = None
        self.yaml = None
        self.hdf5 = None
        self.mongodb = None

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
            raise KeyError()

    def _open_csv_lig(self):
        """ Open the ligand database, populating **self.csv_lig** with a |pd.DataFrame|_ object. """
        if self.csv_lig is None:
            self.csv_lig = pd.read_csv(self.path.csv_lig, index_col=[0, 1],
                                       header=[0, 1], keep_default_na=False)
            self.csv_lig.replace('', np.nan, inplace=True)

    def _open_csv_qd(self):
        """ Open the quantum dot database, populating **self.csv_qd**
        with a |pd.DataFrame|_ object.
        """
        if self.csv_qd is None:
            self.csv_qd = pd.read_csv(self.path.csv_qd, index_col=[0, 1],
                                      header=[0, 1, 2, 3], keep_default_na=False)
            self.csv_qd.replace('', np.nan, inplace=True)

    def open_yaml(self):
        """ Open the job settings database, populating **self.yaml** with a
        |plams.Settings|_ object.
        """
        if self.yaml is None:
            if isfile(self.path.yaml):
                with open(self.path.yaml, 'r') as file:
                    self.yaml = Settings(yaml.load(file))
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
        :parameter bool write: Whether or not **self.csv_lig** / **self.csv_qd** should be
            exported to a .csv file.
        """
        if database == 'ligand':
            self._close_csv_lig(write)
        elif database == 'QD':
            self._close_csv_qd(write)
        else:
            raise KeyError()

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

        :parameter bool write: Whether or not **self.csv_lig** should be exported to a .csv file.
        """
        if self.yaml is not None:
            if write:
                with open(self.path.yaml, 'w') as file:
                    file.write(yaml.dump(self.yaml.as_dict(), default_flow_style=False, indent=4))
            self.yaml = None

    def close_hdf5(self):
        """ Close the .pdb structure database, depopulating **self.hdf5** and closing
        the .hdf5 file.
        """
        if self.hdf5 is not None:
            self.hdf5.close()
            self.hdf5 = None

    """ #################################  Updating the database ############################## """

    def update_csv(self, df, database='ligand', columns=None, overwrite=False, close=True,
                   job_recipe=None):
        """ Update **self.csv_lig** or **self.csv_qd** with
        (potentially) new user provided settings.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter columns: An array with *n* column keys in **df** which
            (potentially) are to be added to **self**. If *None*: Add all columns.
        :type columns: |None|_ or *n*2* |np.ndarray|_ [|np.str_|_]
        :parameter bool overwrite: Whether or not previous entries can be overwritten or not.
        :parameter bool close: If the database should be closed afterwards.
        :parameter job_recipe: A Settings object with settings specific to a job.
        :type job_recipe: |None|_ or |plams.Settings|_ (superclass: |dict|_)
        """
        # Open the database
        self.open_csv(database)

        # Operate on either *self.csv_lig* or *self.csv_qd*
        if database == 'ligand':
            csv = self.csv_lig
        elif database == 'QD':
            csv = self.csv_qd
        else:
            raise TypeError()

        # Update job settings
        if job_recipe is not None:
            job_settings = self.update_yaml(job_recipe)
            for key in job_settings:
                df[key] = job_settings[key]

        # Update columns
        idx_array = np.array(df.index)
        for i in idx_array[np.isin(df.index, csv.columns, invert=True)]:
            csv[i] = np.nan

        # Update index
        if columns is None:
            columns = df.columns
        else:
            columns += list(job_settings.keys())
            columns = np.array(columns)
        for i in columns:
            if i not in csv.index:
                csv.loc[i, :] = None

        # Update hdf5 values
        hdf5_series = self.update_hdf5(df, database=database, overwrite=overwrite, close=False)
        df.update(hdf5_series, overwrite=True)

        # Update csv values
        csv.update(df.T, overwrite=overwrite)

        # Close the database
        if close:
            self.close_csv(database)

    def update_yaml(self, job_recipe, close=True):
        """ Update **self.yaml** with (potentially) new user provided settings.

        :parameter job_recipe: A settings object with one or more settings specific to a job.
        :type job_recipe: |plams.Settings|_ (superclass: |dict|_)
        :parameter bool close: If the job settings database should be closed afterwards.
        :return: A dictionary with the row name of **self.csv_lig** or **self.csv_qd** as keys
            and the key for **self.yaml** as matching values.
        :rtype: |dict|_ (keys: |str|_, values: |str|_).
        """
        if self.yaml is None:
            self.open_yaml()

        ret = {}
        for item in job_recipe:
            # Prepare keys and values
            if not isinstance(job_recipe[item].key, str):
                template = job_recipe[item].template
                key = type_to_string(job_recipe[item].key)
                value = qmflows.get_template(template)['specific'][key]
                value.update(_sanitize_yaml_settings(job_recipe[item].value))
                key = key.rsplit('.', 1)[-1].split("'")[0]
            else:
                key = job_recipe[item].key
                value = job_recipe[item].value
            name = job_recipe[item].name

            # Update ret
            if key not in self.yaml:
                self.yaml[key] = []
            if value in self.yaml[key]:
                ret[name] = key + ' ' + str(self.yaml[key].index(value))
            else:
                self.yaml[key].append(value)
                ret[name] = key + ' ' + str(len(self.yaml[key]) - 1)

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
        :parameter bool close: If the database should be closed afterwards.
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
            pdb_array = as_pdb_array(df['mol'][new], min_size=j)

            # Reshape and update **self.hdf5**
            self.hdf5[database].shape = i + pdb_array.shape[0], pdb_array.shape[1]
            self.hdf5[database][i:i+pdb_array.shape[0]] = pdb_array
            ret = pd.Series(np.arange(i, i + pdb_array.shape[0]),
                            index=new.index, name=('hdf5 index', ''))

        # If **overwrite** is *True*
        if overwrite and old.any():
            self.hdf5[database][old] = as_pdb_array(df['mol'][old], min_size=j)

        if close:
            self.close_hdf5()
            self.close_csv(database)

        return ret

    """ ########################  Pulling results from the database ########################### """

    def from_csv(self, df, database='ligand', columns=None, close=True, inplace=True):
        """ Pull results from **self.csv_lig** or **self.csv_qd**.
        Performs in inplace update of **df**.

        :parameter df: A dataframe of new (potential) database entries.
        :type df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool close: If the database should be closed afterwards.
        :parameter bool inplace: If *True* perform an inplace update of the *mol* column in **df**.
            Otherwise Return a new series of PLAMS molecules.
        :return: If **inplace** = *False*: return a new series of PLAMS molecules pulled
            from **self**.
        :rtype: |pd.Series|_ (index: |str|, values: |plams.Molecule|_, name=|str|_)
        """
        # Open the database
        self.open_csv(database)

        # Operate on either *self.csv_lig* or *self.csv_qd*
        if database == 'ligand':
            csv = self.csv_lig
        elif database == 'QD':
            csv = self.csv_qd
        else:
            raise TypeError()

        # Attempt to update **df** with preexisting .pdb files from **self**
        # Or create and return a new series of PLAMS molecules
        columns = list(df.loc[np.isin(df.index, csv.columns)].index)
        if columns:
            hdf5_idx = csv[columns]['hdf5 index']
            if inplace:
                mol_list = self.from_hdf5(hdf5_idx)
                for i, rdmol in zip(columns, mol_list):
                    df.at[i, 'mol'].from_rdmol(rdmol)
            else:
                mol_list = self.from_hdf5(hdf5_idx, rdmol=False)
                idx = pd.MultiIndex.from_tuples(columns, names=df.index.names)
                ret = pd.Series(mol_list, index=idx, name=('mol', ''))

        # Close the database
        if close:
            self.close_csv(database, write=False)

        # Return a new series if **inplace** = *False*
        if not inplace:
            if columns:
                return ret
            return pd.Series(np.empty((0), dtype=object), name=('mol', ''))

    def from_hdf5(self, index, database='ligand', rdmol=True, close=False):
        """ Import structures from the hdf5 database as RDKit or PLAMS molecules.

        :parameter index: The indices of the to be retrieved structures.
        :type index: |list|_ [|int|_]
        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool rdmol: If *True*, return an RDKit molecule instead of a PLAMS molecule.
        :parameter bool close: If the database should be closed afterwards.
        :return: A list of PLAMS or RDKit molecules.
        :rtype: |list|_ [|plams.Molecule|_] or |list|_ [|rdkit.Chem.Mol|_]
        """
        # Open the database
        self.open_hdf5()

        # Pull entries from the database as a list of RDKit or PLAMS molecules
        pdb_array = self.hdf5[database][index]
        ret = from_pdb_array(pdb_array, rdmol=rdmol)

        # Close the database
        if close:
            self.close_hdf5()
        return ret
