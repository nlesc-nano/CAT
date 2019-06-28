"""A module for holding functions related to the Database class."""

__all__ = ['mol_to_file', 'df_to_mongo_dict']

from os import getcwd
from os.path import (join, isfile, isdir)
from typing import (Dict, Any, List)

import yaml
import h5py
import numpy as np
import pandas as pd
from pymongo import MongoClient, ASCENDING

from scm.plams import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..mol_utils import from_rdmol
from ..utils import (get_time, get_template)


def even_index(df1: pd.DataFrame,
               df2: pd.DataFrame) -> pd.DataFrame:
    """Ensure that ``df2.index`` is a subset of ``df1.index``.

    Parameters
    ----------
    df1 : |pd.DataFrame|_
        A DataFrame whose index is to-be a superset of ``df2.index``.

    df2 : |pd.DataFrame|_
        A DataFrame whose index is to-be a subset of ``df1.index``.

    Returns
    -------
    |pd.DataFrame|_
        A new

    """
    # Figure out if ``df1.index`` is a subset of ``df2.index``
    bool_ar = df2.index.isin(df1.index)
    if bool_ar.all():
        return df1

    # Make ``df1.index`` a subset of ``df2.index``
    nan_row = get_nan_row(df1)
    idx = df2.index[~bool_ar]
    df_tmp = pd.DataFrame(len(idx) * [nan_row], index=idx, columns=df1.columns)
    return df1.append(df_tmp, sort=True)


def get_unflattend(input_dict: dict) -> zip:
    """Flatten a dictionary and return a :class:`zip` instance consisting of keys and values."""
    def _unflatten(input_dict_: dict) -> dict:
        """ """
        ret = Settings()
        for key, value in input_dict_.items():
            s = ret
            for k1, k2 in zip(key[:-1], key[1:]):
                s = s[k1]
            s[key[-1]] = value

        return ret.as_dict()

    return zip(*[(k, _unflatten(v)) for k, v in input_dict.items()])


def df_to_mongo_dict(df: pd.DataFrame) -> List[dict]:
    """Convert a dataframe into a dictionary suitable for MongoDB."""
    keys, ret = get_unflattend(df.T.to_dict())
    idx_names = df.index.names

    for item, idx in zip(ret, keys):
        idx_dict = dict(zip(idx_names, idx))
        item.update(idx_dict)

    return ret


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

    if not mol_format:
        return None

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
            if 'xyz' in mol_format and not isfile(mol_path + '.xyz'):
                mol.write(mol_path + '.xyz')


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
    dtype_dict = {
        np.dtype('int64'): -1,
        np.dtype('float64'): np.nan,
        np.dtype('O'): None,
        np.dtype('bool'): False
    }

    if not isinstance(df.index, pd.MultiIndex):
        return [dtype_dict[df[i].dtype] for i in df]
    else:
        ret = []
        for _, value in df.items():
            try:
                j = dtype_dict[value.dtype]
            except KeyError:  # dtype is neither int, float nor object
                j = None
            ret.append(j)
        return ret


def as_pdb_array(mol_list, min_size=0):  # TODO return a generator instead of an array
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

    columns = pd.MultiIndex.from_tuples(
        [('hdf5 index', ''), ('formula', ''), ('opt', ''), ('settings', 1)],
        names=['index', 'sub index']
    )

    df = pd.DataFrame(None, index=idx, columns=columns)
    df['hdf5 index'] = -1
    df['formula'] = 'str'
    df['settings'] = 'str'
    df['opt'] = False
    df.to_csv(path)


def _create_csv_qd(path):
    """ Create a QD database and and return its absolute path.

    :param str path: The path to the database.
    """
    idx = pd.MultiIndex.from_tuples(
        [('-', '-', '-', '-')],
        names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    )

    columns = pd.MultiIndex.from_tuples(
        [('hdf5 index', ''), ('ligand count', ''), ('opt', ''), ('settings', 1), ('settings', 2)],
        names=['index', 'sub index']
    )

    df = pd.DataFrame(None, index=idx, columns=columns)
    df['hdf5 index'] = -1
    df['ligand count'] = -1
    df['settings'] = 'str'
    df['opt'] = False
    df.to_csv(path)


def _create_hdf5(path, name='structures.hdf5'):
    """ Create a pdb structure database (hdf5 format), populate it with the *core*, *ligand*
    and *QD* datasets and finally return its absolute path.

    :param str path: The path to the database.
    :param str name: The filename of the database (excluding its path)
    :return: The absolute path to the pdb structure database.
    :rtype: |str|_
    """
    # Define arguments for 2D datasets
    path = join(path, name)
    dataset_names = ('core', 'core_no_opt', 'ligand', 'ligand_no_opt', 'QD', 'QD_no_opt', )
    kwarg = {'chunks': True, 'maxshape': (None, None), 'compression': 'gzip'}

    # Create new 2D datasets
    with h5py.File(path, 'a') as f:
        for name in dataset_names:
            if name not in f:
                f.create_dataset(name=name, data=np.empty((0, 1), dtype='S80'), **kwarg)

    # Define arguments for 3D datasets
    dataset_names_3d = ('job_settings_crs', 'job_settings_QD_opt', 'job_settings_BDE')
    kwarg_3d = {'chunks': True, 'maxshape': (None, None, None), 'compression': 'gzip'}

    # Create new 3D datasets
    with h5py.File(path, 'a') as f:
        for name in dataset_names_3d:
            if name not in f:
                f.create_dataset(name=name, data=np.empty((0, 1, 1), dtype='S120'), **kwarg_3d)

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


def _create_mongodb(host: str = 'localhost',
                    port: int = 27017,
                    **kwargs: Dict[str, Any]) -> dict:
    """Create the the MongoDB collections and set their index.

    Paramaters
    ----------
    host : |str|_
        Hostname or IP address or Unix domain socket path of a single mongod or
        mongos instance to connect to, or a mongodb URI, or a list of hostnames mongodb URIs.
        If **host** is an IPv6 literal it must be enclosed in ``"["`` and ``"["`` characters
        following the RFC2732 URL syntax (e.g. ``"[::1]"`` for localhost).
        Multihomed and round robin DNS addresses are not supported.

    port : |str|_
        port number on which to connect.

    kwargs : |dict|_
        Optional keyword argument for `pymongo.MongoClient <http://api.mongodb.com/python/current/api/pymongo/mongo_client.html>`_.  # noqa

    Returns
    -------
    |dict|_
        A dictionary with all supplied keyword arguments.

    Raises
    ------
    ServerSelectionTimeoutError
        Raised if no connection can be established with the host.

    """
    # Open the client
    client = MongoClient(host, port, serverSelectionTimeoutMS=5000, **kwargs)
    client.server_info()  # Raises an ServerSelectionTimeoutError error if the server is inaccesible

    # Open the database
    db = client.cat_database

    # Open and set the index of the ligand collection
    lig_collection = db.ligand_database
    lig_key = 'smiles_1_anchor_1'
    if lig_key not in lig_collection.index_information():
        lig_collection.create_index([
            ('smiles', ASCENDING),
            ('anchor', ASCENDING)
        ], unique=True)

    # Open and set the index of the QD collection
    qd_collection = db.QD_database
    qd_key = 'core_1_core anchor_1_ligand smiles_1_ligand anchor_1'
    if qd_key not in qd_collection.index_information():
        qd_collection.create_index([
            ('core', ASCENDING),
            ('core anchor', ASCENDING),
            ('ligand smiles', ASCENDING),
            ('ligand anchor', ASCENDING)
        ], unique=True)

    # Return all provided keyword argument
    ret = {'host': host, 'port': port}
    ret.update(kwargs)
    return ret
