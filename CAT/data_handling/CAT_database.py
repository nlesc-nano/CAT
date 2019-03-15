""" A module which holds the Database class. """

__all__ = ['Database']

from os import get_cwd
from os.path import (join, isfile)

import yaml
import h5py
import numpy as np
import pandas as pd

from scm.plams.core.settings import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from qmflows.templates import templates as qmflows

from rdkit import Chem

from CAT.utils import (get_time, type_to_string)


def as_pdb_array(mol_list, min_size=0):
    """ Converts a PLAMS molecule into an array of strings, the array consisting of a
    (partially) de-serialized .pdb file.

    :parameter mol_list: A list of PLAMS molecules.
    :type mol_list: |list|_ [|plams.Molecule|_]
    :parameter int min_size: The minimumum length of the pdb_array. The array is padded with empty
        strings if required.
    :return: An array with *m* partially deserialized .pdb files with up to *n* lines each.
    :rtype: *m*n* |np.ndarray|_ [|np.bytes|_ / *|S80*]
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


def from_pdb_array(array):
    """ Converts a an array with a (partially) de-serialized .pdb file into a PLAMS molecule.

    :parameter array: A (partially) de-serialized .pdb file.
    :type array: |np.ndarray|_ [|np.bytes|_ / *|S80*]
    :return: A PLAMS molecule build from **array**.
    :rtype: |plams.Molecule|_
    """
    pdb_str = ''.join([item.decode() + '\n' for item in array if item])
    return molkit.from_rdmol(Chem.MolFromPDBBlock(pdb_str, removeHs=False, proximityBonding=False))


def _sanitize_yaml_settings(s):
    """ Remove a predetermined set of unwanted keys and values from a settings object.

    :param s: A settings object with, potentially, undesired keys and values.
    :type s: |Settings|_
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


def _job_recipe_to_yaml(job_recipe):
    """
    Convert **job_recipe** into a more suitable form for storage in a *Database* object

    :parameter job_recipe: A list of one or more settings specific to a job.
    :type job_recipe: |list|_ [|Settings|_]
    :return: A list of tuples, each tuple containg a key for *Database.yaml* and matching value.
    :rtype: |list|_ [|tuple|_ [|str|_, |Settings|_]]
    """
    file_dict = {
            'ligand.optimize': ('geometry.json', None),
            'qd.optimize': ('geometry.json', 'geometry.json'),
            'qd.dissociate': ('singlepoint.json', 'freq.json')}

    ret = []
    for job_dict in job_recipe:
        if job_dict.job and job_dict.s:
            if job_dict.job_type is not None:
                template_key = type_to_string(job_recipe.job_type)
                s = qmflows.get_template(file_dict[job_dict.job_type])['specific'][template_key]
            else:
                s = Settings()
            s.update(job_dict.s)
            tmp = _sanitize_yaml_settings(s), str(job_dict.job).rsplit('.', 1)[-1].split("'")[0]
            ret.append(tmp)
        else:
            ret.append((None, None))


class Database():
    """ """
    def __init__(self, path=None):
        """ """
        path = path or get_cwd()

        # Attributes which hold the absolute paths to various components of the database
        self.path = Settings()
        self.path.lig_csv = self.create_csv(path, 'ligand')
        self.path.qd_csv = self.create_csv(path, 'QD')
        self.path.yaml = self.create_yaml(path)
        self.path.hdf5 = self.create_hdf5(path)
        self.path.mongodb = None

        # Attributes which hold the actual components of the database (when opened)
        self.lig_csv = None
        self.qd_csv = None
        self.yaml = None
        self.hdf5 = None
        self.mongodb = None

    """ #################################  Update **self.yaml** ############################### """

    def update_yaml(self, job_recipe):
        """ Update **self.yaml** with (potentially) new user provided settings.

        :parameter job_recipe: A list of one or more settings specific to a job.
        :type job_recipe: |list|_ [|Settings|_]
        :return: A tuple of keys, used for linking an entry in **self.csv** to a specific job
            setting in **self.yaml**.
        :rtype: |list|_ [|str|_].
        """
        self.open_yaml()
        yaml_input = _job_recipe_to_yaml()

        # Update the job settings database
        idx = []
        for job, s in yaml_input:
            if None in (job, s):
                idx.append(None)
            else:
                if job not in self.yaml:
                    self.yaml[job] = []
                if s in self.yaml[job]:
                    idx.append(job + ' ' + str(self.yaml[job].index(s)))
                else:
                    self.yaml[job].append(s)
                    idx.append(job + ' ' + str(len(self.yaml[job]) - 1))

        self.close_yaml()
        return idx

    """ #################################  Opening the database ############################### """

    def open_csv(self, database='ligand'):
        """ Open the ligand or quantum dot database, populating **self.lig_csv** or **self.qd_csv**
        with a |pd.DataFrame|_ object.

        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        """
        if database == 'ligand':
            self._open_csv_lig()
        elif database == 'QD':
            self._open_csv_qd()
        raise KeyError()

    def _open_csv_lig(self):
        """ Open the ligand database, populating **self.lig_csv** with a |pd.DataFrame|_ object. """
        self.csv_lig = pd.read_csv(self.path.lig_csv, index_col=[0, 1],
                                   header=[0, 1], keep_default_na=False)

    def _open_csv_qd(self):
        """ Open the quantum dot database, populating **self.qd_csv**
        with a |pd.DataFrame|_ object.
        """
        self.csv_qd = pd.read_csv(self.path.qd_csv, index_col=[0, 1],
                                  header=[0, 1, 2, 3], keep_default_na=False)

    def open_yaml(self):
        """ Open the job settings database, populating **self.yaml** with a |Settings|_ object. """
        self.yaml = Settings(yaml.load(self.path.yaml))

    def open_hdf5(self):
        """ Open the .pdb structure database, populating **self.hdf5**
        with a |h5py.File|_ object.
        """
        self.hdf5 = h5py.File(self.path.hdf5, 'a')

    def open_all(self):
        """ Open all components of the database, populating their respective attributes. """
        self.open_csv(database='ligand')
        self.open_csv(database='QD')
        self.open_yaml()
        self.open_hdf5()

    """ #################################  Closing the database ############################### """

    def close_csv(self, database='ligand', write=True):
        """ Close the ligand or quantum dot database, depopulating **self.lig_csv** or
        **self.qd_csv** and exporting its content to a .csv file.

        :parameter str database: The type of database; accepted values are *ligand* and *QD*.
        :parameter bool write: Whether or not **self.lig_csv** / **self.qd_csv** should be
            exported to a .csv file.
        """
        if database == 'ligand':
            self._close_csv_lig(write)
        elif database == 'QD':
            self._close_csv_qd(write)
        else:
            raise KeyError()

    def _close_csv_lig(self, write=True):
        """ Close the ligand database, depopulating **self.lig_csv** and exporting its content to
        a .csv file.

        :parameter bool write: Whether or not **self.lig_csv** should be exported to a .csv file.
        """
        if write:
            self.csv_lig.to_csv(self.path.csv_lig)
        self.csv_lig = None

    def _close_csv_qd(self, write=True):
        """ Close the quantum dot database, depopulating **self.qd_csv** and exporting its content
        to a .csv file.

        :parameter bool write: Whether or not **self.qd_csv** should be exported to a .csv file.
        """
        if write:
            self.csv_qd.to_csv(self.path.csv_qd)
        self.csv_qd = None

    def close_yaml(self, write=True):
        """ Close the job settings database, depopulating **self.yaml** and exporting
        its content to a .yaml file.

        :parameter bool write: Whether or not **self.lig_csv** should be exported to a .csv file.
        """
        if write:
            with open(self.path.yaml, 'w') as file:
                file.write(yaml.dump(self.yaml, default_flow_style=False, indent=4))
        self.yaml = None

    def close_hdf5(self):
        """ Close the .pdb structure database, depopulating **self.hdf5** and closing
        the .hdf5 file.
        """
        self.hdf5.close()
        self.hdf5 = None

    def close_all(self, write=True):
        """ Close all components of the database, depopulating their respective attributes and
        exporting its content.

        :parameter bool write: Whether or not the database content should be exported.
        """
        self.close_csv(database='ligand', write)
        self.close_csv(database='QD', write)
        self.close_yaml(write)
        self.close_hdf5()

    """ ##################################  Creating the database ############################# """

    def create_csv(self, path, database='ligand'):
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
                  self.dir + ', creating ' + database + ' database')
            if database == 'ligand':
                self._create_csv_lig(path)
            elif database == 'QD':
                self._create_csv_qd(path)
            else:
                raise TypeError()
        return path

    def _create_csv_lig(self, path):
        """ Create a ligand database and and return its absolute path.

        :param str path: The path to the database.
        """
        idx = sorted(['hdf5 index', 'settings1', 'formula'])
        idx = pd.MultiIndex.from_tuples([(i, '') for i in idx], names=['index', 'sub index'])
        columns = pd.MultiIndex.from_tuples([(None, None)], names=['smiles', 'anchor'])
        df = pd.DataFrame(None, index=idx, columns=columns)
        df.to_csv(path)

    def _create_csv_qd(self, path):
        """ Create a QD database and and return its absolute path.

        :param str path: The path to the database.
        """
        idx = sorted(['hdf5 index', 'settings1', 'settings2', 'ligand count'])
        idx = pd.MultiIndex.from_tuples([(i, '') for i in idx], names=['index', 'sub index'])
        columns = pd.MultiIndex.from_tuples(
                [(None, None, None, None)],
                names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
        )
        df = pd.DataFrame(None, index=idx, columns=columns)
        df.to_csv(path)

    def create_hdf5(self, path, name='structures.hdf5'):
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
        databases = 'core', 'ligand', 'QD'
        for name in databases:
            if name not in hdf5:
                hdf5.create_dataset(name=name, data=np.empty((0, 1), dtype='S80'),
                                    chunks=True, maxshape=(None, None), compression='gzip')
        hdf5.close()
        return path

    def create_yaml(self, path, name='job_settings.yaml'):
        """ Create a job settings database (yaml format), if it does not exist, and return
        its absolute path.

        :param str path: The path to the database.
        :param str name: The filename of the database (excluding its path)
        :return: The absolute path to the job settings database.
        :rtype: |str|_.
        """
        path = join(path, name)
        if not isfile(path):
            with open(path, 'w') as file:
                file.write()
        return path
