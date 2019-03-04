""" A module which manages all interactions with the database. """

__all__ = ['ligand_from_database', 'ligand_to_database', 'qd_from_database', 'qd_to_database']

import os
from os.path import (join, isfile)

import h5py
import numpy as np
import pandas as pd

from scm.plams import Molecule
from scm.plams.core.functions import add_to_class
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import get_time


@add_to_class(Molecule)
def as_pdb_array(self):
    """ Converts a PLAMS molecule into an array of strings, the array consisting of a
    (partially) de-serialized .pdb file.

    return <np.ndarray>: An array of strings (dtype: S80).
    """
    pdb_list = Chem.MolToPDBBlock(molkit.to_rdmol(self)).splitlines()
    return np.array(pdb_list, dtype='S80')


def from_pdb_array(array):
    """ Converts a an array with a (partially) de-serialized .pdb file into a PLAMS molecule.

    array <np.ndarray>: An with a (partially) de-serialized .pdb file.
    return <Molecule>: A PLAMS molecule.
    """
    pdb_str = ''
    for item in array:
        if item:
            pdb_str += item.decode() + '\n'
    return molkit.from_rdmol(Chem.MolFromPDBBlock(pdb_str, removeHs=False, proximityBonding=False))


def get_ligand_database(arg):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists.

    arg <dict>: A dictionary containing all (optional) arguments.
    return <pd.DataFrame>: The ligand database as dataframe.
    """
    df_file = join(arg.optional.database.dirname, 'Ligand_database.csv')
    df_index = sorted(['anchor', 'hdf5 index', 'formula', 'settings'])

    # Check if the database exists and has the proper keys; create it if it does not
    if isfile(df_file):
        df = pd.read_csv(df_file, index_col=0)
        assert df_index == sorted(list(df.index))
    else:
        print(get_time() + 'Ligand_database.csv not found in ' +
              arg.optional.database.dirname + ', creating ligand database')
        df = pd.DataFrame(None, index=df_index)

    # Check if the ligand dataset is already available in Structures.hdf5
    hdf5_file = join(arg.optional.database.dirname, 'structures.hdf5')
    hdf5 = h5py.File(hdf5_file, 'a')
    if 'ligand' not in hdf5:
        hdf5.create_dataset(name='ligand', data=np.empty((0, 1), dtype='S80'),
                            chunks=True, maxshape=(None, None), compression='gzip')

    hdf5.close()
    return df


def _anchor_to_idx(string):
    for i, _ in enumerate(string, 1):
        try:
            return int(string[i:])
        except ValueError:
            pass


def ligand_from_database(ligand_list, arg):
    """
    Open the database and check if a ligands is already present in the database based on its SMILES
    string. Try to pull the strucure(s) if it is.

    ligand_list <list> [<plams.Molecule>]: A list of ligands.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list> [<plams.Molecule>]: A database of previous calculations.
    """
    # Open the database the grab the .hdf5 file with all structures
    df = get_ligand_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')

    # Pull structures from the database
    for i, lig in enumerate(ligand_list):
        smiles = lig.properties.smiles
        if smiles in df:  # The ligand is present in the database
            idx = int(df[smiles]['hdf5 index'])
            lig_new = from_pdb_array(hdf5['ligand'][idx])
            lig_new.properties = lig.properties
            lig_new.properties.read = True
            lig_new.properties.dummies = lig_new[_anchor_to_idx(df[smiles]['anchor'])]
            ligand_list[i] = lig_new

    hdf5.close()
    return ligand_list


def ligand_to_database(ligand_list, arg):
    """ Write ligands to the ligand database.

    ligand_list <list> [<plams.Molecule>]: A list of ligands which are (potentially) to be written
        to the database..
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    print(get_time() + 'Updating Ligand_database.csv')

    # Check if previous entries can be overwritten
    overwrite = 'ligand' in arg.optional.database.overwrite

    # A loop which **does not** allow previous entries in the database to be overwritten
    if not overwrite:
        _ligand_to_data(ligand_list, arg)

    # A loop which **does** allow previous entries in the database to be overwritten
    else:
        _ligand_to_data_overwrite(ligand_list, arg)

    # Optional: export molecules to filetypes other than .hdf5
    if arg.optional.database.mol_format:
        for lig in ligand_list:
            path = join(arg.optional.ligand.dirname, lig.properties.name)
            if 'pdb' in arg.optional.database.mol_format:
                molkit.writepdb(lig, path + '.pdb')
            if 'xyz' in arg.optional.database.mol_format:
                lig.write(path + '.xyz')

    print(get_time() + 'Ligand_database.csv has been updated\n')


def _ligand_to_data_overwrite(ligand_list, arg):
    """ Export ligands to the database; overwriting previous entries if necessary. """
    df = get_ligand_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')

    for lig in ligand_list:
        # Resize the hdf5 dataset
        smiles = lig.properties.smiles
        pdb_array = lig.as_pdb_array()
        if smiles not in df:
            df[smiles] = None
            idx = len(hdf5['ligand'])
            shape = (hdf5['ligand'].shape[0] + 1,
                     max(hdf5['ligand'].shape[1], pdb_array.shape[0]))
            hdf5['ligand'].resize(shape)
        else:
            idx = int(df[smiles]['hdf5 index'])

        # Pad the pdb_array and export to the hdf5 dataset
        if pdb_array.shape != hdf5['ligand'][idx].shape:
            shape = hdf5['ligand'][idx].shape[0] - pdb_array.shape[0]
            pdb_array = np.append(pdb_array, np.zeros(shape, dtype='S80'))
        hdf5['ligand'][idx] = pdb_array

        # Update the database
        df[smiles]['anchor'] = lig.properties.anchor
        df[smiles]['hdf5 index'] = idx
        df[smiles]['formula'] = lig.get_formula()
        df[smiles]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'ligand_database.csv')
    df.to_csv(file)
    hdf5.close()


def _ligand_to_data(ligand_list, arg):
    """ Export ligands to the database without overwriting previous entries. """
    df = get_ligand_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')

    for lig in ligand_list:
        smiles = lig.properties.smiles
        if smiles not in df:
            # Resize the hdf5 dataset
            idx = len(hdf5['ligand'])
            pdb_array = lig.as_pdb_array()
            shape = (hdf5['ligand'].shape[0] + 1,
                     max(hdf5['ligand'].shape[1], pdb_array.shape[0]))
            hdf5['ligand'].shape = shape

            # Pad the pdb_array and export to the hdf5 dataset
            if pdb_array.shape != hdf5['ligand'][idx].shape:
                shape = hdf5['ligand'][idx].shape[0] - pdb_array.shape[0]
                pdb_array = np.append(pdb_array, np.zeros(shape, dtype='S80'))
            hdf5['ligand'][idx] = pdb_array

            # Update the database
            df[smiles] = None
            df[smiles]['anchor'] = lig.properties.anchor
            df[smiles]['hdf5 index'] = idx
            df[smiles]['formula'] = lig.get_formula()
            df[smiles]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'ligand_database.csv')
    df.to_csv(file)
    hdf5.close()


def get_qd_database(arg):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists.

    arg <dict>: A dictionary containing all (optional) arguments.
    return <pd.DataFrame>: The quantum dot database as dataframe.
    """
    df_file = join(arg.optional.database.dirname, 'QD_database.csv')
    df_index = sorted(['ligand count', 'ligand anchor', 'ligand',
                       'core', 'hdf5 index', 'formula', 'settings'])

    # Check if the database exists and has the proper keys; create it if it does not
    if isfile(df_file):  # The database exists
        df = pd.read_csv(df_file, index_col=0)
        assert df_index == sorted(list(df.index))
    else:
        print(get_time() + 'QD_database.csv not found in ' +
              arg.optional.database.dirname + ', creating quantum dot database')
        df = pd.DataFrame(None, index=df_index)

    # Check if the quantum dot dataset is already available in structures.hdf5
    hdf5_file = join(arg.optional.database.dirname, 'structures.hdf5')
    hdf5 = h5py.File(hdf5_file, 'a')
    if 'QD' not in hdf5:
        hdf5.create_dataset(name='QD', data=np.empty((0, 1), dtype='S80'),
                            chunks=True, maxshape=(None, None), compression='gzip')
    hdf5.close()
    return df


def qd_from_database(ligand_list, core_list, arg):
    """
    Open the database and check if a quantum dot is already present in the database based on its
    name string. Try to pull the strucure(s) if it is.

    ligand_list <list> [<plams.Molecule>]: A list of ligands.
    core_list <list> [<plams.Molecule>]: A list of cores.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list> [<tuple> [<plams.Molecule>]]: A list of 2 tuples representing unique core+ligand
        combinations. A tuple is subsituted for a fully fledged quantum dot if a match is found
        in the database.
    """
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'))
    df = get_qd_database(arg)
    qd_list = []

    for core in core_list:
        for ligand in ligand_list:
            i = str(len(core.properties.dummies))
            name = core.properties.name + '__' + i + '_' + ligand.properties.name
            if name in df:
                idx = int(df[name]['hdf5 index'])
                qd = from_pdb_array(hdf5['QD'][idx])

                # Set a property
                qd.properties.indices = []
                for i, at in enumerate(qd, 1):
                    if at.properties.pdb_info.ResidueName == 'COR':
                        qd.properties.indices.append(i)
                    elif at.properties.charge != 0:
                        qd.properties.indices.append(i)

                # Set more properties
                qd.properties.read = True
                qd.properties.path = arg.optional.qd.dirname
                qd.properties.core = core.properties.name
                qd.properties.ligand = ligand.properties.smiles
                qd.properties.ligand_anchor = ligand.properties.anchor
                qd.properties.ligand_count = qd[-1].properties.pdb_info.ResidueNumber - 1
                qd.properties.name = core.properties.name + '__'
                qd.properties.name += str(qd.properties.ligand_count) + '_' + ligand.properties.name
                qd_list.append(qd)
            else:
                qd_list.append((core, ligand))

    hdf5.close()
    return qd_list


def qd_to_database(qd_list, arg):
    """ Write quantum dots to the quantum dot database.

    qd_list <list> [<plams.Molecule>]: A list of quantum dots which are (potentially) to be written
        to the database.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    print(get_time() + 'Updating QD_database.csv')

    # Check if previous entries can be overwritten
    overwrite = 'qd' in arg.optional.database.overwrite

    # A loop which **does not** allow previous entries in the database to be overwritten
    if overwrite:
        _qd_to_data_overwrite(qd_list, arg)

    # A loop which **does** allow previous entries in the database to be overwritten
    else:
        _qd_to_data(qd_list, arg)

    # Optional: export molecules to filetypes other than .hdf5
    if arg.optional.database.mol_format:
        for qd in qd_list:
            path = join(arg.optional.qd.dirname, qd.properties.name)
            if 'pdb' in arg.optional.database.mol_format:
                molkit.writepdb(qd, path + '.pdb')
            if 'xyz' in arg.optional.database.mol_format:
                qd.write(path + '.xyz')


def _qd_to_data_overwrite(qd_list, arg):
    """ Export quantum dots to the database; overwriting previous entries if necessary. """
    df = get_qd_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'))

    for qd in qd_list:
        # Resize the hdf5 dataset
        name = qd.properties.name
        pdb_array = qd.as_pdb_array()
        if name not in df:
            df[name] = None
            idx = len(hdf5['QD'])
            shape = (hdf5['QD'].shape[0] + 1, max(hdf5['QD'].shape[1], pdb_array.shape[0]))
            hdf5['QD'].shape = shape
        else:
            idx = int(df[name]['hdf5 index'])

        # Pad the pdb_array and export to the hdf5 dataset
        if pdb_array.shape != hdf5['QD'][idx].shape:
            shape = hdf5['QD'][idx].shape[0] - pdb_array.shape[0]
            pdb_array = np.append(pdb_array, np.zeros(shape, dtype='S80'))
        hdf5['QD'][idx] = pdb_array

        # Update the database
        df[name]['ligand anchor'] = qd.properties.ligand_anchor
        df[name]['ligand count'] = qd.properties.ligand_count
        df[name]['ligand'] = qd.properties.ligand
        df[name]['core'] = qd.properties.core
        df[name]['hdf5 index'] = idx
        df[name]['formula'] = qd.get_formula()
        df[name]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'QD_database.csv')
    df.to_csv(file)
    hdf5.close()


def _qd_to_data(qd_list, arg):
    """ Export quantum dots to the database without overwriting previous entries if necessary. """
    df = get_qd_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'))

    for qd in qd_list:
        name = qd.properties.name
        if name not in df:
            # Resize the hdf5 dataset
            idx = len(hdf5['QD'])
            pdb_array = qd.as_pdb_array()
            shape = (hdf5['QD'].shape[0] + 1,
                     max(hdf5['QD'].shape[1], pdb_array.shape[0]))
            hdf5['QD'].resize(shape)

            # Pad the pdb_array and export to the hdf5 dataset
            if pdb_array.shape != hdf5['QD'][idx].shape:
                shape = hdf5['QD'][idx].shape[0] - pdb_array.shape[0]
                pdb_array = np.append(pdb_array, np.empty(shape, dtype='S80'))
            hdf5['QD'][idx] = pdb_array

            # Update the database
            df[name] = None
            df[name]['ligand anchor'] = qd.properties.ligand_anchor
            df[name]['ligand count'] = qd.properties.ligand_count
            df[name]['ligand'] = qd.properties.ligand
            df[name]['core'] = qd.properties.core
            df[name]['hdf5 index'] = idx
            df[name]['formula'] = qd.get_formula()
            df[name]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'QD_database.csv')
    df.to_csv(file)
    hdf5.close()


def write_database(mol_list, database, path, mol_type='ligand'):
    """
    Write the new database entries to the database.
    New entries are defined by molecules marked as mol.properties.entry = True.

    mol_list <list>[<plams.Molecule>]: A list of ligands or quantum dots.
    database <pd.DataFrame>: A database of previous calculations.
    path <str>: The path to the database.
    mol_type <str>: 'ligand' for ligands and 'qd' for quantum dots.

    return: A .json and .xlsx file.
    """
    database_entries = []
    for mol in mol_list:
        if mol.properties.entry:
            prop = mol.properties
            if mol_type == 'ligand':
                database_entries.append(
                    [prop.name,
                     prop.group,
                     mol.get_formula().split('Xx')[0],
                     os.path.join(prop.path, prop.name.split('@')[0]) + '.pdb',
                     os.path.join(prop.path, prop.name.split('@')[0]) + '.opt.pdb',
                     prop.smiles,
                     prop.energy.E_solv.Acetone,
                     prop.energy.E_solv.Acetonitrile,
                     prop.energy.E_solv.DMF,
                     prop.energy.E_solv.DMSO,
                     prop.energy.E_solv.EtOAc,
                     prop.energy.E_solv.Ethanol,
                     prop.energy.E_solv.Hexane,
                     prop.energy.E_solv.Toluene,
                     prop.energy.E_solv.Water,
                     prop.gamma_solv.Acetone,
                     prop.gamma_solv.Acetonitrile,
                     prop.gamma_solv.DMF,
                     prop.gamma_solv.DMSO,
                     prop.gamma_solv.EtOAc,
                     prop.gamma_solv.Ethanol,
                     prop.gamma_solv.Hexane,
                     prop.gamma_solv.Toluene,
                     prop.gamma_solv.Water])
            elif mol_type == 'qd':
                database_entries.append(
                    [prop.name,
                     mol.get_formula().split('Xx')[0],
                     os.path.join(prop.path, prop.name) + '.pdb',
                     os.path.join(prop.path, prop.name) + '.opt.pdb',
                     prop.energy.E,
                     prop.energy.Eint,
                     prop.energy.Estrain])

    if database_entries:
        database_entries = list(zip(*database_entries))
        if mol_type == 'ligand':
            keys = ('Ligand_name', 'Ligand_group', 'Ligand_formula', 'Ligand_pdb', 'Ligand_opt_pdb',
                    'Ligand_SMILES', 'Gsolv_Acetone', 'Gsolv_Acetonitrile', 'Gsolv_DMF',
                    'Gsolv_DMSO', 'Gsolv_EtOAc', 'Gsolv_Ethanol', 'Gsolv_Hexane', 'Gsolv_Toluene',
                    'Gsolv_Water',
                    'Gamma_Acetone', 'Gamma_Acetonitrile', 'Gamma_DMF', 'Gamma_DMSO',
                    'Gamma_EtOAc', 'Gamma_Ethanol', 'Gamma_Hexane', 'Gamma_Toluene', 'Gamma_Water')
            name = 'Ligand_database'
            sheet_name = 'Ligand'
        if mol_type == 'qd':
            keys = ('Quantum_dot_name', 'Quantum_dot_formula', 'Quantum_dot_pdb',
                    'Quantum_dot_opt_pdb', 'Quantum_dot_E', 'Quantum_dot_Eint',
                    'Quantum_dot_Estrain')
            name = 'QD_database'
            sheet_name = 'Quantum_dot'

        database_entries = pd.DataFrame(dict(zip(keys, database_entries)))

        if not database.empty:
            database = database.append(database_entries, ignore_index=True)
        else:
            database = database_entries

        path = os.path.join(path, name)
        database.to_excel(path + '.xlsx', sheet_name=sheet_name)
        database.to_json(path + '.json')
