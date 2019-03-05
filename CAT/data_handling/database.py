""" A module which manages all interactions with the database. """

__all__ = ['ligand_from_database', 'ligand_to_database', 'qd_from_database', 'qd_to_database']

from os.path import (join, isfile)

import h5py
import numpy as np
import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import get_time


def as_pdb_array(mol_list, min_size=0):
    """ Converts a PLAMS molecule into an array of strings, the array consisting of a
    (partially) de-serialized .pdb file.

    mol_list <list> [<plams.Molecule>]: A list of PLAMS molecules.
    min_size <int>: The minimumum length of the pdb_array. The array is padded with empty
        strings if required.
    return <np.ndarray>: An array of strings (dtype: S80).
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
        ret[i][0:len(item)] = item

    return ret


def from_pdb_array(array):
    """ Converts a an array with a (partially) de-serialized .pdb file into a PLAMS molecule.

    array <np.ndarray>: An with a (partially) de-serialized .pdb file.
    return <Molecule>: A PLAMS molecule.
    """
    pdb_str = ''.join([item.decode() + '\n' for item in array if item])
    return molkit.from_rdmol(Chem.MolFromPDBBlock(pdb_str, removeHs=False, proximityBonding=False))


def get_ligand_database(arg):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists.

    arg <dict>: A dictionary containing all (optional) arguments.
    return <pd.DataFrame>: The ligand database as dataframe.
    """
    df_file = join(arg.optional.database.dirname, 'Ligand_database.csv')

    # Check if the database exists and has the proper keys; create it if it does not
    if isfile(df_file):
        df = pd.read_csv(df_file, index_col=0).T
    else:
        print(get_time() + 'Ligand_database.csv not found in ' +
              arg.optional.database.dirname + ', creating ligand database')
        idx = sorted(['anchor', 'hdf5 index', 'formula', 'settings', 'smiles'])
        idx = pd.MultiIndex.from_tuples([(i, None) for i in idx], names=['index', 'sub index'])
        columns = pd.MultiIndex.from_tuples([(None, None)], names=['smiles', 'anchor'])
        df = pd.DataFrame(None, index=idx, columns=columns)
        df.T.to_csv(df_file)

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
        key = lig.properties.smiles, lig.properties.anchor
        if key in df:  # The ligand is present in the database
            idx = int(df[key]['hdf5 index'])
            lig_new = from_pdb_array(hdf5['ligand'][idx])
            lig_new.properties = lig.properties
            lig_new.properties.read = True
            lig_new.properties.dummies = lig_new[_anchor_to_idx(key[1])]
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

    # A loop which **does** allow previous entries in the database to be overwritten
    if 'ligand' in arg.optional.database.overwrite:
        _ligand_to_data_overwrite(ligand_list, arg)

    # A loop which **does not** allow previous entries in the database to be overwritten
    else:
        _ligand_to_data(ligand_list, arg)

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
    # Open the database
    df = get_ligand_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')
    j = hdf5['ligand'].shape[0]

    # Split the ligand_list into a list of new ligands and a list of to be overridden ligands
    lig_new = []
    lig_old = []
    idx_old = []
    for lig in ligand_list:
        if lig.properties.name in df:
            lig_old.append(lig)
            idx_old.append(int(df[lig.properties.name]['hdf5 index']))
        else:
            lig_new.append(lig)

    # Update the database with new ligands
    if lig_new:
        pdb_new = as_pdb_array(lig_new, min_size=hdf5['ligand'].shape[1])
        hdf5['ligand'].shape = j + len(lig_new), pdb_new.shape[1]
        hdf5['ligand'][j:j+len(lig_new)] = pdb_new
        for i, lig in enumerate(lig_new, j):
            key = lig.properties.smiles, lig.properties.anchor
            df[key] = None
            df[key]['hdf5 index'] = i
            df[key]['formula'] = lig.get_formula()
            df[key]['settings'] = None

    # Update the database with old ligands
    if lig_old:
        pdb_old = as_pdb_array(lig_old, min_size=hdf5['ligand'].shape[1])
        hdf5['ligand'][idx_old] = pdb_old
        for lig in lig_old:
            key = lig.properties.smiles, lig.properties.anchor
            df[key]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'ligand_database.csv')
    df.T.to_csv(file)
    hdf5.close()


def _ligand_to_data(ligand_list, arg):
    """ Export ligands to the database without overwriting previous entries. """
    # Open the database
    df = get_ligand_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')

    # Remove ligand entries from ligand_list if they are already present in the database
    ligand_list = [lig for lig in ligand_list if lig.properties.name not in df]

    # Prepare the pdb array and reshape the database
    j = hdf5['ligand'].shape[0]
    pdb_array = as_pdb_array(ligand_list, min_size=hdf5['ligand'].shape[1])
    hdf5['ligand'].shape = j + pdb_array.shape[0], pdb_array.shape[1]
    hdf5['ligand'][j:len(hdf5['ligand'])] = pdb_array

    # Update the database
    for i, lig in enumerate(ligand_list, j):
        key = lig.properties.smiles, lig.properties.anchor
        df[key] = None
        df[key]['hdf5 index'] = i
        df[key]['formula'] = lig.get_formula()
        df[key]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'ligand_database.csv')
    df.T.to_csv(file)
    hdf5.close()


def ligand_solv_to_database(solv_df, arg):
    """ Export ligand solvation energies and activity voefficients to the database. """
    # Open and update database
    df = get_ligand_database(arg)
    if arg.optional.database.overwrite:
        df = pd.concat(solv_df, df, sort=False)
    else:
        df = pd.concat(df, solv_df, sort=False)

    # Close the database
    file = join(arg.optional.database.dirname, 'ligand_database.csv')
    df.T.to_csv(file)


def get_qd_database(arg):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists.

    arg <dict>: A dictionary containing all (optional) arguments.
    return <pd.DataFrame>: The quantum dot database as dataframe.
    """
    df_file = join(arg.optional.database.dirname, 'QD_database.csv')

    # Check if the database exists and has the proper keys; create it if it does not
    if isfile(df_file):  # The database exists
        df = pd.read_csv(df_file, index_col=0).T
    else:
        print(get_time() + 'QD_database.csv not found in ' +
              arg.optional.database.dirname + ', creating quantum dot database')
        idx = sorted(['ligand count', 'hdf5 index', 'formula', 'settings'])
        idx = pd.MultiIndex.from_tuples([(i, None) for i in idx], names=['index', 'sub index'])
        columns = pd.MultiIndex.from_tuples([(None, None, None)],
                                            names=['core', 'ligand smiles', 'ligand anchor'])
        df = pd.DataFrame(None, index=idx, columns=columns)
        df.T.to_csv(df_file)

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
        for lig in ligand_list:
            key = core.properties.name, lig.properties.smiles, lig.properties.anchor
            if key in df:
                idx = int(df[key]['hdf5 index'])
                qd = from_pdb_array(hdf5['QD'][idx])

                # Collect all core indices
                qd.properties.indices = []
                for i, at in enumerate(qd, 1):
                    if at.properties.pdb_info.ResidueName == 'COR':
                        qd.properties.indices.append(i)
                    else:
                        break

                # Collect all ligand anchor indices
                k = _anchor_to_idx(key[2])
                for j, _ in enumerate(qd[i+k::k]):
                    qd.properties.indices.append(j)

                # Set more properties
                qd.properties.read = True
                qd.properties.path = arg.optional.qd.dirname
                qd.properties.core = key[0]
                qd.properties.ligand = key[1]
                qd.properties.ligand_anchor = key[2]
                qd.properties.ligand_count = qd[-1].properties.pdb_info.ResidueNumber - 1
                qd.properties.name = core.properties.name + '__'
                qd.properties.name += str(qd.properties.ligand_count) + '_' + lig.properties.name
                qd_list.append(qd)
            else:
                qd_list.append((core, lig))

    hdf5.close()
    return qd_list


def qd_to_database(qd_list, arg):
    """ Write quantum dots to the quantum dot database.

    qd_list <list> [<plams.Molecule>]: A list of quantum dots which are (potentially) to be written
        to the database.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    print(get_time() + 'Updating QD_database.csv')

    # A loop which **does** allow previous entries in the database to be overwritten
    if 'qd' in arg.optional.database.overwrite:
        _qd_to_data_overwrite(qd_list, arg)

    # A loop which **does not** allow previous entries in the database to be overwritten
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

    print(get_time() + 'QD_database.csv has been updated\n')


def _qd_to_data_overwrite(qd_list, arg):
    """ Export quantum dots to the database; overwriting previous entries if necessary. """
    # Open the database
    df = get_qd_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')
    j = hdf5['QD'].shape[0]

    # Split the qd_list into a list of new QDs and a list of to be overridden QDs
    qd_new = []
    qd_old = []
    idx_old = []
    for qd in qd_list:
        if qd.properties.name in df:
            qd_old.append(qd)
            idx_old.append(int(df[qd.properties.name]['hdf5 index']))
        else:
            qd_new.append(qd)

    # Update the database with new quantum dots
    if qd_new:
        pdb_new = as_pdb_array(qd_new, min_size=hdf5['QD'].shape[1])
        hdf5['QD'].shape = j + len(qd_new), pdb_new.shape[1]
        hdf5['QD'][j:j+len(qd_new)] = pdb_new
        for i, qd in enumerate(qd_new, j):
            key = qd.properties.core, qd.properties.ligand, qd.properties.ligand_anchor
            df[key] = None
            df[key]['ligand count'] = qd.properties.ligand_count
            df[key]['hdf5 index'] = i
            df[key]['formula'] = qd.get_formula()
            df[key]['settings'] = None

    # Update the database with old quantum dots
    if qd_old:
        pdb_old = as_pdb_array(qd_old, min_size=hdf5['QD'].shape[1])
        hdf5['QD'][idx_old] = pdb_old
        for qd in qd_old:
            df[qd.properties.name]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'QD_database.csv')
    df.T.to_csv(file)
    hdf5.close()


def _qd_to_data(qd_list, arg):
    """ Export quantum dots to the database without overwriting previous entries if necessary. """
    # Open the database
    df = get_qd_database(arg)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'))

    # Remove ligand entries from ligand_list if they are already present in the database
    qd_list = [qd for qd in qd_list if qd.properties.name not in df]

    # Prepare the pdb array and reshape the database
    j = hdf5['QD'].shape[0]
    pdb_array = as_pdb_array(qd_list, min_size=hdf5['QD'].shape[1])
    hdf5['QD'].shape = j + pdb_array.shape[0], pdb_array.shape[1]
    hdf5['QD'][j:len(hdf5['QD'])] = pdb_array

    # Update the database
    for i, qd in enumerate(qd_list, j):
        key = qd.properties.core, qd.properties.ligand, qd.properties.ligand_anchor
        df[key] = None
        df[key]['ligand count'] = qd.properties.ligand_count
        df[key]['hdf5 index'] = i
        df[key]['formula'] = qd.get_formula()
        df[key]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'QD_database.csv')
    df.T.to_csv(file)
    hdf5.close()
