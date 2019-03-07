""" A module which manages all interactions with the database. """

__all__ = ['mol_from_database', 'mol_to_database', 'property_to_database', 'get_empty_columns']

from os.path import (join, isfile)

import h5py
import numpy as np
import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import get_time


""" ############################  Functions for .pdb conversion  ############################## """


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
        ret[i][:len(item)] = item

    return ret


def from_pdb_array(array):
    """ Converts a an array with a (partially) de-serialized .pdb file into a PLAMS molecule.

    array <np.ndarray>: An with a (partially) de-serialized .pdb file.
    return <Molecule>: A PLAMS molecule.
    """
    pdb_str = ''.join([item.decode() + '\n' for item in array if item])
    return molkit.from_rdmol(Chem.MolFromPDBBlock(pdb_str, removeHs=False, proximityBonding=False))


""" ####################  Functions for database creation and retrieval  ###################### """


def get_database(arg, database='ligand'):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists.

    arg <dict>: A dictionary containing all (optional) arguments.
    database <str>: The type of database, accepted values are 'ligand' and 'qd'.
    return <pd.DataFrame>: The ligand or QD database as dataframe.
    """
    database = _sanitize_database_name(database)
    df_file = join(arg.optional.database.dirname, database + '_database.csv')

    # Check if the database exists and has the proper keys; create it if it does not
    if isfile(df_file):
        header_dict = {'ligand': [0, 1], 'QD': [0, 1, 2, 3]}
        header = header_dict[database]
        df = pd.read_csv(df_file, index_col=[0, 1], header=header, keep_default_na=False)
    else:
        print(get_time() + database + '_database.csv not found in ' +
              arg.optional.database.dirname + ', creating ' + database + ' database')
        data_dict = {'ligand': _export_ligand_df, 'QD': _export_qd_df}
        df = data_dict[database](df_file)

    # Check if the dataset is already available in structures.hdf5
    hdf5_file = join(arg.optional.database.dirname, 'structures.hdf5')
    hdf5 = h5py.File(hdf5_file, 'a')
    if database not in hdf5:
        hdf5.create_dataset(name=database, data=np.empty((0, 1), dtype='S80'),
                            chunks=True, maxshape=(None, None), compression='gzip')
    hdf5.close()
    return df


def _sanitize_database_name(database):
    """ """
    if len(database) == 2:
        ret = database.upper()
    else:
        ret = database.lower()
    assert ret in ('QD', 'ligand')
    return ret


def _export_ligand_df(df_file):
    """ """
    idx = sorted(['hdf5 index', 'settings', 'formula'])
    idx = pd.MultiIndex.from_tuples([(i, '') for i in idx], names=['index', 'sub index'])
    columns = pd.MultiIndex.from_tuples([(None, None)], names=['smiles', 'anchor'])
    df = pd.DataFrame(None, index=idx, columns=columns)
    df.to_csv(df_file)
    return df


def _export_qd_df(df_file):
    """ """
    idx = sorted(['hdf5 index', 'settings', 'ligand count'])
    idx = pd.MultiIndex.from_tuples([(i, '') for i in idx], names=['index', 'sub index'])
    columns = pd.MultiIndex.from_tuples(
            [(None, None, None, None)],
            names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    )
    df = pd.DataFrame(None, index=idx, columns=columns)
    df.to_csv(df_file)
    return df


""" ########################  Functions for pulling from the database  ######################## """


def mol_from_database(mol_list, arg, database='ligand', mol_list2=None):
    """ Open the database and check if a mol is already in it

    mol_list <list> [<plams.Molecule>]: A list of molecules.
    arg <dict>: A dictionary containing all (optional) arguments.
    database <str>: The type of database, accepted values are 'ligand' and 'qd'.

    return <list> [<plams.Molecule>]: A database of previous calculations.
    """
    # Open the database the grab the .hdf5 file with all structures
    database = _sanitize_database_name(database)
    df = get_database(arg, database)
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')

    # Pull structures from the database
    data_dict = {'ligand': _get_ligand_list, 'QD': _get_qd_list}
    mol_list = data_dict[database](hdf5, df, arg, mol_list, mol_list2)

    # Close the .hdf5 file and return
    hdf5.close()
    return mol_list

def _get_ligand_list(hdf5, df, arg, ligand_list, mol_list2=None):
    """ Grab ligands from the ligand database. """
    for i, lig in enumerate(ligand_list):
        key = lig.properties.smiles, lig.properties.anchor
        if key in df:  # The ligand is present in the database
            idx = int(df[key]['hdf5 index'])
            lig_new = from_pdb_array(hdf5['ligand'][idx])
            lig_new.properties = lig.properties
            lig_new.properties.read = True
            lig_new.properties.dummies = lig_new[_anchor_to_idx(key[1])]
            ligand_list[i] = lig_new

    return ligand_list


def _get_qd_list(hdf5, df, arg, ligand_list, core_list):
    """ Grab quantum dots from the QD database. """
    qd_list = []
    for core in core_list:
        for lig in ligand_list:
            key = (core.properties.formula, core.properties.anchor,
                   lig.properties.smiles, lig.properties.anchor)
            if key in df:
                idx = int(df[key]['hdf5 index'])
                qd = from_pdb_array(hdf5['QD'][idx])

                # Set more properties
                qd.properties.indices = _get_qd_indices(qd, lig, key)
                qd.properties.read = True
                qd.properties.path = arg.optional.qd.dirname
                qd.properties.core = key[0]
                qd.properties.core_anchor = key[1]
                qd.properties.ligand = key[2]
                qd.properties.ligand_anchor = key[3]
                qd.properties.ligand_count = qd[-1].properties.pdb_info.ResidueNumber - 1
                qd.properties.name = core.properties.name + '__'
                qd.properties.name += str(qd.properties.ligand_count) + '_' + lig.properties.name
                qd_list.append(qd)
            else:
                qd_list.append((core, lig))

    return qd_list


def _get_qd_indices(qd, ligand, key):
    """ """
    # Collect the indices of all atoms in the core
    ret = []
    for i, at in enumerate(qd, 1):
        if at.properties.pdb_info.ResidueName == 'COR':
            ret.append(i)
        else:
            break

    # Collect the indices of all ligand anchors
    i += _anchor_to_idx(key[3]) - 1
    k = len(ligand)
    for j, _ in enumerate(qd.atoms[i::k]):
        idx = i + j*k
        qd[idx].properties.anchor = True
        ret.append(idx)

    return ret


def _anchor_to_idx(string):
    """ """
    for i, _ in enumerate(string, 1):
        try:
            return int(string[i:])
        except ValueError:
            pass


""" ################  Functions for exporting molecules to the database  ###################### """


def mol_to_database(mol_list, arg, database='ligand'):
    """ Write ligands to the ligand database.

    ligand_list <list> [<plams.Molecule>]: A list of ligands which are (potentially) to be written
        to the database..
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    database = _sanitize_database_name(database)
    print(get_time() + 'Updating ' + database + '_database.csv')

    # A loop which **does** allow previous entries in the database to be overwritten
    if database == 'ligand' and database in arg.optional.database.overwrite:
        _ligand_to_data_overwrite(mol_list, arg)
    elif database == 'ligand' and database not in arg.optional.database.overwrite:
        _ligand_to_data(mol_list, arg)
    elif database == 'QD' and database in arg.optional.database.overwrite:
        _qd_to_data_overwrite(mol_list, arg)
    elif database == 'QD' and database not in arg.optional.database.overwrite:
        _qd_to_data(mol_list, arg)

    # Optional: export molecules to filetypes other than .hdf5
    if arg.optional.database.mol_format:
        export_mol(mol_list, arg, database='ligand')

    print(get_time() + database + '_database.csv has been updated\n')


def export_mol(mol_list, arg, database='ligand'):
    """ """
    if database in arg.optional.database.overwrite:
        for mol in mol_list:
            path = join(arg.optional.ligand.dirname, mol.properties.name)
            if 'pdb' in arg.optional.database.mol_format:
                molkit.writepdb(mol, path + '.pdb')
            if 'xyz' in arg.optional.database.mol_format:
                mol.write(path + '.xyz')
    else:
        for mol in mol_list:
            path = join(arg.optional.ligand.dirname, mol.properties.name)
            if 'pdb' in arg.optional.database.mol_format and not isfile(path + '.pdb'):
                molkit.writepdb(mol, path + '.pdb')
            if 'xyz' in arg.optional.database.mol_format and not isfile(path + '.xyz'):
                mol.write(path + '.xyz')


def _ligand_to_data_overwrite(ligand_list, arg):
    """ Export ligands to the database; overwriting previous entries if necessary. """
    # Open the database
    df = get_database(arg, 'ligand')
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')
    j = hdf5['ligand'].shape[0]

    # Split the ligand_list into a list of new ligands and a list of to be overridden ligands
    lig_new = []
    lig_old = []
    idx_old = []
    for lig in ligand_list:
        key = lig.properties.smiles, lig.properties.anchor
        if key in df:
            lig_old.append(lig)
            idx_old.append(int(df[key]['hdf5 index']))
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
        idx_old, lig_old = zip(*[(i, item) for i, item in sorted(zip(idx_old, lig_old))])
        pdb_old = as_pdb_array(lig_old, min_size=hdf5['ligand'].shape[1])
        hdf5['ligand'][idx_old] = pdb_old
        for lig in lig_old:
            key = lig.properties.smiles, lig.properties.anchor
            df[key]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'ligand_database.csv')
    df.to_csv(file)
    hdf5.close()


def _ligand_to_data(ligand_list, arg):
    """ Export ligands to the database without overwriting previous entries. """
    # Open the database
    df = get_database(arg, 'ligand')
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')

    # Remove ligand entries from ligand_list if they are already present in the database
    ligand_list = [lig for lig in ligand_list if
                   (lig.properties.smiles, lig.properties.anchor) not in df]

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
    df.to_csv(file)
    hdf5.close()


def _qd_to_data_overwrite(qd_list, arg):
    """ Export quantum dots to the database; overwriting previous entries if necessary. """
    # Open the database
    df = get_database(arg, 'QD')
    hdf5 = h5py.File(join(arg.optional.database.dirname, 'structures.hdf5'), 'a')
    j = hdf5['QD'].shape[0]

    # Split the qd_list into a list of new QDs and a list of to be overridden QDs
    qd_new = []
    qd_old = []
    idx_old = []
    for qd in qd_list:
        key = (qd.properties.core, qd.properties.core_anchor,
               qd.properties.ligand, qd.properties.ligand_anchor)
        if key in df:
            qd_old.append(qd)
            idx_old.append(int(df[key]['hdf5 index']))
        else:
            qd_new.append(qd)

    # Update the database with new quantum dots
    if qd_new:
        pdb_new = as_pdb_array(qd_new, min_size=hdf5['QD'].shape[1])
        hdf5['QD'].shape = j + len(qd_new), pdb_new.shape[1]
        hdf5['QD'][j:j+len(qd_new)] = pdb_new
        for i, qd in enumerate(qd_new, j):
            key = (qd.properties.core, qd.properties.core_anchor,
                   qd.properties.ligand, qd.properties.ligand_anchor)
            df[key] = None
            df[key]['ligand count'] = qd.properties.ligand_count
            df[key]['hdf5 index'] = i
            df[key]['settings'] = None

    # Update the database with old quantum dots
    if qd_old:
        idx_old, qd_old = zip(*[(i, item) for i, item in sorted(zip(idx_old, qd_old))])
        pdb_old = as_pdb_array(qd_old, min_size=hdf5['QD'].shape[1])
        hdf5['QD'][idx_old] = pdb_old
        for qd in qd_old:
            key = (qd.properties.core, qd.properties.core_anchor,
                   qd.properties.ligand, qd.properties.ligand_anchor)
            df[key]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'Qd_database.csv')
    df.to_csv(file)
    hdf5.close()


def _qd_to_data(qd_list, arg):
    """ Export quantum dots to the database without overwriting previous entries if necessary. """
    # Open the database
    df = get_database(arg, 'QD')
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
        key = (qd.properties.core, qd.properties.core_anchor,
               qd.properties.ligand, qd.properties.ligand_anchor)
        df[key] = None
        df[key]['ligand count'] = qd.properties.ligand_count
        df[key]['hdf5 index'] = i
        df[key]['settings'] = None

    # Export the database
    file = join(arg.optional.database.dirname, 'Qd_database.csv')
    df.to_csv(file)
    hdf5.close()


""" ##################  Functions for exporting non-molecules to the database  ################ """


def get_empty_columns(index, arg, database='ligand'):
    """ Return all columns in **database** where **index** is completely empty.

    index <str>: An indice in the database
    arg <dict>: arg <dict>: A dictionary containing all (optional) arguments.
    return <list>: A list with the column names.
    """
    df = get_database(arg, database)
    df.replace('', np.nan, inplace=True)
    return df.columns[-df.loc[index].isna().all()].tolist()


def property_to_database(df_new, arg, database='ligand'):
    """ Export a property to the ligand or qd database.

    df_new <pd.DataFrame>: A Pandas dataframe with new results.
    arg <dict>: arg <dict>: A dictionary containing all (optional) arguments.
    database <str>: The type of database, accepted values are 'ligand' and 'qd'.
    """
    # Open the database
    database = _sanitize_database_name(database)
    print('\n' + get_time() + 'Updating ' + database + '_database.csv')
    df = get_database(arg, database)

    # Update the database indices
    for i in df_new.index:
        if i not in df.index:
            df.loc[i, :] = None

    # Update the database values
    if database.lower() in arg.optional.database.overwrite:
        df.update(df_new, overwrite=True)
    else:
        df.update(df_new, overwrite=False)

    # Export the database
    file = join(arg.optional.database.dirname, database + '_database.csv')
    df.to_csv(file)
    print(get_time() + database + '_database.csv has been updated\n')
