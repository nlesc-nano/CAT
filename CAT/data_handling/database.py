""" A module which manages all interactions with the database. """

__all__ = ['ligand_from_database', 'ligand_to_database', 'qd_from_database', 'qd_to_database']

import os
from os.path import (join, isfile)

import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit

from ..utils import get_time


def get_ligand_database(arg):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists. """
    file = join(arg.optional.database.dirname, 'Ligand_database.csv')
    df_index = sorted(['anchor', 'pdb', 'formula', 'settings'])

    if isfile(file):  # The database exists
        df = pd.read_csv(file, index_col=0)
        assert df_index == sorted(list(df.index))
    else:
        print(get_time() + 'Ligand_database.csv not found in ' +
              arg.optional.database.dirname + ', creating ligand database')
        df = pd.DataFrame(None, index=df_index)

    return df


def ligand_from_database(ligand_list, arg):
    """
    Open the database and check if a ligands is already present in the database based on its SMILES
    string. Try to pull the strucure(s) if it is.

    ligand_list <list> [<plams.Molecule>]: A list of ligands.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list> [<plams.Molecule>]: A database of previous calculations.
    """
    df = get_ligand_database(arg)
    for i, lig in enumerate(ligand_list):
        smiles = lig.properties.smiles
        if smiles in df:  # The ligand is present in the database
            pdb = join(arg.optional.ligand.dirname, df[smiles]['pdb'])
            try:
                lig_new = molkit.readpdb(pdb, proximityBonding=False)
                lig_new.properties = lig.properties
                lig_new.properties.read = True
                ligand_list[i] = lig_new
            except FileNotFoundError:
                print(get_time() + lig.properties.name + ' was found in the database, ' +
                      'yet its .pdb file is absent from ' + arg.optional.ligand.dirname)

    return ligand_list


def ligand_to_database(ligand_list, arg):
    """ Write ligands to the ligand database. """
    # Check if previous entries can be overwritten
    overwrite = 'ligand' in arg.optional.database.overwrite
    file = join(arg.optional.database.dirname, 'Ligand_database.csv')
    df = get_ligand_database(arg)

    # A loop which **does not** allow previous entries in the database to be overwritten
    if not overwrite:
        for lig in ligand_list:
            smiles = lig.properties.smiles
            if smiles not in df:
                df[smiles] = None
                df[smiles]['anchor'] = lig.properties.group
                df[smiles]['pdb'] = lig.properties.name + '.pdb'
                df[smiles]['formula'] = lig.get_formula()

    # A loop which **does** allow previous entries in the database to be overwritten
    else:
        for lig in ligand_list:
            smiles = lig.properties.smiles
            if smiles not in df:
                df[smiles] = None
            df[smiles]['anchor'] = lig.properties.group
            df[smiles]['pdb'] = lig.properties.name + '.pdb'
            df[smiles]['formula'] = lig.get_formula()
            df[smiles]['settings'] = None

    # Export the database
    df.to_csv(file)


def get_qd_database(arg):
    """ get the database and return it as dataframe.
    Create a new database if no previous database exists. """
    file = join(arg.optional.database.dirname, 'QD_database.csv')
    df_index = sorted(['ligand count', 'ligand anchor', 'ligand',
                       'core', 'pdb', 'formula', 'settings'])

    if isfile(file):  # The database exists
        df = pd.read_csv(file, index_col=0)
        assert df_index == sorted(list(df.index))
    else:
        print(get_time() + 'QD_database.csv not found in ' +
              arg.optional.database.dirname + ', creating quantum dot database')
        df = pd.DataFrame(None, index=df_index)

    return df


def qd_from_database(ligand_list, core_list, arg):
    """
    Open the database and check if a quantum dot is already present in the database based on its
    name string. Try to pull the strucure(s) if it is.

    ligand_list <list> [<plams.Molecule>]: A list of ligands.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list> [<plams.Molecule>]: A database of previous calculations.
    """
    df = get_qd_database(arg)
    qd_list = []

    for core in core_list:
        for ligand in ligand_list:
            i = str(len(core.properties.dummies))
            name = core.properties.name + '__' + i + '_' + ligand.properties.name
            if name in df:
                pdb = join(arg.optional.qd.dirname, df[name]['pdb'])
                try:
                    qd = molkit.readpdb(pdb, proximityBonding=False)
                    qd.properties.indices = []
                    for i, at in enumerate(qd, 1):
                        if at.properties.pdb_info.ResidueName == 'COR':
                            qd.properties.indices.append(i)
                        elif at.properties.charge != 0:
                            qd.properties.indices.append(i)

                    qd.properties.path = arg.optional.qd.dirname
                    qd.properties.core = core.properties.name
                    qd.properties.ligand = ligand.properties.smiles
                    qd.properties.ligand_anchor = ligand.properties.anchor
                    qd.properties.ligand_count = qd[-1].properties.pdb_info.ResidueNumber - 1
                    qd.properties.name = core.properties.name + '__'
                    qd.properties.name += str(qd.properties.ligand_count)
                    qd.properties.name += '_' + ligand.properties.name
                    qd_list.append(qd)
                except FileNotFoundError:
                    print(get_time() + name + ' was found in the database, ' +
                          'yet its .pdb file is absent from ' + arg.optional.qd.dirname)
            else:
                qd_list.append((core, ligand))

    return qd_list


def qd_to_database(qd_list, arg):
    """ Write quantum dots to the quantum dot database. """
    # Check if previous entries can be overwritten
    overwrite = 'qd' in arg.optional.database.overwrite
    file = join(arg.optional.database.dirname, 'QD_database.csv')
    df = get_qd_database(arg)

    # A loop which **does not** allow previous entries in the database to be overwritten
    if not overwrite:
        for qd in qd_list:
            name = qd.properties.name
            if name not in df:
                df[name] = None
                df[name]['ligand anchor'] = qd.properties.ligand_anchor
                df[name]['ligand count'] = qd.properties.ligand_count
                df[name]['ligand'] = qd.properties.ligand
                df[name]['core'] = qd.properties.core
                df[name]['pdb'] = name + '.pdb'
                df[name]['formula'] = qd.get_formula()

    # A loop which **does** allow previous entries in the database to be overwritten
    else:
        for qd in qd_list:
            name = qd.properties.name
            if name not in df:
                df[name] = None
            df[name]['ligand anchor'] = qd.properties.ligand_anchor
            df[name]['ligand count'] = qd.properties.ligand_count
            df[name]['ligand'] = qd.properties.ligand
            df[name]['core'] = qd.properties.core
            df[name]['pdb'] = name + '.pdb'
            df[name]['formula'] = qd.get_formula()
            df[name]['settings'] = None

    # Export the database
    df.to_csv(file)


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
