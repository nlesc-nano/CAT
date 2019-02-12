""" A module handling the interaction with all other modules, functioning as recipe. """

__all__ = ['prep']

import time
from itertools import chain
from os.path import join

from scm.plams.mol.atom import Atom
from scm.plams.core.errors import MoleculeError

from .misc import (check_sys_var, get_time, create_dir)
from .qd_functions import (find_substructure, find_substructure_split)
from .analysis.asa import init_asa
from .analysis.ligand_bde import init_bde
from .analysis.ligand_solvation import init_solv
from .data_handling.database import (read_database, write_database)
from .data_handling.mol_import import read_mol
from .data_handling.sanitize_input import sanitize_arg_dict
from .attachment.qd_opt import init_qd_opt
from .attachment.ligand_opt import optimize_ligand
from .attachment.ligand_attach import ligand_to_qd


def prep(input_ligands, input_cores, path, arg):
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.

    input_ligands <list>[<plams.Molecule>]: A list of all input ligands.
    input_cores <list>[<plams.Molecule>]: A list of all input cores.
    path <str>: The path where all results will be stored.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A list of all quantum dots (core + n*ligands).
    """
    # The start
    time_start = time.time()
    print('\n')

    arg = sanitize_arg_dict(arg)

    # Create the result directories (if they do not exist) and ligand and core lists
    cor_dir, lig_dir, qd_dir = [create_dir(name, path) for name in arg['dir_name_list']]
    ligand_list = read_mol(input_ligands, lig_dir)
    core_list = read_mol(input_cores, cor_dir, is_core=True)

    # Raises an error if mol_list is empty
    if not ligand_list:
        raise IndexError('No valid input ligands were found, aborting run')
    elif not core_list:
        raise IndexError('No valid input cores were found, aborting run')

    # Adds the indices of the core dummy atoms to core.properties.core
    for core in core_list:
        prep_core(core, arg)

    # Optimize the ligands, find functional groups, calculate properties and read/write the results
    ligand_list = prep_ligand_1(ligand_list, path, arg)

    # Combine the core with the ligands, yielding qd, and format the resulting list
    qd_list = list(ligand_to_qd(core, ligand, qd_dir) for core
                   in core_list for ligand in ligand_list)

    # Optimize the quantum dots, perform an activation strain analyses and read/write the results
    qd_list = prep_qd(qd_list, path, arg)

    # The End
    time_end = time.time()
    message = '\n' + get_time()
    message += 'Total elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec'
    print(message)

    return qd_list, core_list, ligand_list


def prep_core(core, arg):
    """
    Function that handles the identification and marking of all core dummy atoms.

    core <plams.Molecule>: The core molecule.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
    dummy = arg['dummy']

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core
    # An additional dummy atom is added at the core center of mass for orientating the ligands
    if not core.properties.dummies:
        core.properties.dummies = [atom for atom in core.atoms if atom.atnum == dummy]
    else:
        core.properties.dummies = [core[index] for index in core.properties.dummies]

    # Delete all core dummy atoms
    for at in reversed(core.properties.dummies):
        core.delete_atom(at)

    # Returns an error if no dummy atoms were found
    if not core.properties.dummies:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')


def prep_ligand_1(ligand_list, path, arg):
    """
    Function that handles ligand operations.
    Read/write the results from/to a ligand database, launch prep_ligand_2() and
        calculate properties with MOPAC + COSMO-RS.

    ligand_list <list>[<plams.Molecule>]: A list of all ligand molecules.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A copy of all ligands for each identified functional group.
    """
    # Open the ligand database and check if the specified ligand(s) is already present
    if arg['use_database']:
        ligand_database = read_database(path, database_name='Ligand_database')
    else:
        ligand_database = None

    # Optimize all ligands and find their functional groups
    ligand_list = list(chain.from_iterable(prep_ligand_2(ligand, ligand_database, arg) for
                                           ligand in ligand_list))
    if not ligand_list:
        raise IndexError('No valid ligand functional groups found, aborting run')

    if arg.ligand_crs:
        check_sys_var()
        for ligand in ligand_list:
            init_solv(ligand, arg.ligand_crs)

    # Write new entries to the ligand database
    if arg['use_database']:
        if not arg['ligand_opt']:
            for ligand in ligand_list:
                ligand.properties.entry = True
        write_database(ligand_list, ligand_database, path, mol_type='ligand')

    return ligand_list


def prep_ligand_2(ligand, database, arg):
    """
    Function that handles ligand operations.
    Optimize the ligand and search for viable user-defined functional groups.

    ligand <plams.Molecule>: The ligand molecule.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A copy of the ligand for each identified functional group.
    """
    split = arg['split']

    # Identify functional groups within the ligand and add a dummy atom to the center of mass.
    if not ligand.properties.dummies:
        ligand_list = find_substructure(ligand, split)
    else:
        if len(ligand.properties.dummies) == 1:
            ligand.properties.dummies = ligand.properties.dummies[0] - 1
            split = False
        elif len(ligand.properties.dummies) == 2:
            ligand.properties.dummies = [i - 1 for i in ligand.properties.dummies]
            split = True
        ligand_list = [find_substructure_split(ligand, ligand.properties.dummies, split)]

    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand_list = [optimize_ligand(ligand, database, arg['ligand_opt']) for
                   ligand in ligand_list if ligand_list]

    return ligand_list


def prep_qd(qd_list, path, arg):
    """
    Function that handles quantum dot (qd, i.e. core + all ligands) operations.
    Optimize the quantum dot, perform and activation strain analyses on the ligands and read/write
        the results from/to a quantom dot database.

    ligand_list <list>[<plams.Molecule>]: A list of all quantom dots.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A list of all optimized quantom dots.
    """
    if not qd_list:
        raise IndexError('No valid quantum dots found, aborting')

    # Open the quantum dot database and check if the specified quantum dot(s) is already present
    if arg['use_database']:
        qd_database = read_database(path, database_name='QD_database')
    else:
        qd_database = None

    # Optimize the qd with the core frozen
    if arg.qd_opt:
        check_sys_var()
        qd_list = list(init_qd_opt(qd, qd_database, arg.qd_opt) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface
    if arg['qd_int']:
        print(get_time() + 'calculating ligand distortion and inter-ligand interaction...')
        qd_list = list(init_asa(qd) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface upon removal of CdX2
    if arg.qd_dissociate:
        # Start the BDE calculation
        print(get_time() + 'calculating ligand dissociation energy...')
        for qd in qd_list:
            qd.properties.energy.BDE = init_bde(qd, arg.qd_dissociate)
            df = qd.properties.energy.BDE
            df.to_excel(join(path, qd.properties.name + '_BDE.xlsx'))

    # Write the new quantum dot results to the quantum dot database
    if arg['use_database']:
        if not arg['qd_opt']:
            for qd in qd_list:
                qd.properties.entry = True
        write_database(qd_list, qd_database, path, mol_type='qd')

    return qd_list
