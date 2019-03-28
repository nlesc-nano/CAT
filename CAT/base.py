""" A module handling the interaction with all other modules, functioning as recipe. """

__all__ = ['prep']

import time

import numpy as np
import pandas as pd

from scm.plams.mol.atom import Atom
from scm.plams.core.errors import MoleculeError

from .utils import (check_sys_var, get_time)

from .analysis.asa import init_asa
from .analysis.ligand_bde import init_bde
from .analysis.ligand_solvation import init_solv

from .data_handling.mol_import import read_mol
from .data_handling.input_sanitizer import (sanitize_path, sanitize_input_mol, sanitize_optional)

from .attachment.qd_opt import init_qd_opt
from .attachment.ligand_opt import init_ligand_opt
from .attachment.ligand_attach import init_qd_construction
from .attachment.ligand_anchoring import init_ligand_anchoring


def prep(arg, return_mol=True):
    """ function that handles all tasks related to prep_core, prep_ligand and prep_qd.

    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_
    :parameter bool return_mol: If qd_df, core_df & ligand_df should be returned or not.
    :return: If ``return=True``, return a dataframe with quantum dots, cores and ligands.
        Molecules are stored in the *mol* column.
    :rtype: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    """
    # The start
    time_start = time.time()
    print('\n')

    # Interpret and extract the input settings
    ligand_df, core_df = prep_input(arg)

    # Adds the indices of the core dummy atoms to core.properties.core
    core_df = prep_core(core_df, arg)

    # Optimize the ligands, find functional groups, calculate properties and read/write the results
    ligand_df = prep_ligand(ligand_df, arg)

    # Combine the cores and ligands; analyze the resulting quantum dots
    qd_df = prep_qd(ligand_df, core_df, arg)

    # The End
    message = get_time() + 'Total elapsed time:\t\t' + '%.4f' % (time.time() - time_start) + ' sec'
    print(message)

    if return_mol:
        return qd_df, core_df, ligand_df


def prep_input(arg):
    """ Interpret and extract the input settings. Returns a list of ligands and a list of cores.

    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_
    :return: A dataframe of ligands and a dataframe of cores.
    :rtype: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    """
    # Interpret arguments
    arg.update(sanitize_path(arg))
    arg.update(sanitize_optional(arg))
    arg.update(sanitize_input_mol(arg))

    # Read the input ligands and cores
    lig_list = read_mol(arg.input_ligands)
    core_list = read_mol(arg.input_cores)
    del arg.input_ligands
    del arg.input_cores

    # Raises an error if lig_list or core_list is empty
    if not lig_list:
        raise MoleculeError('No valid input ligands were found, aborting run')
    elif not core_list:
        raise MoleculeError('No valid input cores were found, aborting run')

    # Store the molecules in dataframes
    columns = pd.MultiIndex.from_tuples([('mol', '')], names=['index', 'sub index'])
    ligand_df = pd.DataFrame(index=np.arange(len(lig_list)), columns=columns)
    ligand_df['mol'] = lig_list
    core_df = pd.DataFrame(index=np.arange(len(core_list)), columns=columns)
    core_df['mol'] = core_list

    return ligand_df, core_df


def prep_core(core_df, arg):
    """ Function that handles the identification and marking of all core dummy atoms.

    :parameter core_df: A dataframe of core molecules. Molecules are stored in the *mol* column.
    :type core_df: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_)
    :return: A dataframe of cores with all dummy/anchor atoms removed.
    :rtype: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
    """
    formula_list = []
    anchor_list = []

    for i, core in enumerate(core_df['mol']):
        # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
        dummy = arg.optional.core.dummy
        formula_list.append(core.get_formula())

        # Returns the indices and Atoms of all dummy atom ligand placeholders in the core
        if core.properties.dummies is None:
            idx, dummies = zip(*[(j, atom) for j, atom in enumerate(core.atoms, 1) if
                                 atom.atnum == dummy])
        else:
            idx, dummies = zip(*[(j, core[i]) for j in core.properties.dummies])
        core.properties.dummies = dummies
        anchor_list.append(''.join([' ' + str(i) for i in sorted(idx)])[1:])

        # Delete all core dummy atoms
        for at in reversed(core.properties.dummies):
            core.delete_atom(at)

        # Returns an error if no dummy atoms were found
        if not core.properties.dummies:
            raise MoleculeError(Atom(atnum=dummy).symbol +
                                ' was specified as dummy atom, yet no dummy atoms were found')

    # Create and return a new dataframe
    idx_tuples = list(zip(formula_list, anchor_list))
    idx = pd.MultiIndex.from_tuples(idx_tuples, names=['formula', 'anchor'])
    ret = core_df.reindex(idx)
    ret['mol'] = core_df['mol'].values
    return ret


def prep_ligand(ligand_df, arg):
    """ Function that handles all ligand operations:
        - Ligand function group identification
        - Ligand geometry optimization
        - Ligand COSMO-RS calculations

    :parameter ligand_df: A dataframe of ligand molecules. Molecules are stored in the *mol* column.
    :type ligand_df: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_)
    """
    # Identify functional groups within the ligand.
    ligand_df = init_ligand_anchoring(ligand_df, arg)

    # Check if any valid functional groups were found
    if not ligand_df['mol'].any():
        raise MoleculeError('No valid functional groups found in any of the ligands, aborting run')

    # Optimize the ligands
    if arg.optional.ligand.optimize:
        init_ligand_opt(ligand_df, arg)

    # Perform a COSMO-RS calculation on the ligands
    if arg.optional.ligand.crs:
        check_sys_var()
        init_solv(ligand_df, arg)

    return ligand_df


def prep_qd(ligand_df, core_df, arg):
    """ Function that handles all quantum dot (qd, i.e. core + all ligands) operations:
        - Constructing the quantum dots
        - Optimizing the quantum dots
        - Peforming activation strain analyses
        - Dissociating ligands on the quantum dot surface

    :parameter ligand_df: A dataframe of ligand molecules. Molecules are stored in the *mol* column.
    :type ligand_df: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    :parameter core_df: A dataframe of core molecules. Molecules are stored in the *mol* column.
    :type core_df: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_)
    :return: A dataframe of quantum dots molecules. Molecules are stored in the *mol* column.
    :rtype: |pd.DataFrame|_ (columns: |str|_, index: |int|_, values: |plams.Molecule|_)
    """
    # Construct the quantum dots
    qd_df = init_qd_construction(ligand_df, core_df, arg)
    if not qd_df['mol'].any():
        raise MoleculeError('No valid quantum dots found, aborting')

    # Optimize the qd with the core frozen
    if arg.optional.qd.optimize:
        check_sys_var()
        init_qd_opt(qd_df, arg)

    # Calculate the interaction between ligands on the quantum dot surface
    if arg.optional.qd.activation_strain:
        print(get_time() + 'calculating ligand distortion and inter-ligand interaction')
        init_asa(qd_df, arg)

    # Calculate the interaction between ligands on the quantum dot surface upon removal of CdX2
    if arg.optional.qd.dissociate:
        # Start the BDE calculation
        print(get_time() + 'calculating ligand dissociation energy')
        init_bde(qd_df, arg)

    return qd_df
