"""A module handling the interaction with all other modules, functioning as recipe."""

import time
from typing import (Optional, Tuple)

import pandas as pd

from scm.plams import (Atom, Settings)
from scm.plams.core.errors import MoleculeError

from .properties_dataframe import PropertiesDataFrame

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

__all__ = ['prep']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')


def prep(arg: Settings,
         return_mol: bool = True) -> Optional[Tuple[PropertiesDataFrame]]:
    """Function that handles all tasks related to the three prep functions.

    * :func:`.prep_core`
    * :func:`.prep_ligand`
    * :func:`.prep_qd`

    Parameters
    ----------
    arg : |plams.Settings|_
        A settings object containing all (optional) arguments.

    return_mol : bool
        If qd_df, core_df & ligand_df should be returned or not.

    Returns
    -------
    |CAT.PropertiesDataFrame|_
        Optional: If ``return_mol=True`` return the three QD, core and ligand dataframes.

    """
    # The start
    time_start = time.time()
    print('\n')

    # Interpret and extract the input settings
    ligand_df, core_df = prep_input(arg)

    # Adds the indices of the core dummy atoms to core.properties.core
    core_df = prep_core(core_df)

    # Optimize the ligands, find functional groups, calculate properties and read/write the results
    ligand_df = prep_ligand(ligand_df)

    # Combine the cores and ligands; analyze the resulting quantum dots
    qd_df = prep_qd(ligand_df, core_df)

    # The End
    message = get_time() + 'Total elapsed time:\t\t' + '%.4f' % (time.time() - time_start) + ' sec'
    print(message)

    if return_mol:
        return qd_df, core_df, ligand_df
    return None


def prep_input(arg: Settings) -> Tuple[PropertiesDataFrame, PropertiesDataFrame]:
    """Interpret and extract the input settings. Returns a list of ligands and a list of cores.

    Parameters
    ----------
    |plams.Settings|_
        A settings object containing all (optional) arguments.

    Returns
    -------
    |tuple|_ [|CAT.PropertiesDataFrame|_]
        A tuple containing the ligand and core dataframe.

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
    columns = pd.MultiIndex.from_tuples([MOL], names=['index', 'sub index'])

    ligand_df = PropertiesDataFrame(index=pd.RangeIndex(len(lig_list)),
                                    columns=columns.copy(),
                                    properties=arg)
    core_df = PropertiesDataFrame(index=pd.RangeIndex(len(core_list)),
                                  columns=columns.copy(),
                                  properties=arg)

    ligand_df[MOL] = lig_list
    core_df[MOL] = core_list

    return ligand_df, core_df


def prep_core(core_df: PropertiesDataFrame) -> PropertiesDataFrame:
    """Function that handles the identification and marking of all core dummy atoms.

    Parameters
    ----------
    core_df : |CAT.PropertiesDataFrame|_
        A dataframe of core molecules. Molecules are stored in the *mol* column.

    Returns
    -------
    |CAT.PropertiesDataFrame|_
        A dataframe of cores with all dummy/anchor atoms removed.

    """
    # Unpack arguments
    dummy = core_df.properties.optional.core.dummy

    formula_list = []
    anchor_list = []
    for i, core in enumerate(core_df[MOL]):
        # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
        formula_list.append(core.get_formula())

        # Returns the indices and Atoms of all dummy atom ligand placeholders in the core
        if core.properties.dummies is None:
            idx, dummies = zip(*[(j, atom) for j, atom in enumerate(core.atoms, 1) if
                                 atom.atnum == dummy])
        else:
            idx, dummies = zip(*[(j, core[j]) for j in core.properties.dummies])
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
    ret[MOL] = core_df[MOL].values
    return ret


def prep_ligand(ligand_df: PropertiesDataFrame) -> PropertiesDataFrame:
    """Function that handles all ligand operations.

    * Ligand function group identification
    * Ligand geometry optimization
    * Ligand COSMO-RS calculations

    Parameters
    ----------
    ligand_df : |CAT.PropertiesDataFrame|_
        A dataframe of ligand molecules. Molecules are stored in the *mol* column.

    Returns
    -------
    |CAT.PropertiesDataFrame|_
        A new dataframe containing only valid ligands.

    """
    # Unpack arguments
    optimize = ligand_df.properties.optional.ligand.optimize
    crs = ligand_df.properties.optional.ligand.crs

    # Identify functional groups within the ligand.
    ligand_df = init_ligand_anchoring(ligand_df)

    # Check if any valid functional groups were found
    if not ligand_df[MOL].any():
        raise MoleculeError('No valid functional groups found in any of the ligands, aborting run')

    # Optimize the ligands
    if optimize:
        init_ligand_opt(ligand_df)

    # Perform a COSMO-RS calculation on the ligands
    if crs:
        check_sys_var()
        init_solv(ligand_df)

    return ligand_df


def prep_qd(ligand_df: PropertiesDataFrame,
            core_df: PropertiesDataFrame) -> PropertiesDataFrame:
    """Function that handles all quantum dot (qd, i.e. core + all ligands) operations.

    * Constructing the quantum dots
    * Optimizing the quantum dots
    * Peforming activation strain analyses
    * Dissociating ligands on the quantum dot surface

    Parameters
    ----------
    ligand_df : |CAT.PropertiesDataFrame|_
        A dataframe of ligand molecules. Molecules are stored in the *mol* column.

    core_df : |CAT.PropertiesDataFrame|_
        A dataframe of core molecules. Molecules are stored in the *mol* column.

    Returns
    -------
    |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots molecules. Molecules are stored in the *mol* column.

    """
    # Unpack arguments
    optimize = ligand_df.properties.arg.optional.qd.optimize
    dissociate = ligand_df.properties.arg.optional.qd.dissociate
    activation_strain = ligand_df.properties.optional.qd.activation_strain

    # Construct the quantum dots
    qd_df = init_qd_construction(ligand_df, core_df)
    if not qd_df[MOL].any():
        raise MoleculeError('No valid quantum dots found, aborting')

    # Optimize the qd with the core frozen
    if optimize:
        check_sys_var()
        init_qd_opt(qd_df)

    # Calculate the interaction between ligands on the quantum dot surface
    if activation_strain:
        print(get_time() + 'calculating ligand distortion and inter-ligand interaction')
        init_asa(qd_df)

    # Calculate the interaction between ligands on the quantum dot surface upon removal of CdX2
    if dissociate:
        # Start the BDE calculation
        print(get_time() + 'calculating ligand dissociation energy')
        init_bde(qd_df)

    return qd_df
