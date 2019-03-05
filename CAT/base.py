""" A module handling the interaction with all other modules, functioning as recipe. """

__all__ = ['prep']

import time
from os.path import join

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
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.

    arg <dict>: A dictionary containing all (optional) arguments.
    return_mol <bool>: If qd_list, core_list & ligand_list should be returned or not.
    return <list>[<plams.Molecule>]: If return_mol=True, return qd_list, core_list & ligand_list.
    """
    # The start
    time_start = time.time()
    print('\n')

    # Interpret arguments
    arg.update(sanitize_path(arg))
    arg.update(sanitize_optional(arg))
    arg.update(sanitize_input_mol(arg))

    # Read the input ligands and cores
    ligand_list = read_mol(arg.input_ligands)
    core_list = read_mol(arg.input_cores)
    del arg.input_ligands
    del arg.input_cores

    # Raises an error if mol_list is empty
    if not ligand_list:
        raise MoleculeError('No valid input ligands were found, aborting run')
    elif not core_list:
        raise MoleculeError('No valid input cores were found, aborting run')

    # Adds the indices of the core dummy atoms to core.properties.core
    prep_core(core_list, arg)

    # Optimize the ligands, find functional groups, calculate properties and read/write the results
    ligand_list = prep_ligand(ligand_list, arg)

    # Combine the cores and ligands; analyze the resulting quantum dots
    qd_list = prep_qd(ligand_list, core_list, arg)

    # The End
    message = get_time() + 'Total elapsed time:\t\t' + '%.4f' % (time.time() - time_start) + ' sec'
    print(message)

    if return_mol:
        return qd_list, core_list, ligand_list


def prep_core(core_list, arg):
    """
    Function that handles the identification and marking of all core dummy atoms.

    core_list <core>[<plams.Molecule>]: A list of core molecules.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    for core in core_list:
        # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
        dummy = arg.optional.core.dummy

        # Returns the indices (integer) of all dummy atom ligand placeholders in the core
        if core.properties.dummies is None:
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


def prep_ligand(ligand_list, arg):
    """
    Function that handles all ligand operations.

    ligand_list <list>[<plams.Molecule>]: A list of all ligand molecules.
    arg <dict>: A dictionary containing all (optional) arguments.
    return <list>[<plams.Molecule>]: A copy of all ligands for each identified functional group.
    """
    # Identify functional groups within the ligand.
    ligand_list = init_ligand_anchoring(ligand_list, arg)

    # Check if any valid functional groups were found
    if not ligand_list:
        raise MoleculeError('No valid functional groups found in any of the ligands, aborting run')

    # Optimize the ligands
    if arg.optional.ligand.optimize:
        ligand_list = init_ligand_opt(ligand_list, arg)

    # Perform a COSMO-RS calculation on the ligands
    if arg.optional.ligand.crs:
        check_sys_var()
        init_solv(ligand_list, arg)

    return ligand_list


def prep_qd(ligand_list, core_list, arg):
    """
    Function that handles all quantum dot (qd, i.e. core + all ligands) operations.

    ligand_list <list>[<plams.Molecule>]: A list of all ligands.
    core_list <list>[<plams.Molecule>]: A list of all cores.
    arg <dict>: A dictionary containing all (optional) arguments.
    return <list>[<plams.Molecule>]: A list of all optimized quantom dots.
    """
    # Construct the quantum dots
    qd_list = init_qd_construction(ligand_list, core_list, arg)
    if not qd_list:
        raise MoleculeError('No valid quantum dots found, aborting')

    # Optimize the qd with the core frozen
    if arg.optional.qd.optimize:
        check_sys_var()
        init_qd_opt(qd_list, arg)

    # Calculate the interaction between ligands on the quantum dot surface
    if arg.optional.qd.int:
        print(get_time() + 'calculating ligand distortion and inter-ligand interaction')
        qd_list = list(init_asa(qd) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface upon removal of CdX2
    if arg.optional.qd.dissociate:
        # Start the BDE calculation
        print(get_time() + 'calculating ligand dissociation energy')
        init_bde(qd_list, arg)

    return qd_list
