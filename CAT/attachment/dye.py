""" A modules for combining (organic) molecules. """

__all__ = ['bob_ligand', 'bob_core', 'substitution', 'multi_substitution']

from itertools import chain

import numpy as np

from scm.plams import Molecule, Settings

from .ligand_attach import rot_mol_angle, rot_mol_axis


def connect_ligands_to_core(lig_dict, core):
    """ Attaches multiple ligands to multiple copies of a single core.
    Returns a list of cores with attached ligands, each with the properties.min_distance attribute
    containing the smallest distance between ligand and core.

    ligand_list: An iterable container consisting of PLAMS molecules, each with the properties.lig_h & .lig_other attributes (PLAMS Atoms)
    core: A PLAMS molecule with the properties.core_h & .core_other attributes (PLAMS Atoms)
    """
    # Unpack the ligand dictionary
    lig_list = lig_dict['lig_list']
    lig_idx = lig_dict['lig_idx']
    lig_vec = lig_dict['lig_vec']

    # Return if core.properties.vec has been exhausted
    try:
        core_vec = core.properties.vec[0]
    except IndexError:
        return

    # Construct keyword arguments
    kwarg1, kwarg2 = get_args(core, lig_idx)

    # Allign the ligands with the core; perform the ration check
    lig_array = rot_mol_angle(lig_list, lig_vec, core_vec, **kwarg1)
    lig_array, min_dist_array = rot_mol_axis(lig_array, core_vec, **kwarg2)

    # Combine the rotated ligands and core into new molecules
    ret = []
    for lig, xyz, min_dist in zip(lig_list, lig_array, min_dist_array):
        lig_cp = lig.copy()
        lig_cp.from_array(xyz)
        lig_cp.properties = Settings()
        lig_cp += core.copy()
        lig_cp.properties.name = core.properties.name + "_" + lig.properties.name
        lig_cp.properties.min_distance = min_dist
        ret.append(lig_cp)

    return ret


def get_args(core, lig_idx):
    # Extract the various arguments from core and ligand_list
    core_other = core.properties.coords_other[0]
    core_h = core.properties.coords_h[0]

    # Roll through the arguments of core; delete core_h
    core.properties.vec = core.properties.vec[1:]
    core.properties.coords_other = core.properties.coords_other[1:]
    core.properties.coords_h = core.properties.coords_h[1:]
    core.delete_atom(core.closest_atom(core_h))

    # Create a dictionary of arguments
    kwarg1 = {'atoms_other': core_other, 'idx': lig_idx, 'bond_length': 1.5}
    kwarg2 = {'atoms_other': core, 'dist_to_self': False, 'idx': lig_idx, 'ret_min_dist': True}

    return kwarg1, kwarg2


def bob_core(mol):
    """
    Marks a PLAMS molecule with the .properties.h & .properties.other attributes.
    mol <plams.Molecule>: An input molecule with the plams_mol.properties.comment attribute.
    """
    # Read the comment in the second line of the xyz file
    comment = mol.properties.comment
    comment = comment.split()

    idx = np.array(comment, dtype=int)
    at_h = [mol[int(i)] for i in idx]
    at_other = [at.bonds[0].other_end(at) for at in at_h]

    dummy = Molecule()
    mol.properties.coords_h = dummy.as_array(atom_subset=at_h)
    mol.properties.coords_other = dummy.as_array(atom_subset=at_other)
    mol.properties.vec = mol.properties.coords_other - mol.properties.coords_h


def bob_ligand(mol):
    """
    Marks a PLAMS molecule with the .properties.h & .properties.other attributes.
    mol <plams.Molecule>: An input molecule with the plams_mol.properties.comment attribute.
    """
    # Read the comment in the second line of the xyz file
    comment = mol.properties.comment

    # Set an attirbute with indices and vectors
    mol.properties.vec = np.empty(3)

    # Fill the vector array, delete the marked atom in the ligand .xyz file
    at = mol[int(comment)]
    at_other = at.bonds[0].other_end(at)
    mol.properties.vec = np.array(at.coords) - np.array(at_other.coords)

    mol.delete_atom(at)
    mol.properties.idx_other = mol.atoms.index(at_other)


def substitution(input_ligands, input_cores, rep=False):
    """
    To every list of cores one type of ligand is added.
    Mono_subs contaions of key = name of molecule, value = (coordinates of new molecule,
        shortest distance between core and ligand after its connection).
    """
    lig_idx = np.array([lig.properties.idx_other for lig in input_ligands])
    lig_vec = np.array([lig.properties.vec for lig in input_ligands])
    lig_dict = {'lig_list': input_ligands, 'lig_idx': lig_idx, 'lig_vec': lig_vec}

    if not rep:
        ret = [connect_ligands_to_core(lig_dict, core) for core in input_cores]
    else:
        ret = (connect_ligands_to_core(lig_dict, core)[i:] for i, core in enumerate(input_cores))
    return list(chain.from_iterable(ret))


def multi_substitution(input_ligands, input_cores, n=1):
    """ Attach ligands to cores; repeat this process n times. """
    ret = []
    mol_list = input_cores
    for _ in range(n):
        mol_list = substitution(input_ligands, mol_list)
        ret.append(mol_list)
    return list(chain.from_iterable(ret))
