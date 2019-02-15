""" A modules for combining (organic) molecules. """

__all__ = ['bob', 'monosubstitution']

import numpy as np

from scm.plams import Molecule

from .ligand_attach import rot_mol_angle, rot_mol_axis


def prep_for_rot(ligand_list, core):
    """ Attaches multiple ligands to multiple copies of a single core.
    Returns a list of cores with attached ligands, each with the properties.min_distance attribute
    containing the smallest distance between ligand and core.

    ligand_list: An iterable container consisting of PLAMS molecules, each with the properties.lig_h & .lig_other attributes (PLAMS Atoms)
    core: A PLAMS molecule with the properties.core_h & .core_other attributes (PLAMS Atoms)
    """
    # Extract various atoms from the ligand and core properties attribute
    lig_h = [lig.properties.h for lig in ligand_list]
    lig_other = [lig.properties.other for lig in ligand_list]
    core_h = core.properties.h
    core_other = core.properties.other

    # Define ligand and core vector
    lig_vec = Molecule().as_array(atom_subset=lig_other) - Molecule().as_array(atom_subset=lig_h)
    core_vec = np.array(core_h.vector_to(core_other))

    # Create an array with the indices of lig_other
    idx_array = np.array([lig.atoms.index(at) for lig, at in zip(ligand_list, lig_other)])

    # Create an array with bond lengths of the to be formed bonds
    bond_dict = {'CC': 1.54, 'CN': 1.469}
    bond_length = np.array([get_bond_lengths(lig_at, core_other, bond_dict) for lig_at in lig_other])

    kwarg1 = {'atoms_other': core_other, 'idx': idx_array, 'bond_length': bond_length}
    kwarg2 = {'atoms_other': core, 'dist_to_self': False, 'idx': idx_array, 'ret_min_dist': True}

    # Allign the lig_vec & core_vec; perform the ration check
    lig_array = rot_mol_angle(ligand_list, lig_vec, core_vec, **kwarg1)
    lig_array, min_dist_array = rot_mol_axis(lig_array, core_vec, **kwarg2)

    ret = []
    for lig, xyz, min_dist in zip(ligand_list, lig_array, min_dist_array):
        lig_cp = lig.copy()
        lig_cp.from_array(xyz)
        lig_cp += core.copy()

        lig_cp.properties.name = core.properties.name + "_" + lig.properties.name
        lig_cp.properties.min_distance = min_dist
        ret.append(lig_cp)

    return ret


def get_bond_lengths(at1, at2, length_dict, length=1.50):
    """ Take two PLAMS atoms, **at1** and **at2**, and grab a matching bond length from length_dict.
    If no bond lengths are available for this particaular atom pair, set the bond length to
    **length**. """
    # Define the keys
    at1_at2 = at1.symbol + at2.symbol
    at2_at1 = at2.symbol + at1.symbol

    # Search **length_dict** for a matching bond length
    if length_dict.get(at1_at2):
        return length_dict[at1_at2]
    elif length_dict.get(at2_at1):
        return length_dict[at2_at1]
    else:
        return length


def bob(mol):
    """
    Marks a PLAMS molecule with the .properties.h & .properties.other attributes.
    mol <plams.Molecule>: An input molecule with the plams_mol.properties.comment attribute.
    """
    comment = mol.properties.comment
    comment = comment.split()
    del mol.properties.comment
    try:
        mol.properties.h = int(comment[0])
        mol.properties.other = int(comment[1])
    except (IndexError, ValueError) as ex:
        print(ex)
        pass


def monosubstitution(input_ligands, input_cores):
    """
    To every list of cores one type of ligand is added.
    Mono_subs contaions of key = name of molecule, value = (coordinates of new molecule,
        shortest distance between core and ligand after its connection).
    """
    return [prep_for_rot(input_ligands, core) for core in input_cores]


"""
def mono_di_substitution(input_ligands, input_cores, dist_limit):

    mono_subs = monosubstitution(input_ligands, input_cores, dist_limit)

    # Takes monosubstituted core and adds ligands, without duplicas
    di_subs = []
    for d in range(len(input_ligands)):
        ligand = input_ligands[d]
        for c in range(d*len(input_cores), len(mono_subs)):
            core = mono_subs[c][0]
            di_subs.append(connect_two_molecules(core, ligand, dist_limit))

    return mono_subs, di_subs
"""

