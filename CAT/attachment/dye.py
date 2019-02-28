""" A modules for combining (organic) molecules. """

__all__ = ['bob_ligand', 'bob_core', 'substitution', 'multi_substitution']

from itertools import chain

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams import Molecule, Settings
from scm.plams.tools.geometry import rotation_matrix

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
        return []

    # Construct keyword arguments
    kwarg1, kwarg2 = get_args(core, lig_list, lig_idx)

    # Allign the ligands with the core; perform the ration check
    lig_array = rot_mol_angle(lig_list, lig_vec, core_vec, **kwarg1)
    lig_array, min_dist_array = rot_mol_axis(lig_array, core_vec, **kwarg2)

    # Combine the rotated ligands and core into new molecules
    ret = []
    for lig, xyz, min_dist in zip(lig_list, lig_array, min_dist_array):
        # Copy and manipulate ligand
        lig_cp = lig.copy()
        lig_cp.from_array(xyz)
        lig_cp.properties = Settings()

        # Copy and manipulate core
        core_h = core.properties.the_h
        core_cp = core.copy()
        core_cp.delete_atom(core_cp.closest_atom(core_h))

        # Merge core and ligand
        lig_cp += core_cp
        lig_cp.properties.name = core.properties.name + "_" + lig.properties.name
        lig_cp.properties.min_distance = min_dist
        ret.append(lig_cp)

    return ret


def get_args(core, lig_list, lig_idx):
    # Extract the various arguments from core and ligand_list
    core_other = core.properties.coords_other[0]
    core_other_atom = core.properties.coords_other_atom[0]

    lig_other = [lig[lig.properties.idx_other+1] for lig in lig_list]
    # core_h = core.properties.coords_h[0]

    bond_length = np.array([core_other_atom.radius + lig.radius for lig in lig_other])

    # Roll through the arguments of core
    core.properties.the_h = core.properties.coords_h[0]
    at_h = core.closest_atom(core.properties.the_h)
    core.properties.vec = core.properties.vec[1:]
    core.properties.coords_other = core.properties.coords_other[1:]
    core.properties.coords_h = core.properties.coords_h[1:]

    # Convert core into an array and mask core_h with np.nan
    core_array = core.as_array()
    idx = core.atoms.index(at_h)
    core_array[idx] = np.nan

    # Create a dictionary of arguments
    kwarg1 = {'atoms_other': core_other, 'idx': lig_idx, 'bond_length': bond_length}
    kwarg2 = {'atoms_other': core_array, 'dist_to_self': False, 'idx': lig_idx, 'ret_min_dist': True}

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

    mol.properties.coords_h_atom = at_h
    mol.properties.coords_other_atom = at_other


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


def find_equivalent_atoms(mol, idx=None, idx_substract=0):
    """ Take a molecule, **mol**, and return the indices of all symmetry equivalent atoms.
    The implemented function is based on finding duplicates in the (sorted) distance matrix,
    as symmetry equivalent atoms have identical inter-atomic distances.

    mol <Molecule> or <np.ndarray>: A PLAMS Molecule or a numpy array.
    idx <None>, <int> or <list> [<int>]: An iterable consisting of atomic indices. ALl atoms in
        **mol** will be examined if *None*.
    idx_substract <int>: Substract a constant from all values in **idx**; usefull for
        interconverting between 1-based and 0-based indices.
    return <list>[<tuple>[<int>]]: A list with tuples of symmetry equivalent atomic indices.
    """
    # Convert a PLAMS molecule to an array
    if isinstance(mol, Molecule):
        mol = mol.as_array()

    # If **idx** is *None*, investigate all atoms in **mol**
    if idx is None:
        idx = slice(0, len(mol))
    else:
        idx = np.array(idx, dtype=int) - idx_substract

    # Create a distance matrix, round it to 2 decimals and isolate all unique rows
    dist = np.around(cdist(mol, mol), decimals=2)
    dist.sort(axis=1)
    unique_at = np.unique(dist[idx], axis=0)

    # Find and return the indices of all duplicate rows in the distance matrix
    return [tuple(j for j, ar2 in enumerate(dist) if (ar1 == ar2).all()) for ar1 in unique_at]


def substitution(input_ligands, input_cores, rep=False):
    """
    To every list of cores one type of ligand is added.
    Mono_subs contaions of key = name of molecule, value = (coordinates of new molecule,
        shortest distance between core and ligand after its connection).
    """
    lig_idx = np.array([lig.properties.idx_other for lig in input_ligands])
    lig_vec = np.array([lig.properties.vec for lig in input_ligands])
    lig_dict = {'lig_list': input_ligands, 'lig_idx': lig_idx, 'lig_vec': lig_vec}

    ret = (connect_ligands_to_core(lig_dict, core) for core in input_cores)

    return list(chain.from_iterable(ret))


def multi_substitution(input_ligands, input_cores, n=1):
    """ Attach ligands to cores; repeat this process n times. """
    ret = []
    mol_list = input_cores
    for _ in range(n):
        mol_list = substitution(input_ligands, mol_list)
        ret.append(mol_list)
    return list(chain.from_iterable(ret))


def get_rotmat_axis(rot_range, axis='x'):
    """ Calculate the rotation matrix for rotating a vector along an axis.
    A rotation matrix is constructed for each angle in **rot_range**.

    rot_range <np.ndarray>: An array of rotations in radian. """
    ret = np.zeros((len(rot_range), 3, 3))
    ret[:, 0, 0] = ret[:, 1, 1] = ret[:, 2, 2] = np.ones(len(rot_range))

    if axis.lower() == 'x':  # Rotation around x axis
        ret[:, 1, 1] = ret[:, 2, 2] = np.cos(rot_range)
        ret[:, 2, 1] = ret[:, 1, 2] = np.sin(rot_range)
        ret[:, 2, 1] *= -1
    elif axis.lower() == 'y':  # Rotation around y axis
        ret[:, 0, 0] = ret[:, 2, 2] = np.cos(rot_range)
        ret[:, 2, 0] = ret[:, 0, 2] = np.sin(rot_range)
        ret[:, 2, 0] *= -1
    elif axis.lower() == 'z':  # Rotation around z axis
        ret[:, 0, 0] = ret[:, 1, 1] = np.cos(rot_range)
        ret[:, 0, 1] = ret[:, 1, 0] = np.sin(rot_range)
        ret[:, 0, 1] *= -1

    return ret

def reset_origin(mol, at1):
    """ Reset the origin of a molecule by means of translations and rotations.

    mol <plams.Molecule> or <np.ndarray>: A PLAMS molecule or array of xyz coordinates.
    idx <idx>: The atomic index of the atom that will be alligned to the x-axis.
    return <plams.Molecule>: The rotated and translated PLAMS molecule.
    """
    if isinstance(mol, Molecule):
        mol = mol.as_array()

    # Translate
    mol -= mol.mean(axis=0)

    # Rotate
    vec1 = mol[at1] - np.zeros(3)
    vec2 = np.array([1, 0, 0])
    rotmat = rotation_matrix(vec1, vec2)
    mol = mol@rotmat

    return mol


def get_symmetry(mol, decimals=2):
    """ Return the number of equivalent atoms under a number of symmetry operations.

    mol <plams.Molecule> or <np.ndarray>: A PLAMS molecule or array of xyz coordinates.
    decimals <int>: The error marigin (number of decimals) in determining equivalent atoms.
    return <pd.DataFrame>: A Pandas dataframe with the number of equivalent atoms per axis
    per operation.
    """
    if isinstance(mol, Molecule):
        mol = mol.as_array()

    # Prepare the dataframe
    columns = ['x', 'y', 'z']
    index = ['2pi / ' + str(i) for i in range(1,9)] + ['reflection', 'inversion']
    data = np.empty((len(index), 3))
    df = pd.DataFrame(data=data, index=index, columns=columns)
    df.index.name = 'Rotations, reflections and inversions'
    df.columns.name = 'Identical atoms'

    # Fill the dataframe
    rot_range = (2 * np.pi) / np.arange(1, 9)
    for i, j in enumerate(columns):
        # Rotate
        rotmat = get_rotmat_axis(rot_range, axis=j)
        mol_new = np.empty((10, mol.shape[0], 3))
        mol_new[0:8] = mol@rotmat

        # Mirror
        mol_new[8] = mol.copy()
        mol_new[8][i] *= -1

        # Invert
        mol_new[9] = -1 * mol.copy()

        # Calculate the number of equivalent atoms
        dist_mat = np.array([cdist(i, mol) for i in mol_new]).round(decimals)
        df[j] = np.bincount(np.where(dist_mat == 0)[0])

    return df


"""
mol = None
idx = find_equivalent_atoms(mol)
if idx:
    at1, at2 = idx[0][0], idx[0][1]
    mol = reset_origin(mol, at1)
    df = get_symmetry(mol)
"""
