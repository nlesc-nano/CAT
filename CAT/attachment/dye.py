""" A modules for combining (organic) molecules. """

__all__ = ['bob_ligand', 'bob_core', 'substitution', 'multi_substitution']

from itertools import chain

import numpy as np

from scipy.spatial.distance import cdist

from rdkit.Chem import AllChem

from scm.plams import Molecule, Settings

from scm.plams.interfaces.molecule.rdkit import to_rdmol

from .ligand_attach import rot_mol


def uff_constrained_opt(mol, constrain=[]):
    """ Perform a constrained UFF optimization on a PLAMS molecule.
        
        :parameter mol: A PLAMS molecule.
        :parameter constrain: A list of indices of to-be frozen atoms.
        """
    # Chem.SanitizeMol(rdkit_mol)
    rdkit_mol = to_rdmol(mol)
    ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, ignoreInterfragInteractions=False)
    for f in constrain:
        ff.AddFixedPoint(f)
    ff.Minimize()
    mol.from_rdmol(rdkit_mol)
    return mol


def connect_ligands_to_core(lig_dict, core, user_min_dist):
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
    lig_IDs = lig_dict['lig_ID']

    # Return if core.properties.vec has been exhausted
    try:
        core_vec = core.properties.vec[0]
    except IndexError:
        return []

    # Construct keyword arguments
    kwarg = get_args(core, lig_list, lig_idx)

    # Allign the ligands with the core; perform the ration check
    lig_array, min_dist_array = rot_mol(lig_list, lig_vec, core_vec, **kwarg)

    # Reading list of ligands that are already attached to the core
    old_ligID = core.properties.ligID

    # Combine the rotated ligands and core into new molecules
    ret = []
    for lig, xyz, min_dist, ligID in zip(lig_list, lig_array, min_dist_array, lig_IDs):
        # Copy and manipulate ligand
        lig_cp = lig.copy()
        lig_cp.guess_bonds()
        lig_cp.from_array(xyz)
        lig_cp.properties = Settings()

        # Adding new ligand to the list
        new_ligID = str(old_ligID) + str(ligID)

        # Copy and manipulate core
        core_h = core.properties.the_h
        core_cp = core.copy()
        core_cp.guess_bonds()
        core_cp.delete_atom(core_cp.closest_atom(core_h))
        core_cp.properties.ligID = new_ligID

        # Merge core and ligand
        lig_cp += core_cp
        lig_cp.properties.name = core.properties.name + "_" + lig.properties.name
        lig_cp.properties.min_distance = min_dist
        
        # If the minimal distance cirteria is not fulfilled ligands possitions are optimised with UFF
        if user_min_dist > min_dist:
           # UFF for molecule with frozen core, only ligands coordinates are effected
            the_h = lig_cp.closest_atom(lig_cp.properties.the_h)        
            core_other = lig_cp.closest_atom(lig_cp.properties.coords_other_arrays[len(new_ligID)-1])
            lig_cp.add_bond(the_h, core_other)
            
            h_gone = len(lig_cp.properties.coords_h_atom) - len(lig_cp.properties.coords_h)
            frozen_ind_rdkit = list(range(len(lig_cp)-core.properties.core_len + h_gone, len(lig_cp))) # list rdkit standards, counts from 0
            lig_cp = uff_constrained_opt(lig_cp, constrain=frozen_ind_rdkit[:-1])
            
            # New distance between core and ligands
            frozen_ind_plams = (np.array(frozen_ind_rdkit)+1).tolist()
            frozen_atoms = np.array([lig_cp[c].coords for c in frozen_ind_plams])            
            only_ligands = np.array([lig_cp[l].coords for l in range(1,len(lig_cp)+1) if l not in frozen_ind_plams])
            min_dist = np.nanmin(cdist(only_ligands, frozen_atoms))
            lig_cp.properties.min_distance = min_dist
        
        ret.append(lig_cp)


    return ret


def get_args(core, lig_list, lig_idx):
    # Extract the various arguments from core and ligand_list
    core_other_atom = core.properties.coords_other_atom[0]
    lig_other = [lig[lig.properties.idx_other+1] for lig in lig_list]
    bond_length = np.array([core_other_atom.radius + lig.radius for lig in lig_other])

    # Roll through the arguments of core
    core.properties.the_h = core.properties.coords_h[0]
    at_h = core.closest_atom(core.properties.the_h)
    core.properties.vec = core.properties.vec[1:]
    core.properties.coords_other_atom = core.properties.coords_other_atom[1:]
    core.properties.coords_other = core.properties.coords_other[1:]
    core.properties.coords_h = core.properties.coords_h[1:]

    # Convert core into an array and mask core_h with np.nan
    core_array = core.as_array()
    idx = core.atoms.index(at_h)
    core_array[idx] = np.nan
    
    return {'atoms_core': core_other_atom, 'atoms_other': core_array,
            'bond_length': bond_length, 'dist_to_self': False, 'idx': lig_idx, 'ret_min_dist': True}


def bob_core(mol):
    """
    Reads molecule's comment and marks a PLAMS molecule with the .properties.h & .properties.other attributes.
    mol <plams.Molecule>: An input molecule with the plams_mol.properties.comment attribute.
    """
    # Read the comment in the second line of the xyz file
    comment = mol.properties.comment
    comment = comment.split()
    mol.properties.core_len = len(mol)

    idx = np.array(comment, dtype=int)
    at_h = [mol[int(i)] for i in idx]
    at_other = [at.bonds[0].other_end(at) for at in at_h]

    dummy = Molecule()
    mol.properties.coords_h = dummy.as_array(atom_subset=at_h)
    mol.properties.coords_other = dummy.as_array(atom_subset=at_other)
    mol.properties.vec = mol.properties.coords_other - mol.properties.coords_h

    mol.properties.coords_h_atom = at_h
    mol.properties.coords_other_atom = at_other
    mol.properties.coords_other_arrays = [np.array(i.coords) for i in at_other]
    mol.guess_bonds()


def bob_ligand(mols):
    """
    Reads molecule's comment and marks a PLAMS molecule with the .properties.h & .properties.other attributes.
    mol <plams.Molecule>: An input molecule with the plams_mol.properties.comment attribute.
    """
    ligID = list(range(len(mols)))
    
    for mol, id in zip(mols, ligID):
        mol.properties.ligID = id
        mol.guess_bonds()


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
        mol.guess_bonds()


def substitution(input_ligands, input_cores,min_dist, rep=False):
    """
    To every list of cores one type of ligand is added.
    Mono_subs contaions of key = name of molecule, value = (coordinates of new molecule,
        shortest distance between core and ligand after its connection).
    """
    lig_idx = np.array([lig.properties.idx_other for lig in input_ligands])
    lig_vec = np.array([lig.properties.vec for lig in input_ligands])
    lig_ID = [lig.properties.ligID for lig in input_ligands]
    
    lig_dict = {'lig_list': input_ligands, 'lig_idx': lig_idx, 'lig_vec': lig_vec, 'lig_ID' : lig_ID}
    ret = (connect_ligands_to_core(lig_dict, core, min_dist) for core in input_cores)
    return list(chain.from_iterable(ret))


def multi_substitution(input_ligands, input_cores, n=1):
    """ Attach ligands to cores; repeat this process n times. """
    ret = []
    mol_list = input_cores
    for _ in range(n):
        mol_list = substitution(input_ligands, mol_list)
        ret.append(mol_list)
    return list(chain.from_iterable(ret))


