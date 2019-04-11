""" Substitution symmetry """

__all__ = ['del_equiv_structures']

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams import Molecule, Settings
from scm.plams.tools.geometry import rotation_matrix




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



def supstitution_symmetry(mol):
    """ Returns atomic symbols of substituted atoms (or first conection of non diatomic ligand)
    	Writes type of substitution symetry at the molecular properties

    	mol <plams.Molecule>: A PLAMS molecule
        """
    dataframe, at_simbols,type_of_symetry = [], [], []

    hatoms = mol.properties.coords_h_atom
    substituted_atoms = [mol.closest_atom(i) for i in hatoms]
    not_hatoms = [at for at in substituted_atoms if at.symbol != 'H']
    
    # Defining C atoms conected to substituents and making Molecule object (cmol) out of them
    catoms = mol.properties.coords_other_arrays
    cmol = Molecule()
    for at in catoms:
        cmol.add_atom(mol.closest_atom(at))
    
    # If substitution is linear, simple combinations without repetition can be applied
    if len(not_hatoms) <= 2:
        if len(not_hatoms) == 1:
            print ("One does not simply ask for subsymmetry of one atom!")
            pass
        elif len(not_hatoms) == 0:
            print ("What the hell is happening?!")
            pass
        else:
            atomic_symbols = [at.symbol for at in not_hatoms]
            mol.properties.subsymmetry = 'linear'
    
    else:
        # Defining tree, four or more substituents and caring about their symbols
        # Getting non zero row indices from data frame - defines symmetry type
        dataframe = get_symmetry(cmol,decimals=2)
        atomic_symbols = [at.symbol for at in not_hatoms]
        type_of_symetry = np.unique(dataframe.to_numpy().nonzero()[0])

        # Assign type of symetry and atomic symbols
        if list(type_of_symetry) == [0, 1, 8, 9]:
            mol.properties.subsymmetry = 'D2h'
        else:
            print ("Well, Jelena made me to recognize only rectangles and this is not rectangle!")

    return atomic_symbols


def get_symmetry(mol, decimals=2):
    """ Returns the number of equivalent atoms under a number of symmetry operations.

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
        try:
            df[j] = np.bincount(np.where(dist_mat == 0)[0])
        except (ValueError):
            print("Something went wrong:Length of values does not match length of index ")
    return df


def del_equiv_structures(mols):
    """ 
	Returnes list of molecules wihout duplicats 	

    mols <plams.Molecule>: A list of PLAMS molecules
        """
    notuniques=[]
   
    for mol in mols:
        atomic_symbols = supstitution_symmetry(mol)
        #lig_positions = np.unique(atomic_simbols, return_counts=True)

        mol.properties.subsymbols = atomic_symbols
        # order has two arrays, repeated symbol and number of repetition
        subsymmetry = mol.properties.subsymmetry 
        symbol_combination = subsymm_combinations(atomic_symbols, subsymmetry)
        #print ("combos", symbol_combination)
        notuniques.append(symbol_combination)

    scos = [sorted(sc) for sc in notuniques]
    u, indices = np.unique(scos, return_index=True, axis=0)
    #print ("unique", u)
    picked_mols = [mols[i] for i in list(indices)]

    return picked_mols

def subsymm_combinations(listy,subsymmetry):

    def swap_neighbours(j):
        # swaping neighbours 1,2,3,4 becomes 2,1,4,3
        j[1::2], j[::2] = j[::2],j[1::2]
        return j
    def swap_pairs(j):
        # swaping pairs 1,2,3,4 becomes 3,4,1,2
        j[:2], j[2:] = j[2:], j[:2]
        return j
    
    # Making a string out of molecule combinations
    if subsymmetry == 'D2h':
        a = ''.join(listy)
        b = ''.join(swap_neighbours(listy))
        c = ''.join(swap_pairs(listy))
        d = ''.join(swap_neighbours(listy))
        final = []
        final = [a,b,c,d]
    if subsymmetry == 'linear':
        a = ''.join(listy)
        b = ''.join(swap_neighbours(listy))
        final = []
        final = [a,b]

    return final







