""" Substitution symmetry """

__all__ = ['del_equiv_structures']

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from scm.plams import Molecule
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
    dataframe,type_of_symetry = [], []
    ligand_identity = mol.properties.ligID
 
    # Defining C atoms conected to substituents and making Molecule object (cmol) out of them
    catoms = mol.properties.coords_other_arrays
    cmol = Molecule()
    for at in catoms:
        cmol.add_atom(mol.closest_atom(at))

    # If substitution is linear, simple combinations without repetition can be applied
    if len(ligand_identity) <= 2:
        if len(ligand_identity) == 1:
            print ("One does not simply ask for subsymmetry of one atom!")
            pass
        elif len(ligand_identity) == 0:
            print ("What the hell is happening?!")
            pass
        else:
            subsymmetry = 'linear'
    else:
        # Getting non zero row indices from data frame - defines symmetry type
        
        dataframe = get_symmetry(cmol,decimals=2)
        type_of_symetry = np.unique(dataframe.to_numpy().nonzero()[0])
        
        # Assign type of symetry and atomic symbols
        if list(type_of_symetry) == [0, 1, 8, 9]:
            subsymmetry = 'D2h'
        else:
            print ("Well, Jelena made me to recognize only rectangles and this is not rectangle!")

    return subsymmetry


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


def del_equiv_structures(mols, subsymmetry=None):
    """ 
	Returnes list of molecules wihout duplicats 	

    mols <plams.Molecule>: A list of PLAMS molecules
        """
    notunique=[]
   
    for mol in mols:
        if subsymmetry == None:
            subsymmetry = supstitution_symmetry(mol)
        
        ligID = list(mol.properties.ligID)
        
        all_permutations = symm_permutations(subsymmetry, ligID) 
        
        notunique.append(all_permutations)

    scos = [sorted(sc) for sc in notunique]
    u, indices = np.unique(scos, return_index=True, axis=0)
    
    unique_molecules = [mols[i] for i in list(indices)]

    return unique_molecules

def symm_permutations(condition, elements): 
    """ For given list of elements, makes permutations taking in account symmetry condition
    <condition>: string 
    <elements>: list
    """
    def swap_neighbours(j):
        """ swaping neighbours inside a list: 1,2,3,4 becomes 2,1,4,3
        <j>: a list
        """
        j[1::2], j[::2] = j[::2],j[1::2]
        return j
    def swap_pairs(j):
        """ swaping pairs inside a list: 1,2,3,4 becomes 3,4,1,2
        <j>: a list
        """
        j[:2], j[2:] = j[2:], j[:2]
        return j
    def rotate_list(l,n):
        return l[n:] + l[:n]
    def swap_two(j):
        j[0], j[1] = j[1], j[0]
        return j

    # Making list of all permutations
    if condition == 'D2h':
        a = ''.join(elements)
        b = ''.join(swap_neighbours(elements))
        c = ''.join(swap_pairs(elements))
        d = ''.join(swap_neighbours(elements))
        final = []
        final = [a,b,c,d]
    if condition == 'linear':
        a = ''.join(elements)
        b = ''.join(swap_neighbours(elements))
        final = []
        final = [a,b]
    if condition == 'triangle':
        a = ''.join(elements)
        b = ''.join(swap_two(elements))
        #c,d = ''.join(rotate_list(a,1)), ''.join(rotate_list(a,2))
        #e,f = ''.join(rotate_list(b,1)), ''.join(rotate_list(b,2))
        final = []
        final = [a,b]#,c,d,e,f]
        print (final)
    return final







