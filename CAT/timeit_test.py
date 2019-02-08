import itertools
from itertools import chain
import os
import numpy as np
from scipy.spatial.distance import cdist
import time
import math
import timeit
import matplotlib.pyplot as plt
import pandas as pd
import pdb
# import qmflows.components.qd_functions
import copy

from scm.plams import Atom, Molecule, Bond
from scm.plams.core.functions import add_to_class
from scm.plams.core.errors import MoleculeError
from scm.plams.tools.units import Units
from scm.plams.core.settings import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms


@add_to_class(Molecule)
def merge_mol(self, mol_list):
    if isinstance(mol_list, Molecule):
        mol_list = [mol_list]

    for mol in mol_list:
        self.properties.soft_update(mol.properties)

    atom_list = list(itertools.chain.from_iterable(mol.atoms for mol in mol_list))
    bond_list = list(itertools.chain.from_iterable(mol.bonds for mol in mol_list))
    for atom in atom_list:
        atom.mol = self
    for bond in bond_list:
        bond.mol = self
    self.atoms += atom_list
    self.bonds += bond_list


@add_to_class(Molecule)
def from_array(self, array):
    array = array.T
    for at1, x, y, z in zip(self, array[0], array[1], array[2]):
        at1.coords = (x, y, z)


@add_to_class(Molecule)
def to_array(self):
    x, y, z = zip(*[atom.coords for atom in self])
    return np.array((x, y, z)).T


def get_2d(mol):
    def create_rotmat_2d(vec1, vec2):
        a_2d = np.array(vec1) / np.linalg.norm(vec1)
        b_2d = np.array(vec2) / np.linalg.norm(vec2)
        v1, v2, v3 = np.cross(a_2d, b_2d)
        M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])
        return np.identity(3) + M + np.dot(M, M)/(1+np.dot(a_2d, b_2d))

    array = mol.to_array()
    vec1 = array[0] - array[1]
    vec2 = np.array([0, 0, 1])
    rotmat = create_rotmat_2d(vec1, vec2)
    array_out = rotmat.dot(array.T).T
    array_out -= array_out[0] + np.array([0, 0, 1])

    mol_2d = mol.copy()
    mol_2d.from_array(array_out)
    return mol_2d


def get_3d(mol, depth=6):
    def create_rotmat_3d(tens1, tens2):
        """
        Create a 3d rotation matrix of dimensions 3*3*n
        """
        a_3d = (np.array(tens1) / np.linalg.norm(tens1, axis=0)).T
        b_3d = (np.array(tens2) / np.linalg.norm(tens2, axis=0)).T
        v1, v2, v3 = np.cross(a_3d, b_3d).T
        M_3d = np.zeros((3, 3, depth), dtype=np.float64)
        M_3d[0, 1], M_3d[0, 2] = v3, -v2
        M_3d[1, 0], M_3d[1, 2] = -v3, v1
        M_3d[2, 0], M_3d[2, 1] = v2, -v1
        M_3d = M_3d.T
        return np.identity(3) + M_3d + (np.matmul(M_3d, M_3d).T /
                                        (1 + np.einsum('ij,ij->i', a_3d, b_3d))).T

    # Create a m*3*n array representing n copies of mol with m atoms
    array = np.repeat(mol.to_array()[:, :, np.newaxis], depth, axis=2)

    # Rotate n ligands
    tens1 = array[0] - array[1]
    tens2 = np.array(np.random.rand(3, depth))
    rotmat = create_rotmat_3d(tens1, tens2)
    array_out = np.empty((len(mol), 3, depth))
    for i in range(depth):
        array_out[:, :, i] = rotmat[i].dot(array[:, :, i].T).T

    # Translate n ligands
    array_out -= array_out[0] + tens2*3

    # Convert the array into a new output molecule
    mol_3d = mol.copy()
    mol_3d.merge_mol([mol.copy() for i in range(depth - 1)])
    mol_3d.from_array(np.concatenate(array_out.T, axis=1).T)
    return mol_3d


mol = molkit.from_smiles('[O-]CCCCCC')


def function_timeit(mol):
    ret = []
    x = []
    for i in range(6):
        # mol_cp = mol.copy()
        # mol_cp.translate_np((2.0, 2.0, 2.0))
        # mol += mol_cp
        # mol2 = mol.copy()
        mol_list = [mol.copy() for j in range(2**i)]
        x.append(len(list(itertools.chain(mol_list))) + len(mol))
        # array1, array2 = mol.to_array(), mol2.to_array()
        loops = 1000
        functions = [timeit.timeit(lambda: mol.merge_mol(mol_list), number=loops),
                     timeit.timeit(lambda: mol.merge_mol2(mol_list), number=loops)
                     ]
        ret.append([i/loops for i in functions])

        message = '\n\naverage time per loop for ' + mol.get_formula() + ' over '
        message += str(loops) + ' loops:'
        print(message)
        for i, item in enumerate(functions):
            message = 'function ' + str(i).zfill(2) + ':\t'
            message += "{0:.3f}".format(item / loops * 1000**0) + ' s;\t'
            message += "{0:.3f}".format(item / loops * 1000**1) + ' ms;\t'
            message += "{0:.3f}".format(item / loops * 1000**2) + ' us'
            print(message)
        if len(functions) == 2:
            message = '00/01 ratio:\t' + "{0:.4f}".format(functions[0] / functions[1]) + ' / 1.0'
            print(message)
    a, b = zip(*ret)
    return x, a, b


def get_all_dist(mol):
    def to_array2(mol):
        x, y, z = zip(*[atom.coords for atom in atom_gen])
        return np.array((x, y, z)).T

    atom_gen = (at for at in mol if at.properties.pdb_info.ResidueNumber != 1)
    xyz_array = to_array2(atom_gen)
    dist_array = cdist(xyz_array, xyz_array)

    a = mol[-1].properties.pdb_info.ResidueNumber - 1
    b = len(xyz_array) // a

    idx = np.empty((a, 2), dtype=float)
    E = np.empty((a), dtype=float)
    for i in range(a):
        c1, c2 = i*b, (i+1)*b
        dist_array[c1:c2, c1:c2] = np.nan
        idx[i, 0:2] = np.unravel_index(np.nanargmin(dist_array[c1:c2]), dist_array.shape)
        idx[i, 0] += c1
        E[i] = np.nanmin(dist_array[c1:c2])
    idx += (len(mol) - len(xyz_array)) + 1
    idx = np.hstack((idx, E[:, None]))

    return xyz_array, dist_array, idx[np.argsort(idx[:, 2])]


def get_all_dist2(mol):
    def to_array2(mol):
        x, y, z = zip(*[atom.coords for atom in atom_gen])
        return np.array((x, y, z)).T

    atom_gen = (at for at in mol if at.properties.pdb_info.ResidueNumber != 1)
    xyz_array = to_array2(atom_gen)
    dist_array = cdist(xyz_array, xyz_array)

    # a: number of ligand residues; b: number of atoms per ligand residue
    a = mol[-1].properties.pdb_info.ResidueNumber - 1
    b = len(xyz_array) // a

    for i in range(a):
        c1, c2 = i*b, (i+1)*b
        dist_array[c1:c2, c1:c2] = np.nan

    return dist_array


vec = (4, 4, 4)
rot = (0, -1, 0, 1, 0, 0, 0, 0, 1)
path = r'/Users/basvanbeek/Documents/CdSe/Week_5/QD/Core_Cd176Se147__58_Ligand_OC[C[CCC][CCC]CCCCCCCCCCCC]=O@O1.opt.pdb'
mol = molkit.readpdb(path)
# molkit.writepdb(mol2, '/Users/basvanbeek/Documents/CdSe/Week_5/test.pdb')
start = time.time()
dist_array = get_all_dist2(mol)

print('\n', time.time() - start)
name = 'np_benchmark.xlsx'
func = 'distance_to_mol'


"""
# print(mol1)
x, a, b = function_timeit(mol1)
plt.plot(x, a, 'r', label='numpy')
plt.plot(x, b, 'b', label='no numpy')
# print('filter more bonds:\t', "{0:.2f}".format((time.time() - start)*1000), '\tms')
"""

"""
new = pd.DataFrame({func + ' np': y, func: z})
spreadsheet = os.path.join(path, name)
if os.path.exists(spreadsheet):
    df = pd.read_excel(spreadsheet)
else:
    df = pd.DataFrame({'X': x})
for i in new:
    df[i] = new[i]
df.to_excel(spreadsheet)
"""
