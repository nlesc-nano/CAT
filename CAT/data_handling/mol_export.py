""" A module related to the exporting of molecules. """

__all__ = ['export_mol']

from os.path import join

import h5py

import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import get_time


def export_mol(mol, message='Mol:\t\t\t\t'):
    """
    Write results to a .pdb and .xyz file
    """
    name = mol.properties.name
    mol_path = join(mol.properties.path, name)
    pdb = Chem.MolToPDBBlock(molkit.to_rdmol(mol))
    mol.write(mol_path + '.xyz')
    print(get_time() + str(message) + str(name) + '.pdb')
