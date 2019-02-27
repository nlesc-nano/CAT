""" A module related to the exporting of molecules. """

__all__ = ['export_mol']

from os.path import join

import scm.plams.interfaces.molecule.rdkit as molkit

from ..utils import get_time


def export_mol(mol, message='Mol:\t\t\t\t'):
    """
    Write results to a .pdb and .xyz file
    """
    name = mol.properties.name
    mol_path = join(mol.properties.path, name)
    molkit.writepdb(mol, mol_path + '.pdb')
    mol.write(mol_path + '.xyz')
    print(get_time() + str(message) + str(name) + '.pdb')
