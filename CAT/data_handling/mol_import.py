""" A module related to the importing of molecules. """

__all__ = ['read_mol', 'set_prop']

import os
import itertools

from scm.plams.mol.molecule import Molecule
from scm.plams.core.errors import PlamsError
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import get_time
from ..data_handling.input_sanitizer import (sanitize_mol_type, get_mol_defaults)


def read_mol(input_mol):
    """
    Checks the filetypes of the input molecules, sets their properties and
    returns a list of plams molecules.
    """
    # Creates a dictionary of file extensions
    extension_dict = {'xyz': read_mol_xyz,
                      'pdb': read_mol_pdb,
                      'mol': read_mol_mol,
                      'smiles': read_mol_smiles,
                      'folder': read_mol_folder,
                      'txt': read_mol_txt,
                      'plams_mol': read_mol_plams,
                      'rdmol': read_mol_rdkit}

    # Create a list of PLAMS molecules
    mol_list = []
    for mol_dict in input_mol:
        try:
            read_mol = extension_dict[mol_dict.type]
        except KeyError as ex:
            print(get_time() + str(type(ex).__name__) + ':\t' + str(ex) + '\n')

        # Convert mol into either a list or PLAMS molecule
        if read_mol:
            mol = read_mol(mol_dict)
            if mol:
                # if mol is a list
                if isinstance(mol, list):
                    mol_list += mol

                # if mol is a PLAMS molecule
                else:
                    if mol_dict.guess_bonds:
                        mol.guess_bonds()
                    set_prop(mol, mol_dict)
                    mol_list.append(mol)

    return mol_list


def read_mol_xyz(mol):
    """ Read an .xyz file """
    try:
        return Molecule(mol.mol, inputformat='xyz')
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_xyz.__code__, ex, mol.mol)


def read_mol_pdb(mol):
    """ Read a .pdb file """
    try:
        return molkit.readpdb(mol.mol)
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_pdb.__code__, ex, mol.mol)


def read_mol_mol(mol, mol_dict):
    """ Read a .mol file """
    try:
        return molkit.from_rdmol(Chem.MolFromMolFile(mol.mol, removeHs=False))
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_mol.__code__, ex, mol.mol)


def read_mol_smiles(mol):
    """ Read a SMILES string """
    try:
        return molkit.from_smiles(mol.mol)
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_smiles.__code__, ex, mol.mol)


def read_mol_plams(mol):
    """ Read a PLAMS molecule """
    try:
        return mol.mol
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_plams.__code__, ex, mol.mol)


def read_mol_rdkit(mol):
    """ Read a RDKit molecule """
    try:
        return molkit.from_rdmol(mol.mol)
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_rdkit.__code__, ex, mol.mol)


def read_mol_folder(mol):
    """ Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder """
    try:
        file_list = [file for file in os.listdir(mol.mol)]
        input_mol = get_mol_defaults(file_list, path=mol.path, core=mol.is_core)
        input_mol = sanitize_mol_type(input_mol)
        return read_mol(input_mol)
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_folder.__code__, ex, mol.mol)


def read_mol_txt(mol):
    """ Read a plain text file containing one or more SMILES strings """
    try:
        with open(mol.mol, 'r') as file:
            file_list = file.read().splitlines()
        file_list = [file.split()[mol.column] for file in file_list[mol.row:] if file]
        input_mol = get_mol_defaults(file_list, path=mol.path, core=mol.is_core)
        input_mol = sanitize_mol_type(input_mol)
        return read_mol(input_mol)
    except (Exception, PlamsError) as ex:
        print_exception(read_mol_txt.__code__, ex, mol.mol)


def set_prop(mol, mol_dict):
    """
    Set molecular and atomic properties
    """
    if mol_dict.is_core:
        residue_name = 'COR'
        mol.properties.name = mol.get_formula()
    else:
        residue_name = 'LIG'
        mol.properties.name = mol_dict.name

    mol.properties.dummies = mol_dict.indices
    mol.properties.path = mol_dict.path

    # Prepare a list of letters for pdb_info.Name
    alphabet = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphabet = [i + j for i in alphabet for j in alphabet]

    # Create a dictionary of elements and their formal atomic charge
    elements_dict = set_prop_dict()

    # Set the atomic properties
    for i, atom in enumerate(mol):
        set_prop_atom(atom, alphabet[i], residue_name, elements_dict)

    if not mol.properties.smiles:
        mol.properties.smiles = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))
        mol.properties.smiles = Chem.CanonSmiles(mol.properties.smiles)


def set_prop_atom(atom, alphabet, residue_name, elements_dict):
    """
    Set atomic properties.
    """
    symbol = atom.symbol + alphabet + '  '

    # Add a number of properties to atom
    atom.properties.pdb_info.ResidueName = residue_name
    atom.properties.pdb_info.Occupancy = 1.0
    atom.properties.pdb_info.TempFactor = 0.0
    atom.properties.pdb_info.ResidueNumber = 1
    atom.properties.pdb_info.Name = symbol[:4]
    atom.properties.pdb_info.ChainId = 'A'

    # Changes hydrogen and carbon from heteroatom to atom
    if atom.symbol in ('H', 'C'):
        atom.properties.pdb_info.IsHeteroAtom = False
    else:
        atom.properties.pdb_info.IsHeteroAtom = True

    # Sets the formal atomic charge
    if atom.properties.charge:
        return

    if atom.symbol in elements_dict:
        total_bonds = int(sum([bond.order for bond in atom.bonds]))
        default_charge = elements_dict[atom.symbol]
        sign = int(-1 * default_charge / abs(default_charge))
        atom.properties.charge = default_charge + sign*total_bonds

        # Update formal atomic charges for hypervalent atoms
        abs_charge = abs(default_charge)
        if total_bonds > abs_charge:
            if total_bonds is abs_charge + 2:
                atom.properties.charge -= sign*2
            elif total_bonds is abs_charge + 4:
                atom.properties.charge -= sign*4
            elif total_bonds >= abs_charge + 6:
                atom.properties.charge -= sign*6
    else:
        atom.properties.charge = 0


def set_prop_dict():
    """
    Create a dictionary of elements and their formal atomic charge.
    """
    # Create a list of atomic charges and elements
    charges = [1, 2, -3, -2, -1, 2]
    group01 = ['Li', 'Na', 'K', 'Rb', 'Cs']       # Alkaline metals
    group02 = ['Be', 'Mg', 'Ca', 'Sr', 'Ba']      # Alkaline earth metals
    group15 = ['N', 'P', 'As', 'Sb', 'Bi']        # Pnictogens
    group16 = ['O', 'S', 'Se', 'Te', 'Po']        # Chalcogens
    group17 = ['H', 'F', 'Cl', 'Br', 'I', 'At']   # Halogens
    misc = ['Cd', 'Pb']                        # Misc

    # Combine the elements and atomic charges into a dictionary
    elements = [group01, group02, group15, group16, group17, misc]
    charges = [charges[i] for i, column in enumerate(elements) for element in column]
    elements = list(itertools.chain(*elements))

    return dict(zip(elements, charges))


def print_exception(func, ex, name):
    """
    Manages the printing of exceptions upon failing to import a molecule.
    """
    extension_dict = {'read_mol_xyz': '.xyz file', 'read_mol_pdb': '.pdb file',
                      'read_mol_mol': '.mol file', 'read_mol_smiles': 'SMILES string',
                      'read_mol_folder': 'folder', 'read_mol_txt': '.txt file',
                      'read_mol_excel': '.xlsx file', 'read_mol_plams': 'PLAMS molecule',
                      'read_mol_rdkit': 'RDKit molecule'}
    print(get_time() + str(type(ex).__name__), str(ex))
    print(get_time() + 'Warning:', name, 'not recognized as a valid',
          extension_dict[func.co_name], '\n')
    return []
