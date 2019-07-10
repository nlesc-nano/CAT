"""
CAT.data_handling.mol_import
============================

A module related to the importing of molecules.

Index
-----
.. currentmodule:: CAT.data_handling.mol_import
.. autosummary::
    read_mol
    read_mol_xyz
    read_mol_pdb
    read_mol_mol
    read_mol_smiles
    read_mol_plams
    read_mol_rdkit
    read_mol_folder
    read_mol_txt
    get_charge_dict
    set_mol_prop
    set_atom_prop
    print_exception

API
---
.. autofunction:: CAT.data_handling.mol_import.read_mol
.. autofunction:: CAT.data_handling.mol_import.read_mol_xyz
.. autofunction:: CAT.data_handling.mol_import.read_mol_pdb
.. autofunction:: CAT.data_handling.mol_import.read_mol_mol
.. autofunction:: CAT.data_handling.mol_import.read_mol_smiles
.. autofunction:: CAT.data_handling.mol_import.read_mol_plams
.. autofunction:: CAT.data_handling.mol_import.read_mol_rdkit
.. autofunction:: CAT.data_handling.mol_import.read_mol_folder
.. autofunction:: CAT.data_handling.mol_import.read_mol_txt
.. autofunction:: CAT.data_handling.mol_import.get_charge_dict
.. autofunction:: CAT.data_handling.mol_import.set_mol_prop
.. autofunction:: CAT.data_handling.mol_import.set_atom_prop
.. autofunction:: CAT.data_handling.mol_import.print_exception

"""

import os
import itertools
from string import ascii_letters
from typing import (Dict, Iterable, List, Callable, Sequence)

from scm.plams import (Molecule, Atom, Settings)
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..utils import get_time
from ..data_handling.validate_mol import validate_mol

__all__ = ['read_mol', 'set_mol_prop']


def read_mol(input_mol: Iterable[Settings]) -> List[Molecule]:
    """Checks the filetypes of the input molecules.

    Sets the molecules' properties and returns a list of plams molecules.

    Parameters
    ----------
    input_mol : |list|_ [|Settings|_]
        An iterable consisting of dictionaries with input settings per mol.

    Returns
    -------
    |plams.Molecule|_
        A list of plams Molecules.

    """
    # Creates a dictionary of file extensions
    extension_dict = {
        'xyz': read_mol_xyz,
        'pdb': read_mol_pdb,
        'mol': read_mol_mol,
        'smiles': read_mol_smiles,
        'folder': read_mol_folder,
        'txt': read_mol_txt,
        'plams_mol': read_mol_plams,
        'rdmol': read_mol_rdkit
    }

    # Create a list of PLAMS molecules
    mol_list = []
    for mol_dict in input_mol:
        try:
            read_mol = extension_dict[mol_dict.type]
        except KeyError as ex:
            print(get_time() + f'{ex.__class__.__name__}:\t{ex}\n')
            continue

        mol = read_mol(mol_dict)
        if not mol:  # Failed to import any molecules
            continue

        if isinstance(mol, Molecule):  # if mol is a PLAMS molecule
            if mol_dict.guess_bonds:
                mol.guess_bonds()
            set_mol_prop(mol, mol_dict)
            mol_list.append(mol)
        else:  # if mol is a list of molecules
            mol_list += mol

    return mol_list


def read_mol_xyz(mol: Settings) -> Molecule:
    """Read an .xyz file."""
    try:
        return Molecule(mol.mol, inputformat='xyz')
    except Exception as ex:
        print_exception(read_mol_xyz.__code__, ex, mol.mol)


def read_mol_pdb(mol: Settings) -> Molecule:
    """Read a .pdb file."""
    try:
        return molkit.readpdb(mol.mol)
    except Exception as ex:
        print_exception(read_mol_pdb.__code__, ex, mol.mol)


def read_mol_mol(mol: Settings) -> Molecule:
    """Read a .mol file."""
    try:
        return molkit.from_rdmol(Chem.MolFromMolFile(mol.mol, removeHs=False))
    except Exception as ex:
        print_exception(read_mol_mol.__code__, ex, mol.mol)


def read_mol_smiles(mol: Settings) -> Molecule:
    """Read a SMILES string."""
    try:
        return molkit.from_smiles(mol.mol)
    except Exception as ex:
        print_exception(read_mol_smiles.__code__, ex, mol.mol)


def read_mol_plams(mol: Settings) -> Molecule:
    """Read a PLAMS molecule."""
    try:
        return mol.mol
    except Exception as ex:
        print_exception(read_mol_plams.__code__, ex, mol.mol)


def read_mol_rdkit(mol: Settings) -> Molecule:
    """Read a RDKit molecule."""
    try:
        return molkit.from_rdmol(mol.mol)
    except Exception as ex:
        print_exception(read_mol_rdkit.__code__, ex, mol.mol)


def read_mol_folder(mol: Settings) -> Molecule:
    """Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder."""
    try:
        file_list = [file for file in os.listdir(mol.mol)]
        mol_type = 'input_cores' if mol.is_core else 'input_ligands'
        input_mol = validate_mol(file_list, mol.path, mol_type)
        return read_mol(input_mol)
    except Exception as ex:
        print_exception(read_mol_folder.__code__, ex, mol.mol)


def read_mol_txt(mol: Settings) -> Molecule:
    """Read a plain text file containing one or more SMILES strings."""
    try:
        with open(mol.mol, 'r') as file:
            file_list = file.read().splitlines()
        file_list = [file.split()[mol.column] for file in file_list[mol.row:] if file]
        mol_type = 'input_cores' if mol.is_core else 'input_ligands'
        input_mol = validate_mol(file_list, mol.path, mol_type)
        return read_mol(input_mol)
    except Exception as ex:
        print_exception(read_mol_txt.__code__, ex, mol.mol)


def get_charge_dict() -> Dict[str, int]:
    """Create a dictionary of elements and their formal atomic charge."""
    # Create a list of atomic charges and elements
    charges = (1, 2, -3, -2, -1, 2)
    elements = (
        ('Li', 'Na', 'K', 'Rb', 'Cs'),      # Group 1: Alkali metals
        ('Be', 'Mg', 'Ca', 'Sr', 'Ba'),     # Group 2: Alkaline earth metals
        ('N', 'P', 'As', 'Sb', 'Bi'),       # Group 15: pnictogens
        ('O', 'S', 'Se', 'Te', 'Po'),       # Group 16: Chalcogens
        ('H', 'F', 'Cl', 'Br', 'I', 'At'),  # Group 17: Halogens
        ('Cd', 'Pb')                        # Misc
    )

    # Combine the elements and atomic charges into a dictionary
    values = itertools.chain.from_iterable([i] * len(j) for i, j in zip(charges, elements))
    keys = itertools.chain(*elements)

    return dict(zip(keys, values))


# Create a dictionary of elements and their formal atomic charge
charge_dict: Dict[str, int] = get_charge_dict()


def set_mol_prop(mol: Molecule,
                 mol_dict: Settings) -> None:
    """Set molecular and atomic properties."""
    if mol_dict.is_core:
        residue_name = 'COR'
        mol.properties.name = mol.get_formula()
    else:
        residue_name = 'LIG'
        mol.properties.name = mol_dict.name

    mol.properties.dummies = mol_dict.indices
    mol.properties.path = mol_dict.path
    mol.properties.job_path = []

    # Prepare a generator of letters for pdb_info.Name
    alphabet = itertools.combinations(ascii_letters, 2)

    # Set the atomic properties
    for atom, i in zip(mol, alphabet):
        set_atom_prop(atom, i, residue_name)

    if not mol.properties.smiles:
        tmp = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))
        mol.properties.smiles = Chem.CanonSmiles(tmp)


def set_atom_prop(atom: Atom,
                  i: Sequence[str],
                  residue_name: str) -> None:
    """Set atomic properties."""
    symbol = '{:4}'.format(atom.symbol + ''.join(i))

    # Add a number of properties to atom
    atom.properties.pdb_info.ResidueName = residue_name
    atom.properties.pdb_info.Occupancy = 1.0
    atom.properties.pdb_info.TempFactor = 0.0
    atom.properties.pdb_info.ResidueNumber = 1
    atom.properties.pdb_info.Name = symbol
    atom.properties.pdb_info.ChainId = 'A'

    # Changes hydrogen and carbon from heteroatom to atom
    if atom.symbol in ('H', 'C'):
        atom.properties.pdb_info.IsHeteroAtom = False
    else:
        atom.properties.pdb_info.IsHeteroAtom = True

    # Sets the formal atomic charge
    if not atom.properties.charge:
        if atom.symbol in charge_dict:
            total_bonds = int(sum([bond.order for bond in atom.bonds]))
            default_charge = charge_dict[atom.symbol]
            abs_charge = abs(default_charge)
            sign = -1 * int(default_charge / abs_charge)
            atom.properties.charge = default_charge + sign*total_bonds

            # Update formal atomic charges for hypervalent atoms
            if total_bonds > abs_charge:
                if total_bonds is abs_charge + 2:
                    atom.properties.charge += 2 * sign
                elif total_bonds is abs_charge + 4:
                    atom.properties.charge += 4 * sign
                elif total_bonds >= abs_charge + 6:
                    atom.properties.charge += 6 * sign
        else:
            atom.properties.charge = 0


def print_exception(func: Callable,
                    ex: Exception,
                    name: str) -> None:
    """Manages the printing of exceptions upon failing to import a molecule."""
    extension_dict = {'read_mol_xyz': '.xyz file', 'read_mol_pdb': '.pdb file',
                      'read_mol_mol': '.mol file', 'read_mol_smiles': 'SMILES string',
                      'read_mol_folder': 'folder', 'read_mol_txt': '.txt file',
                      'read_mol_excel': '.xlsx file', 'read_mol_plams': 'PLAMS molecule',
                      'read_mol_rdkit': 'RDKit molecule'}
    print(get_time() + str(type(ex).__name__), str(ex))
    print(get_time() + 'Warning:', name, 'not recognized as a valid',
          extension_dict[func.co_name], '\n')
