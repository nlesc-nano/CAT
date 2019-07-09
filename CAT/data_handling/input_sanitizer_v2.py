"""
CAT.data_handling.input_sanitizer
=================================

A module designed for sanitizing and interpreting the input file.

"""

import os
from os.path import (join, isfile, isdir, exists, basename)
from collections import abc
from typing import (Sequence, Any, Union, Optional)

import yaml
import numpy as np
import schema

from rdkit import Chem
from scm.plams import Settings, Molecule
import scm.plams.interfaces.molecule.rdkit as molkit

from CAT.utils import get_time

_INT = (int, np.integer)
_FLOAT = (int, float, np.integer, np.float)


def validate_path(path: Optional[str]) -> str:
    """Validate a provided directory path.

    Parameters
    ----------
    path : str
        Optional: A path to a directory.
        Will default to the current working directory if ``None``.

    Results
    -------
    Returns either the provided **path** parameter or the current working directory.

    Raises
    ------
    FileNotFoundError
        Raised if **path** cannot be found.

    NotADirectoryError
        Raised if **path** is not a directory.

    """
    if path is None:
        return os.getcwd()
    elif isdir(path):
        return path
    elif not exists(path):
        raise FileNotFoundError(get_time() + f"'{path}' not found")
    elif isfile(path):
        raise NotADirectoryError(get_time() + f"'{path}' is not a directory")


#: Schema for validating input molecules.
mol_schema = schema.Schema({
    schema.Optional('guess_bonds', default=False):
        bool,

    schema.Optional('is_core'):
        bool,

    schema.Optional('column'):
        schema.And(_INT, schema.Use(int), lambda n: n >= 0),

    schema.Optional('row'):
        schema.And(_INT, schema.Use(int), lambda n: n >= 0),

    schema.Optional('indices'):
        schema.And(abc.Iterable, lambda n: all(isinstance(i, _INT) for i in n)),

    schema.Optional('type'):
        object,

    schema.Optional('name'):
        str,

    schema.Optional('path'): schema.Use(validate_path)
})


def santize_smiles(smiles: str) -> str:
    """Sanitize a SMILES string: turn it into a valid filename."""
    name = smiles.replace('(', '[').replace(')', ']')
    cis_trans = [item for item in smiles if item in ('/', '\\')]
    if cis_trans:
        cis_trans = [item + cis_trans[i*2+1] for i, item in enumerate(cis_trans[::2])]
        cis_trans_dict = {'//': 'trans-', '/\\': 'cis-'}
        for item in cis_trans[::-1]:
            name = cis_trans_dict[item] + name
        name = name.replace('/', '').replace('\\', '')

    return name


def parse_mol_type(mol_dict: Settings):
    """ Sanitize and return the (file) type of the input molecule (SMILES, .xyz, dir, etc...). """
    mol = mol_dict.mol
    if isinstance(mol, str):
        if isfile(mol):
            mol_dict.type = mol.rsplit('.', 1)[-1]
            mol_dict.name = basename(mol.rsplit('.', 1)[0])
        elif isdir(mol):
            mol_dict.type = 'folder'
            mol_dict.name = basename(mol)
        else:
            mol_dict.type = 'smiles'
            mol_dict.mol = basename(mol)
            mol_dict.name = santize_smiles(basename(mol))

    elif isinstance(mol, Molecule):
        mol_dict.type = 'plams_mol'
        if not mol.properties.name:
            mol_dict.name = Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol))))
        else:
            mol_dict.name = mol.properties.name

    elif isinstance(mol, Chem.rdchem.Mol):
        mol_dict.type = 'rdmol'
        mol_dict.name = Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol.mol)))

    else:
        raise TypeError(f"mol_dict['mol'] expects an instance of 'str', 'Molecule' or 'Mol'; "
                        "observed type: {mol.__class__.__name__}")


def _parse_mol_type(mol_type: str) -> bool:
    """Parse the **mol_type** parameter of :func:`.validate_mol`."""
    if mol_type.lower() == 'input_cores':
        return True
    elif mol_type.lower() == 'input_ligands':
        return False
    else:
        raise ValueError(f"accepted values for mol_type are 'input_cores' and input_ligands; "
                         f"observed value: {repr(mol_type)}")


def validate_mol(args: Sequence[Union[Any, Settings]],
                 mol_type: str,
                 path: Optional[str] = None) -> None:
    r"""Validate the ``"input_ligands"`` or ``"input_cores"`` blocks in the CAT input.

    Performs an inpalce update of **args**.

    Examples
    --------

    An example using a list of .xyz files as input

    .. code:: python

        >>> print(args1)  # A list of .xyz files
        ['mol1.xyz', 'mol2.xyz']

        >>> validate_mol(args1, 'input_ligands')
        >>> print(args1[0], '\n', args1[1])
        is_core:        False
        mol:    mol1.xyz
        name:   mol1.xyz
        path:   /Users/bvanbeek/Documents/GitHub/CAT/CAT/data_handling
        type:   smiles

        is_core:        False
        mol:    mol2.xyz
        name:   mol2.xyz
        path:   /Users/bvanbeek/Documents/GitHub/CAT/CAT/data_handling
        type:   smiles


    An example using a list of .xyz-containing dictionaries as input

    .. code :: python

        >>> print(args2)  # A list of Settings instances with .xyz files
        [mol3.pdb:
             guess_bonds:       True
         mol4.pdb:
             guess_bonds:       True
        ]

        >>> validate_mol(args2, 'input_cores')
        >>> print(args2[0], '\n', args2[1])
        guess_bonds:    True
        is_core:        True
        mol:    mol3.pdb
        name:   mol3.pdb
        path:   /Users/bvanbeek/Documents/GitHub/CAT/CAT/data_handling
        type:   smiles

        guess_bonds:    True
        is_core:        True
        mol:    mol4.pdb
        name:   mol4.pdb
        path:   /Users/bvanbeek/Documents/GitHub/CAT/CAT/data_handling
        type:   smiles


    Parameters
    ----------
    args : |list|_ [|Settings|_]
        A list of input molecules.
        Accepts strings, PLAMS molecules and RDKit molecules.
        Additional arguments can be provided by putting above-mentioned molecules in a dictionary.

    mol_type : str
        The type of molecule.
        Accepted values are ``"input_ligands"`` and ``"input_cores"``.

    path : str
        Optional: The path to the molecule-containing directory.

    Raises
    ------
    FileNotFoundError
        Raised if **path** cannot be found.

    NotADirectoryError
        Raised if **path** is not a directory.

    ValueError
        Raised if the **mol_type** parameter is neither ``"input_cores"`` nor ``"input_ligands"``.

    SchemaError
        Raised if invalid input settings are found in while validating **args**.

    """
    # Validate arguments
    is_core = _parse_mol_type(mol_type)
    _path = validate_path(path)

    for i, item in enumerate(args):
        if not isinstance(item, dict):  # No optional arguments provided
            mol = item
            value = Settings({'path': _path, 'is_core': is_core})
        else:  # Pptional arguments have been provided: parse and validate them
            mol, value = next(iter(item.items()))
            value.setdefault('is_core', is_core)
            value = mol_schema.validate(value)
            value.setdefault('path', _path)

        if isinstance(mol, str) and not isfile(mol):
            mol = join(value.path, mol)
        value.mol = mol

        parse_mol_type(value)
        args[i] = value


args = Settings(yaml.load("""
    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False
""", Loader=yaml.FullLoader))

validate_mol(args.input_cores, 'input_cores', None)
