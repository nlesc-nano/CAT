""" A module related to performing activation strain analyses. """

__all__ = ['init_asa']

import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit.Chem import AllChem


def init_asa(plams_mol):
    """ Perform an activation-strain analyses (RDKit UFF) on the ligands in the absence of the core.

    plams_mol <plams.Molecule>: A PLAMS molecule.
    return <plams.Molecule>: A PLAMS molecule with the int and int_mean properties.
    """
    mol_copy = plams_mol.copy()
    uff = AllChem.UFFGetMoleculeForceField

    # Calculate the total energy of all perturbed ligands in the absence of the core
    for atom in reversed(mol_copy.atoms):
        if atom.properties.pdb_info.ResidueName == 'COR':
            mol_copy.delete_atom(atom)

    rdmol = molkit.to_rdmol(mol_copy)
    E_no_frag = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate the total energy of the isolated perturbed ligands in the absence of the core
    mol_frag = mol_copy.separate()
    E_frag = 0.0
    for mol in mol_frag:
        rdmol = molkit.to_rdmol(mol)
        E_frag += uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate the total energy of the optimized ligand
    uff(rdmol, ignoreInterfragInteractions=False).Minimize()
    E_opt = uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

    # Calculate E, Eint and Estrain
    plams_mol.properties.energy.Eint = float(E_no_frag - E_frag)
    plams_mol.properties.energy.Estrain = float(E_frag - (E_opt * len(mol_frag)))
    plams_mol.properties.energy.E = plams_mol.properties.energy.Eint + plams_mol.properties.energy.Estrain

    return plams_mol
