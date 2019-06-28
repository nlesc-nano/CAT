""" A module related to performing activation strain analyses. """

__all__ = ['init_asa']

import numpy as np

from scm.plams.core.settings import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

import rdkit
from rdkit.Chem import AllChem

from data_CAT import Database


def init_asa(qd_df, arg):
    """ Initialize the activation-strain analyses (RDKit UFF level) on the ligands in the
    absence of the core.

    :parameter qd_df: A dataframe of quantum dots.
    :type qd_df: |pd.DataFrame|_ (columns: |str|_, index=|int|_, values=|plams.Molecule|_)
    """
    data = Database(arg.optional.database.dirname)
    overwrite = 'qd' in arg.optional.database.overwrite

    # Prepare columns
    columns = [('ASA', 'E_int'), ('ASA', 'E_strain'), ('ASA', 'E')]
    for i in columns:
        qd_df[i] = np.nan

    # Fill columns
    qd_df['ASA'] = get_asa_energy(qd_df['mol'])

    # Calculate E_int, E_strain and E
    if 'qd' in arg.optional.database.write:
        recipe = Settings()
        recipe['ASA 1'] = {'key': 'RDKit_' + rdkit.__version__, 'value': 'UFF'}
        data.update_csv(qd_df, columns=[('settings', 'ASA 1')]+columns,
                        job_recipe=recipe, database='QD', overwrite=overwrite)


def get_asa_energy(mol_series):
    """ Calculate the interaction, strain and total energy in the framework of the
    activation-strain analysis (ASA).
    The ASA is performed on all ligands in the absence of the core at the UFF level (RDKit).

    :parameter mol_series: A series of PLAMS molecules.
    :type mol_series: |pd.Series|_ (index=|str|_, values: |plams.Molecule|_)
    :return: An array containing E_int, E_strain and E for all *n* molecules in **mol_series**.
    :rtype: *n*3* |np.ndarray|_ [|np.float64|_]
    """
    ret = np.zeros((len(mol_series), 4))

    for i, mol in enumerate(mol_series):
        mol_cp = mol.copy()
        rd_uff = AllChem.UFFGetMoleculeForceField

        # Calculate the total energy of all perturbed ligands in the absence of the core
        for atom in reversed(mol_cp.atoms):
            if atom.properties.pdb_info.ResidueName == 'COR':
                mol_cp.delete_atom(atom)
        rdmol = molkit.to_rdmol(mol_cp)
        E_no_frag = rd_uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

        # Calculate the total energy of the isolated perturbed ligands in the absence of the core
        mol_frag = mol_cp.separate()
        E_frag = 0.0
        for plams_mol in mol_frag:
            rdmol = molkit.to_rdmol(plams_mol)
            E_frag += rd_uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

        # Calculate the total energy of the optimized ligand
        rd_uff(rdmol, ignoreInterfragInteractions=False).Minimize()
        E_opt = rd_uff(rdmol, ignoreInterfragInteractions=False).CalcEnergy()

        # Update ret with the new activation strain terms
        ret[i] = E_no_frag, E_frag, E_opt, len(mol_frag)

    # Post-process and return
    ret[:, 0] -= ret[:, 1]
    ret[:, 1] -= ret[:, 2] * ret[:, 3]
    ret[:, 2] = ret[:, 0] + ret[:, 1]
    return ret[:, 0:3]
