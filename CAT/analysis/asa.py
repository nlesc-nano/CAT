"""A module related to performing activation strain analyses."""

import numpy as np
import pandas as pd

from scm.plams.core.settings import Settings
import scm.plams.interfaces.molecule.rdkit as molkit

import rdkit
from rdkit.Chem import AllChem

from data_CAT import Database

from ..properties_dataframe import PropertiesDataFrame

__all__ = ['init_asa']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')
SETTINGS1 = ('settings', 'ASA 1')


def init_asa(qd_df: PropertiesDataFrame) -> None:
    """Initialize the activation-strain analyses (ASA).

    The ASA (RDKit UFF level) is conducted on the ligands in the absence of the core.

    Parameters
    ----------
    |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    """
    # Unpack arguments
    overwrite = 'qd' in qd_df.properties.optional.database.overwrite
    write = 'qd' in qd_df.properties.optional.database.write
    data = Database(qd_df.properties.optional.database.dirname)

    # Prepare columns
    columns = [ASA_INT, ASA_STRAIN, ASA_E]
    for i in columns:
        qd_df[i] = np.nan

    # Fill columns
    qd_df['ASA'] = get_asa_energy(qd_df[MOL])

    # Calculate E_int, E_strain and E
    if write:
        recipe = Settings()
        recipe['ASA 1'] = {'key': 'RDKit_' + rdkit.__version__, 'value': 'UFF'}
        data.update_csv(qd_df,
                        columns=[SETTINGS1]+columns,
                        job_recipe=recipe,
                        database='QD',
                        overwrite=overwrite)


def get_asa_energy(mol_series: pd.Series) -> np.ndarray:
    """Perform an activation strain analyses (ASA).

    The ASA calculates the interaction, strain and total energy.
    The ASA is performed on all ligands in the absence of the core at the UFF level (RDKit).

    Parameters
    ----------
    mol_series : |pd.Series|_
        A series of PLAMS molecules.

    Returns
    -------
    :math:`n*3` |np.ndarray|_ [|np.float64|_]
        An array containing E_int, E_strain and E for all *n* molecules in **mol_series**.

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
