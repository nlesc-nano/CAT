"""A module designed for calculating solvation energies."""

import os
import shutil
from itertools import product
from os.path import (join, dirname)
from typing import (Optional, Sequence, Callable, Container, Tuple, List, Iterable)

import numpy as np
import pandas as pd

from scm.plams import (Settings, Molecule, Results)
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.adf import ADFJob

import qmflows
from data_CAT import Database

from .crs import CRSJob
from .. import utils as CAT
from ..utils import (get_time, type_to_string)
from ..properties_dataframe import PropertiesDataFrame

__all__ = ['init_solv']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
JOB_SETTINGS_CRS = ('job_settings_crs', '')
SETTINGS1 = ('settings', 'solv 1')
SETTINGS2 = ('settings', 'solv 2')


def init_solv(ligand_df: PropertiesDataFrame,
              solvent_list: Optional[Sequence[str]] = None) -> None:
    """Initialize the ligand solvation energy calculation.

    Performs an inplace update of **ligand_df**, creating 2 sets of columns (*E_solv* & *gamma*)
    to hold all solvation energies and activity coefficients, respectively.

    Parameters
    ----------
    ligand_df : |CAT.PropertiesDataFrame|_
        A dataframe of ligands.

    solvent_list : |list|_ [|str|_]
        Optional: A list of paths to the .t21 or .coskf files of solvents.
        If ``None``, use the default .coskf files distributed with CAT (see :mod:`CAT.data.coskf`).

    """
    # Unpack arguments
    overwrite = 'ligand' in ligand_df.properties.optional.database.overwrite
    read = 'ligand' in ligand_df.properties.optional.database.read
    write = 'ligand' in ligand_df.properties.optional.database.write
    data = Database(path=ligand_df.properties.optional.database.dirname)

    # Prepare the job settings and solvent list
    solvent_list = get_solvent_list(solvent_list)

    # Update the columns of **ligand_df**
    columns = update_columns(ligand_df, solvent_list)

    # Check if the calculation has been done already
    if not overwrite and read:
        data.from_csv(ligand_df, database='ligand', get_mol=False)

    # Run COSMO-RS
    idx = ligand_df[['E_solv', 'gamma']].isna().all(axis='columns')
    if idx.any():
        start_crs_jobs(ligand_df)
        ligand_df[JOB_SETTINGS_CRS] = get_job_settings()
    else:
        return None  # No new molecules here; move along

    # Update the database
    if write:
        with pd.option_context('mode.chained_assignment', None):
            _ligand_to_db(ligand_df, idx, columns)
    return None


def start_crs_jobs(ligand_df: PropertiesDataFrame,
                   idx: pd.Series,
                   solvent_list: Iterable[str]) -> None:
    """Loop over all molecules in ``ligand_df.loc[idx]`` and perform COSMO-RS calculations."""
    # Unpack arguments
    path = ligand_df.properties.optional.ligand.dirname
    j1 = ligand_df.properties.optional.ligand.crs.job1
    j2 = ligand_df.properties.optional.ligand.crs.job2
    s1 = ligand_df.properties.optional.ligand.crs.s1
    s2 = ligand_df.properties.optional.ligand.crs.s2

    # Start the main loop
    init(path=path, folder='ligand_solvation')
    for i, mol in ligand_df[MOL][idx].iteritems():
        mol.properties.job_path = []

        # Calculate the COSMO surface
        coskf = get_surface_charge(mol, job=j1, s=s1)

        # Perform the actual COSMO-RS calculation
        e_and_gamma = get_solv(mol, solvent_list, coskf, job=j2, s=s2)
        ligand_df.loc[i, 'E_solv'], ligand_df.loc[i, 'gamma'] = e_and_gamma
    finish()


def update_columns(ligand_df: PropertiesDataFrame,
                   solvent_list: Iterable[str]) -> List[Tuple[str, str]]:
    """Add all COSMO-RS related columns to **ligand_df**."""
    clm_tups = [i.rsplit('.', 1)[0].rsplit('/', 1)[-1] for i in solvent_list]
    columns = list(product(('E_solv', 'gamma'), clm_tups))
    for item in columns:
        ligand_df[item] = np.nan
    return columns


def get_solvent_list(solvent_list: Optional[Sequence[str]]) -> Sequence[str]:
    """Construct the list of solvents."""
    if solvent_list is None:
        coskf_path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')
        solvent_list = [join(coskf_path, solv) for solv in os.listdir(coskf_path) if
                        solv not in ('__init__.py', 'README.rst')]
    solvent_list.sort()
    return solvent_list


def get_job_settings(ligand_df: PropertiesDataFrame) -> List[str]:
    """Create a nested list of input files for each molecule in **ligand_df**."""
    job_settings = []
    for mol in ligand_df[MOL]:
        try:
            job_settings.append(mol.properties.pop('job_path'))
        except KeyError:
            job_settings.append([])
    return job_settings


def _ligand_to_db(ligand_df: PropertiesDataFrame,
                  idx: pd.Series,
                  columns: Sequence) -> None:
    """Export all COSMO-RS results to the database."""
    data = Database(path=ligand_df.properties.optional.database.dirname)
    overwrite = 'ligand' in ligand_df.properties.optional.database.overwrite
    j1 = ligand_df.properties.optional.ligand.crs.job1
    j2 = ligand_df.properties.optional.ligand.crs.job2
    s1 = ligand_df.properties.optional.ligand.crs.s1
    s2 = ligand_df.properties.optional.ligand.crs.s2

    value1 = qmflows.singlepoint['specific'][type_to_string(j1)].copy()
    value1.update(s1)
    recipe = Settings()
    recipe['solv 1'] = {'key': j1, 'value': value1}
    recipe['solv 2'] = {'key': j2, 'value': s2}

    data.update_csv(
        ligand_df.loc[idx],
        database='ligand',
        columns=[SETTINGS1, SETTINGS2, JOB_SETTINGS_CRS]+columns,
        overwrite=overwrite,
        job_recipe=recipe
    )


def get_surface_charge(mol: Molecule,
                       job: Callable,
                       s: Settings) -> Optional[str]:
    """Construct the COSMO surface of the **mol**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS Molecule.

    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    s : |plams.Settings|_
        The settings for **job**.

    Returns
    -------
    |plams.Settings|_
        Optional: The path+filename of a file containing COSMO surface charges.

    """
    # Special procedure for ADF jobs
    # Use the gas-phase electronic structure as a fragment for the COSMO single point
    if job is ADFJob:
        s = get_surface_charge_adf(mol, job, s)

    s.runscript.post = '$ADFBIN/cosmo2kf "mopac.cos" "mopac.coskf"'
    results = mol.job_single_point(job, s, ret_results=True)
    results.wait()
    return get_coskf(results)


def get_solv(mol: Molecule,
             solvent_list: Iterable[str],
             coskf: str,
             job: Callable,
             s: Settings,
             keep_files: bool = True) -> Tuple[List[float], List[float]]:
    """Calculate the solvation energy of *mol* in various *solvents*.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS Molecule.

    solvent_list : |List|_ [|str|_]
        A list of solvent molecules (*i.e.* .coskf files).

    coskf : str
        The path+filename of the .coskf file of **mol**.

    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    s : |plams.Settings|_
        The settings for **job**.

    Returns
    -------
    |list|_ [|float|_] & |list|_ [|float|_]
        A list of solvation energies and gammas.

    """
    # Return 2x np.nan if no coskf is None (i.e. the COSMO-surface construction failed)
    if coskf is None:
        return np.nan, np.nan

    # Prepare a list of job settings
    s.input.Compound._h = coskf
    s.ignore_molecule = True
    s_list = []
    for solv in solvent_list:
        s_tmp = s.copy()
        s_tmp.name = solv.rsplit('.', 1)[0].rsplit('/', 1)[-1]
        s_tmp.input.compound._h = solv
        s_list.append(s_tmp)

    # Run the job
    jobs = [CRSJob(settings=s, name=s.name) for s in s_list]
    results = [j.run(jobrunner=JobRunner(parallel=True)) for j in jobs]

    # Extract solvation energies and activity coefficients
    E_solv = []
    Gamma = []
    for result in results:
        result.wait()
        try:
            E_solv.append(result.get_energy())
            Gamma.append(result.get_activity_coefficient())
        except ValueError:
            print(get_time() + 'WARNING: Failed to retrieve COSMO-RS results of ' +
                  results.job.name)
            E_solv.append(np.nan)
            Gamma.append(np.nan)

    # Delete all mopac and cosmo-rs files if keep_files=False
    if not keep_files:
        mopac = dirname(s.input.Compound._h)
        shutil.rmtree(mopac)
        for job in jobs:
            shutil.rmtree(job.path)

    if 'job_path' not in mol.properties:
        mol.properties.job_path = []
    mol.properties.job_path += [join(j.path, j.name + '.in') for j in jobs]

    # Return the solvation energies and activity coefficients as dict
    return E_solv, Gamma


def get_surface_charge_adf(mol: Molecule,
                           job: Callable,
                           s: Settings) -> Settings:
    """Perform a gas-phase ADF single point and return settings for a COSMO-ADF single point.

    The previous gas-phase calculation as moleculair fragment.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS Molecule.

    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    s : |plams.Settings|_
        The settings for **job**.

    Returns
    -------
    |plams.Settings|_
        A new Settings intance, constructed from **s**, suitable for DFT COSMO-RS calculations.

    """
    s.input.allpoints = ''
    s.input.charge = sum([at.properties.charge for at in mol])
    results = mol.job_single_point(job, s, ret_results=True)
    coskf = get_coskf(results)

    for at in mol:
        at.properties.adf.fragment = 'gas'
    s.update(CAT.get_template('qd.yaml')['COSMO-ADF'])
    s.input.fragments.gas = coskf

    return s


def get_coskf(results: Results,
              extensions: Container[str] = ['.coskf', '.t21']) -> Optional[str]:
    """Return the file in **results** containing the COSMO surface.

    Parameters
    ----------
    results : |plams.Results|_
        A Results instance.

    extensions : |list|_ [|str|_]
        Valid filetypes which can contain COSMO surfaces.

    Returns
    -------
        Optional: The path+filename of a file containing COSMO surface charges.

    """
    for file in results.files:
        for ext in extensions:
            if ext in file:
                return results[file]
    print(get_time() + 'WARNING: Failed to retrieve COSMO surface charges of ' + results.job.name)
    return None
