""" A module designed for calculating solvation energies. """

__all__ = ['init_solv']

import os
from itertools import product
from os.path import (join, dirname)

import numpy as np

from scm.plams.core.settings import Settings
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.adf import ADFJob

from .crs import CRSJob
from .. import utils as CAT
from ..utils import get_time
from ..data_handling.CAT_database import Database


def init_solv(ligand_df, arg, solvent_list=None):
    """ Initialize the ligand solvation energy calculation.
    Performs an inplace update of **ligand_df**, creating 2 sets of columns (*E_solv* & *gamma*)
    to hold all solvation energies and activity coefficients, respectively.

    :parameter ligand_df: A dataframe of ligands.
    :type ligand_df: |pd.DataFrame|_ (columns: |str|_, index=|int|_, values=|plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
    :parameter solvent_list: A list of paths to the .t21 or .coskf files of solvents. If *None*,
        use the default .coskf files distributed with CAT (see CAT.data.coskf).
    :type solvent_list: |None|_ or |list|_ [|str|_].
    """
    data = Database(path=arg.optional.database.dirname)

    # Prepare the job settings and solvent list
    j1, j2 = arg.optional.ligand.crs.job1, arg.optional.ligand.crs.job2
    s1, s2 = arg.optional.ligand.crs.s1, arg.optional.ligand.crs.s2
    if solvent_list is None:
        path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')
        solvent_list = [join(path, solv) for solv in os.listdir(path) if
                        solv not in ('__init__.py', 'README.rst')]
    solvent_list.sort()

    # Update the columns of **ligand_df**
    columns = [i.rsplit('.', 1)[0].rsplit('/', 1)[-1] for i in solvent_list]
    columns = list(product(('E_solv', 'gamma'), columns))
    for item in columns:
        ligand_df[item] = np.nan

    # Check if the calculation has been done already
    overwrite = 'ligand' in arg.optional.database.overwrite
    if not overwrite and 'ligand' in arg.optional.database.read:
        data.from_csv(ligand_df, database='ligand', get_mol=False)

    # Run COSMO-RS
    idx = ligand_df[['E_solv', 'gamma']].isna().all(axis='columns')
    if idx.any():
        init(path=ligand_df['mol'][0].properties.path, folder='ligand_solvation')
        for i, mol in ligand_df['mol'][idx].iteritems():
            coskf = get_surface_charge(mol, job=j1, s=s1)
            ligand_df.loc[i, 'E_solv'], ligand_df.loc[i, 'gamma'] = get_solv(mol, solvent_list,
                                                                             coskf, job=j2, s=s2)
        finish()
        print('')

    # Update the database
    if 'ligand' in arg.optional.database.write:
        recipe = Settings()
        recipe.settings1 = {'name': 'solv 1', 'key': j1, 'value': s1,
                            'template': 'singlepoint.json'}
        recipe.settings2 = {'name': 'solv 2', 'key': j2, 'value': s2}
        data.update_csv(ligand_df, database='ligand',
                        columns=[('settings', 'solv 1'), ('settings', 'solv 2')]+columns,
                        overwrite=overwrite, job_recipe=recipe)


def get_surface_charge(mol, job=None, s=None):
    """ Construct the COSMO surface of the *mol*. """
    # Special procedure for ADF jobs
    # Use the gas-phase electronic structure as a fragment for the COSMO single point
    if job is ADFJob:
        s = get_surface_charge_adf(mol, job, s)

    results = mol.job_single_point(job, s, ret_results=True)
    results.wait()
    return get_coskf(results)


def get_solv(mol, solvent_list, coskf, job=None, s=None):
    """ Calculate the solvation energy of *mol* in various *solvents*. """
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
    results = [job.run(jobrunner=JobRunner(parallel=True)) for job in jobs]

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

    # Return the solvation energies and activity coefficients as dict
    return E_solv, Gamma


def get_surface_charge_adf(mol, job, s):
    """ Perform a gas-phase ADF single point and return settings for a
    COSMO-ADF single point, using the previous gas-phase calculation as moleculair fragment. """
    s.input.allpoints = ''
    s2 = s.copy()
    del s2.input.solvation

    results = mol.job_single_point(job, s, ret_results=True)
    coskf = results[get_coskf(results)]

    for at in mol:
        at.properties.adf.fragment = 'gas'
    s.update(CAT.get_template('qd.json')['COSMO-ADF'])
    s.input.fragments.gas = coskf

    return s


def get_coskf(results, extensions=['.coskf', '.t21']):
    """ Return the file in results containing the COSMO surface. """
    for file in results.files:
        for ext in extensions:
            if ext in file:
                return results[file]
    print(get_time() + 'WARNING: Failed to retrieve COSMO surface charges of ' + results.job.name)
    return None
