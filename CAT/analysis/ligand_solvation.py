""" A module designed for calculating solvation energies. """

__all__ = ['init_solv']

import os
from os.path import (join, dirname)

import numpy as np
import pandas as pd

from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.adf import ADFJob

from .crs import CRSJob
from .. import utils as CAT
from ..utils import get_time
from ..data_handling.database import (property_to_database, check_index)


def init_solv(mol_list, arg, solvent_list=None):
    """ Initialize the solvation energy calculation. """
    # Prepare the job settings and solvent list
    job_recipe = arg.optional.ligand.crs
    job_recipe.s1.input.mopac.Mozyme = 'Yes'
    if solvent_list is None:
        path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')
        solvent_list = [join(path, solv) for solv in os.listdir(path) if
                        solv not in ('__init__.py', 'README.rst')]
    solvent_list.sort()

    # Prepare the dataframe
    solv = ['E_solv', 'gamma'], [i.rsplit('.', 1)[0].rsplit('/', 1)[-1] for i in solvent_list]
    idx = pd.MultiIndex.from_product(solv, names=['index', 'sub index'])
    columns = pd.MultiIndex.from_tuples([(None, None)], names=['smiles', 'anchor'])
    df = pd.DataFrame(index=idx, columns=columns)

    # Check if the calculation has been donealready
    if 'ligand' not in arg.optional.database.overwrite:
        previous_entries = check_index('E_solv', arg, database='ligand')
        mol_list = [mol for mol in mol_list if
                    (mol.properties.smiles, mol.properties.anchor) not in previous_entries]
        if not mol_list:
            return

    # Run COSMO-RS
    init(path=mol_list[0].properties.path, folder='ligand_solvation')
    for mol in mol_list:
        coskf = get_surface_charge(mol, job=job_recipe.job1, s=job_recipe.s1)
        E_solv, gamma = get_solv(mol, solvent_list, coskf, job=job_recipe.job2, s=job_recipe.s2)

        key = mol.properties.smiles, mol.properties.anchor
        df[key] = None
        df[key]['E_solv'] = E_solv
        df[key]['gamma'] = gamma
    finish()

    # Update the database
    if 'ligand' in arg.optional.database.write:
        del df[(np.nan, np.nan)]
        property_to_database(df, arg, database='ligand')


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
    # Return 2x None if no coskf is None
    if coskf is None:
        return np.nan, np.nan

    # Prepare the job settings
    s.input.Compound._h = coskf
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
