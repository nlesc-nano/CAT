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


def init_solv(mol_list, arg, solvent_list=None):
    """ Initialize the solvation energy calculation. """
    # Prepare the job settings and solvent list
    job_recipe = arg.optional.ligand.crs
    if solvent_list is None:
        path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')
        solvent_list = [join(path, solv) for solv in os.listdir(path) if
                        solv not in ('__init__.py', 'README.rst')]
    solvent_list.sort()

    # Prepare the dataframe
    index = ['E_solv ' + i.rsplit('.', 1)[0].rsplit('/', 1)[-1] for i in solvent_list]
    index += ['Gamma ' + i.rsplit('.', 1)[0].rsplit('/', 1)[-1] for i in solvent_list]
    df = pd.DataFrame(index=index)

    # Run COSMO-RS
    init(path=mol_list[0].properties.path, folder='ligand_solvation')
    for mol in mol_list:
        coskf = get_surface_charge(mol, job=job_recipe.job1, s=job_recipe.s1)
        E_solv = get_solv(mol, solvent_list, coskf, job=job_recipe.job2, s=job_recipe.s2)
        df[mol.properties.name] = E_solv
    finish()

    # Update the database
    if arg.properties.database.write:
        pass


def get_surface_charge(mol, job=None, s=None):
    """ Construct the COSMO surface of the *mol*. """
    # Special procedure for ADF jobs
    # Use the gas-phase electronic structure as a fragment for the COSMO single point
    if job is ADFJob:
        s = get_surface_charge_adf(mol, job, s)

    results = mol.job_single_point(job, s, ret_results=True)
    results.wait()
    return results[get_coskf(results)]


def get_solv(mol, solvent_list, coskf, job=None, s=None):
    """ Calculate the solvation energy of *mol* in various *solvents*. """
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
        E_solv.append(result.get_energy())
        Gamma.append(result.get_activity_coefficient())

    # Return the solvation energies and activity coefficients as dict
    return np.array(E_solv + Gamma)


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
                return file
