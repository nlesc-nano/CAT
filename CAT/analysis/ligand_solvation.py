""" A module designed for calculating solvation energies. """

__all__ = ['init_solv']

import os
from os.path import (join, dirname)

from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob

from qmflows.templates.templates import get_template

from .crs import CRSJob


def init_solv(mol, solvent_list=None, job1=None, s1=None, job2=None, s2=None):
    """ Initialize the solvation energy calculation. """
    if solvent_list is None:
        path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')
        solvent_list = [join(path, solv) for solv in os.listdir(path) if not
                        '__init__.py' or not 'README.rst']

    coskf = get_surface_charge(mol, job=job1, s=s1)
    solv_dict = get_solv(mol, solvent_list, coskf, job=job2, s=s2)
    mol.properties.energy.E_solv = solv_dict


def get_surface_charge(mol, job=None, s=None):
    """ Construct the COSMO surface of the *mol*. """
    init(path=mol.properties.path, folder='coskf')

    # Switch to default settings if no job & s are <None>
    if job is None and s is None:
        job = AMSJob
        s = get_template('qd.json')['COSMO-MOPAC']
    elif job is None or s is None:
        finish()
        raise TypeError('job & s should neither or both be None')

    # Special procedure for ADF jobs
    # Use the gas-phase electronic structure as a fragment for the COSMO single point
    if job is ADFJob:
        s = get_surface_charge_adf(mol, job, s)

    results = mol.job_single_point(job, s, ret_results=True)
    finish()

    return results[get_coskf(results)]


def get_solv(mol, solvent_list, coskf, job=None, s=None):
    """ Calculate the solvation energy of *mol* in various *solvents*. """
    init(path=mol.properties.path, folder='crs')

    # Switch to default settings if no job & s are <None>
    if job is None and s is None:
        job = CRSJob
        s = get_template('qd.json')['COSMO-RS activity coefficient']
        s.input.update(get_template('crs.json')['MOPAC PM6'])
    elif job is None or s is None:
        finish()
        raise TypeError('job & s should neither or both be None')

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

    # Extract solvation energies
    solv_dict = {}
    for result in results:
        result.wait()
        solv_dict[result.job.settings.name] = result.get_energy()
    finish()

    return solv_dict


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
    s.update(get_template('qd.json')['COSMO-ADF'])
    s.input.fragments.gas = coskf

    return s


def get_coskf(results, extensions=['.coskf', '.t21']):
    """ Return the file in results containing the COSMO surface. """
    for file in results.files:
        for ext in extensions:
            if ext in file:
                return file
