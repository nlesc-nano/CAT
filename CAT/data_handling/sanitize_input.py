""" A module designed for sanitizing and interpreting the input file. """

__all__ = ['sanitize_arg_dict']

from itertools import chain
from schema import (Schema, Or, And, Use)

from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.uff import UFFJob
from scm.plams.interfaces.adfsuite.band import BANDJob
from scm.plams.interfaces.adfsuite.dftb import DFTBJob
from scm.plams.interfaces.adfsuite.mopac import MOPACJob
from scm.plams.interfaces.adfsuite.reaxff import ReaxFFJob

from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob
from scm.plams.interfaces.thirdparty.dftbplus import DFTBPlusJob

from scm.plams.core.basejob import Job
from scm.plams.core.settings import Settings
from scm.plams.tools.periodic_table import PeriodicTable

from qmflows.templates.templates import get_template

from ..analysis.crs import CRSJob
from ..misc import get_time
from ..qd_functions import to_atnum


str_to_class = {'adf': ADFJob, 'adfjob': ADFJob,
               'ams': AMSJob, 'amsjob': AMSJob,
               'uff': UFFJob, 'uffjob': UFFJob,
               'band': BANDJob, 'bandjob': BANDJob,
               'dftb': DFTBJob, 'dftbjob': DFTBJob,
               'mopac': MOPACJob, 'mopacjob': MOPACJob,
               'reaxff': ReaxFFJob, 'reaxffjob': ReaxFFJob,

               'cp2k': Cp2kJob, 'cp2kjob': Cp2kJob,
               'orca': ORCAJob, 'orca': ORCAJob,
               'dirac': DiracJob, 'diracjob': DiracJob,
               'gamess': GamessJob, 'gamessjob': GamessJob,
               'dftbplus': DFTBPlusJob, 'dftbplusjob': DFTBPlusJob,
               'crs': CRSJob, 'cosmo-rs': CRSJob, 'crsjob': CRSJob}


def init_mol(mol_list):
    """ Validate the input of *mol_list*; returns a <Settings> object with validated input. """
    ret = Settings()
    return ret


def get_default_arg_dict():
    return Settings({'dir_name_list': ('core', 'ligand', 'QD'),
                     'dummy': 'Cl',
                     'use_database': True,
                     'ligand_opt': True,
                     'ligand_crs': False,
                     'qd_opt': False,
                     'qd_int': False,
                     'qd_dissociate': False,
                     'split': True,
                     })


def sanitize_arg_dict(arg_dict):
    """ Validate the input of *argument_dict*; returns a <Settings> object with validated input. """
    ret = Settings()
    argument_dict = get_default_arg_dict()
    argument_dict.update(arg_dict)

    # Validate arguments consisting of booleans, integers, strings and/or iterable
    ret.dir_name_list = val_dir_name_list(argument_dict['dir_name_list'])
    ret.dummy = val_atnum(argument_dict['dummy'])
    ret.use_database = val_bool(argument_dict['use_database'])
    ret.ligand_opt = val_bool(argument_dict['ligand_opt'])
    ret.split = val_bool(argument_dict['split'])
    ret.qd_int = val_bool(argument_dict['qd_int'])

    # Validate arguments containing job recipes
    s_crs = get_template('qd.json')['COSMO-RS activity coefficient']
    s_crs.update(get_template('crs.json')['MOPAC PM6'])
    ret.ligand_crs = val_job(argument_dict['ligand_crs'],
                             job1=AMSJob,
                             job2=CRSJob,
                             s1=get_template('qd.json')['COSMO-MOPAC'],
                             s2=s_crs)

    ret.qd_opt = val_job(argument_dict['qd_opt'],
                         job1=AMSJob,
                         job2=AMSJob,
                         s1=get_template('qd.json')['UFF'],
                         s2=get_template('qd.json')['UFF'])

    ret.qd_dissociate = val_job(argument_dict['qd_dissociate'],
                                job1=AMSJob,
                                job2=AMSJob,
                                s1=get_template('qd.json')['MOPAC'],
                                s2=get_template('qd.json')['UFF'])

    return ret


def val_dir_name_list(dir_name_list):
    """ Validate an iterable of length 3; returns a <tuple> consisting of 3 <str>. """
    schema = Schema(And([str, str, str], Use(tuple)))
    return schema.validate(list(dir_name_list))


def val_atnum(atnum):
    """ Validate an atomic number or symbol; returns an atomic number <int>. """
    at_gen = chain.from_iterable([[i, j[0]] for i, j in enumerate(PeriodicTable.data)])
    schema = Schema(And(Or(int, str), lambda n: n in at_gen, Use(to_atnum)))
    return schema.validate(atnum)


def val_bool(my_bool):
    """ Validate a boolean; returns a <bool>. """
    return Schema(bool).validate(my_bool)


def val_job(job, job1=None, job2=None, s1=None, s2=None):
    """ Validate a job recipe.
    Returns a dictionary: {'job1': <Job>, 'job2': <Job>, 's1': <Settings>, 's2': <Settings>}. """
    # Validate the object type
    Schema(Or(bool, dict)).is_valid(job)
    if isinstance(job, bool):
        if job is False:
            return job
        job = {'job1': None, 's1': None,
               'job2': None, 's2': None}

    # Validate the object types of the various elements
    schema = Schema({'job1': Or(None, Job, str),
                     's1': Or(None, dict),
                     'job2': Or(None, Job, str),
                     's2': Or(None, dict)})
    schema.is_valid(job)

    # Assign proper default settings
    str_to_def = {'job1': job1, 'job2': job2, 's1': s1, 's2': s2}
    for key in job:
        if job[key] is None:
            job[key] = str_to_def[key]
        elif isinstance(job[key], str):
            try:
                job[key] = str_to_class[job[key].lower()]
            except KeyError:
                raise KeyError(get_time() + 'No Job-derived object exists for the string:', job[key]
                               + ', please provide the actual <Job> object instead of <str>')
        elif isinstance(job[key], Job):
            pass
        elif isinstance(job[key], dict):
            job[key] = Settings(job[key])
        else:
            raise TypeError(get_time() + str(type(job[key])), 'is an unspported object type')
    return job
