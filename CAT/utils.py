"""A module with miscellaneous functions."""

import os
import time
import yaml
import pkg_resources as pkg
from os.path import join
from typing import Callable

from scm.plams.core.settings import Settings

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob

__all__ = ['check_sys_var', 'dict_concatenate', 'get_time', 'get_template']

_job_dict = {
    ADFJob: 'adf',
    AMSJob: 'ams',
    DiracJob: 'dirac',
    Cp2kJob: 'cp2k',
    GamessJob: 'gamess',
    ORCAJob: 'orca'
}


def type_to_string(job: Callable) -> str:
    """Turn a :class:`type` instance into a string."""
    try:
        return _job_dict[job]
    except KeyError:
        print(get_time() + 'WARNING: No default settings available for ' + str(job))
        return ''


def get_time() -> str:
    """Return the current time as string."""
    return '[{}] '.format(time.strftime('%H:%M:%S'))


def check_sys_var() -> None:
    """Validate all ADF environment variables.

    Raises
    ------
    EnvironmentError
        Raised if one or more of the following environment variables are absent:
        * ``'ADFBIN'``
        * ``'ADFHOME'``
        * ``'ADFRESOURCES'``
        * ``'SCMLICENSE'``

    ImportError
        Raised if an ADF version prior to 2019 is found.

    """
    sys_var = ['ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE']
    sys_var_exists = [item in os.environ for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            err = 'WARNING: The environment variable {} has not been set'
            print(get_time() + err.format(sys_var[i]))

    if not all(sys_var_exists):
        raise EnvironmentError(get_time() + 'One or more ADF environment variables have '
                               'not been set, aborting ADF job.')

    if '2019' not in os.environ['ADFHOME']:
        error = get_time() + 'No ADF/2019 detected in ' + os.environ['ADFHOME']
        error += ', aborting ADF job.'
        raise ImportError(error)


def dict_concatenate(dic: dict) -> dict:
    """Concatenates a list of dictionaries."""
    ret = {}
    for item in dic:
        ret.update(item)
    return ret


def get_template(template_name: str,
                 from_cat_data: bool = True) -> Settings:
    """Grab a yaml template and return it as Settings object."""
    if from_cat_data:
        path = join('data/templates', template_name)
        xs = pkg.resource_string('CAT', path)
        return Settings(yaml.load(xs.decode(), Loader=yaml.FullLoader))
    else:
        with open(template_name, 'r') as file:
            return Settings(yaml.load(file, Loader=yaml.FullLoader))
