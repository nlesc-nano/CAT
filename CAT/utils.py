"""
CAT.utils
=========

A module with miscellaneous functions.

Index
-----
.. currentmodule:: CAT.utils
.. autosummary::
    type_to_string
    get_time
    check_sys_var
    dict_concatenate
    get_template

API
---
.. autofunction:: type_to_string
.. autofunction:: get_time
.. autofunction:: check_sys_var
.. autofunction:: dict_concatenate
.. autofunction:: get_template

"""

import os
import time
import yaml
import pkg_resources as pkg
from shutil import rmtree
from typing import (Callable, Iterable, Optional)
from os.path import (join, isdir, isfile, exists)

from scm.plams import (JobManager, init, config, Settings)
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
    """Turn a :class:`type` instance into a :class:`str`."""
    try:
        return _job_dict[job]
    except KeyError:
        err = 'WARNING: No default settings available for {}'
        print(get_time() + err.format(repr(job.__class__.__name__)))
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
    sys_var = ('ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE')
    sys_var_exists = [item in os.environ and os.environ[item] for item in sys_var]
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


def dict_concatenate(dict_list: Iterable[dict]) -> dict:
    """Concatenates a list of dictionaries."""
    ret = {}
    for item in dict_list:
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


def validate_path(path: Optional[str]) -> str:
    """Validate a provided directory path.

    Parameters
    ----------
    path : str
        Optional: A path to a directory.
        Will default to the current working directory if ``None``.

    Results
    -------
    |str|_
        Returns either **path** or the current working directory.

    Raises
    ------
    FileNotFoundError
        Raised if **path** cannot be found.

    NotADirectoryError
        Raised if **path** is not a directory.

    """
    if path in (None, '.', ''):
        return os.getcwd()
    elif isdir(path):
        return path
    elif not exists(path):
        raise FileNotFoundError(get_time() + f"'{path}' not found")
    elif isfile(path):
        raise NotADirectoryError(get_time() + f"'{path}' is not a directory")


def restart_init(path: str,
                 folder: str,
                 hashing: Optional[str] = 'input') -> None:
    """A wrapper around the plams.init_ function; used for importing one or more previous jobs.

    All pickled .dill files in **path**/**folder**/ will be loaded into the :class:`JobManager` instance
    initiated by :func:`init`.

    .. _plams.init: https://www.scm.com/doc/plams/components/functions.html#scm.plams.core.functions.init

    Paramaters
    ----------
    path : str
        The path to the PLAMS workdir.

    folder : str
        The name of the PLAMS workdir.

    hashing : str
        Optional: The type of hashing used by the PLAMS :class:`JobManager`.
        Accepted values are: ``"input"``, ``"runscript"``, ``"input+runscript"`` and ``None``.

    """  # noqa
    attr_dict = {
        'path': path,
        'foldername': folder,
        'workdir': join(path, folder),
        'logfile': join(path, folder, 'logfile'),
        'input': join(path, folder, hashing)
    }

    # Create a job manager
    settings = Settings({'counter_len': 3, 'hashing': hashing, 'remove_empty_directories': True})
    manager = JobManager(settings)
    for k, v in attr_dict.items():
        setattr(manager, k, v)

    # Change the default job manager
    init()
    rmtree(config.default_jobmanager.workdir)
    config.default_jobmanager = manager

    # Update the default job manager with previous Jobs
    for folder in os.listdir(manager.workdir):
        dill_file = join(manager.workdir, folder, folder + '.dill')
        if not isfile(dill_file):  # Not a .dill file; move along
            continue

        # Update JobManager.hashes
        job = manager.load_job(dill_file)

        # Update JobManager.jobs
        manager.jobs.append(job)

        # Grab the job name
        name, num = folder.rsplit('.', 1)
        try:
            int(num)
        except ValueError:  # Jobname is not appended with a number
            name = folder

        # Update JobManager.names
        try:
            manager.names[name] += 1
        except KeyError:
            manager.names[name] = 1
