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
import sys
import yaml
import pkg_resources as pkg
from shutil import rmtree
from typing import (Callable, Iterable, Optional, Union)
from os.path import (join, isdir, isfile, exists)

from scm.plams import (config, Settings, Molecule, MoleculeError, PeriodicTable, init)
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob

from .logger import logger
from .mol_utils import to_atnum
from .gen_job_manager import GenJobManager

__all__ = ['check_sys_var', 'dict_concatenate', 'get_template']

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
        logger.warn(f"No default settings available for type: '{job.__class__.__name__}'")
        return ''


def check_sys_var() -> None:
    """Validate all ADF environment variables.

    Raises
    ------
    EnvironmentError
        Raised if an ADF version prior to 2019 is found or if one or more of the
        following environment variables are absent:
        * ``'ADFBIN'``
        * ``'ADFHOME'``
        * ``'ADFRESOURCES'``
        * ``'SCMLICENSE'``

    """
    sys_var = ('ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE')
    sys_var_exists = [item in os.environ and os.environ[item] for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            logger.warn(f"The environment variable '{sys_var[i]}' has not been set")

    if not all(sys_var_exists):
        raise EnvironmentError('One or more ADF environment variables have not been set, '
                               'aborting ADF job')

    if '2019' not in os.environ['ADFHOME']:
        raise EnvironmentError(f"ADF/2019 not detected in {os.environ['ADFHOME']}")


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
        raise FileNotFoundError(f"'{path}' not found")
    elif isfile(path):
        raise NotADirectoryError(f"'{path}' is not a directory")


def validate_core_atom(atom: Union[str, int]) -> Union[Molecule, int]:
    """Parse and validate the ``["optional"]["qd"]["dissociate"]["core_atom"]`` argument."""
    # Potential atomic number or symbol
    if isinstance(atom, int) or atom in PeriodicTable.symtonum:
        return to_atnum(atom)

    # Potential SMILES string
    try:
        mol = from_smiles(atom)
    except Exception as ex:
        raise ex.__class__(f'Failed to recognize {repr(atom)} as a valid atomic number, '
                           f'atomic symbol or SMILES string\n\n{ex}')

    # Double check the SMILES string:
    charge_dict = {}
    for at in mol:
        charge = at.properties.charge
        try:
            charge_dict[charge] += 1
        except KeyError:
            charge_dict[charge] = 1
    if 0 in charge_dict:
        del charge_dict[0]

    # Only a single charged atom is allowed
    if len(charge_dict) > 1:
        charge_count = sum([v for v in charge_dict.values()])
        raise MoleculeError(f'The SMILES string {repr(atom)} contains more than one charged atom: '
                            f'charged atom count: {charge_count}')
    return mol


class SupressPrint:
    """A context manager for supressing :func:`print` calls."""

    def __init__(self) -> None:
        """Initialize a :class:`SupressPrint` instance."""
        self.stdout = None

    def __enter__(self) -> None:
        """Enter the :class:`SupressPrint` context manager."""
        self.stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args) -> None:
        """Exit the :class:`SupressPrint` context manager."""
        sys.stdout.close()
        sys.stdout = self.stdout
        self.stdout = None


def restart_init(path: str,
                 folder: str,
                 hashing: Optional[str] = 'input') -> None:
    """A wrapper around the plams.init_ function; used for importing one or more previous jobs.

    All pickled .dill files in **path**/**folder**/ will be loaded into the
    :class:`GenJobManager` instance initiated by :func:`init`.

    .. _plams.init: https://www.scm.com/doc/plams/components/functions.html#scm.plams.core.functions.init

    Note
    ----
    Previous jobs are stored in a more generator-esque manner in :attr:`GenJobManager.hashes` and
    :class:`Job` instances are thus created on demand rather than
    permanently storing them in memory.

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

    # Create a job manager
    settings = Settings({'counter_len': 3, 'hashing': hashing, 'remove_empty_directories': True})
    manager = GenJobManager(settings, path, folder, hashing)

    # Change the default job manager
    with SupressPrint():
        init()
    rmtree(config.default_jobmanager.workdir)
    config.default_jobmanager = manager
    config.log.file = 3
    config.log.stdout = 0

    workdir = manager.workdir
    if not isdir(workdir):
        os.mkdir(workdir)
        return

    # Update the default job manager with previous Jobs
    for f in os.listdir(workdir):
        job_dir = join(workdir, f)
        if not isdir(job_dir):  # Not a directory; move along
            continue

        dill_file = join(job_dir, f + '.dill')
        if isfile(dill_file):  # Update JobManager.hashes
            manager.load_job(dill_file)

        # Grab the job name
        try:
            name, num = f.rsplit('.', 1)
            num = int(num)
        except ValueError:  # Jobname is not appended with a number
            name = f
            num = 1

        # Update JobManager.names
        try:
            manager.names[name] = max(manager.names[name], num)
        except KeyError:
            manager.names[name] = num
