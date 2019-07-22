"""
CAT.jobs
========

A module designed for running Jobs.

Index
-----
.. currentmodule:: CAT.jobs
.. autosummary::
    get_main_molecule
    get_energy
    job_single_point
    job_geometry_opt
    job_freq

API
---
.. automethod:: get_main_molecule
.. automethod:: get_energy
.. autofunction:: job_single_point
.. autofunction:: job_geometry_opt
.. autofunction:: job_freq

"""

from os.path import join
from typing import (Optional, Callable)

import numpy as np

from scm.plams import (Molecule, Settings, Results, config)
from scm.plams.core.functions import add_to_class
from scm.plams.tools.units import Units

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kResults

import qmflows

from .logger import logger
from .thermo_chem import get_thermo
from .utils import type_to_string
from .mol_utils import (adf_connectivity, from_mol_other)

__all__ = ['job_single_point', 'job_geometry_opt', 'job_freq']


@add_to_class(Cp2kResults)
def get_main_molecule(self) -> Optional[Molecule]:
    for file in self.files:
        if '.xyz' in file:
            return Molecule(join(self.job.path, file))
    return None


@add_to_class(Cp2kResults)
def get_energy(self, index: int = 0,
               unit: str = 'Hartree') -> float:
    """Returns last occurence of 'Total energy:' in the output."""
    energy = self._get_energy_type('Total', index=index)
    return Units.convert(energy, 'Hartree', unit)


def _get_name(name) -> str:
    manager = config.default_jobmanager
    if name in manager.names:
        return name + '.' + str(1 + manager.names[name]).zfill(manager.settings.counter_len)
    else:
        return name


@add_to_class(Molecule)
def job_single_point(self, job: Callable,
                     settings: Settings,
                     name: str = 'Single_point',
                     ret_results: bool = False) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies.

    Paramaters
    ----------
    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = qmflows.singlepoint['specific'][type_to_string(job)].copy()
    s.update(settings)
    if job == AMSJob:
        s.input.ams.system.bondorders._1 = adf_connectivity(self)
        if 'uff' not in s.input:
            s.input.ams.system.charge = sum([at.properties.charge for at in self if
                                             'charge' in at.properties])

    # Run the job; extract energies
    my_job = job(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    logger.info(f'{my_job.__class__.__name__}: {self.properties.name} single point ({_name}) '
                'has started')
    results = my_job.run()

    try:
        self.properties.energy.E = energy = results.get_energy(unit='kcal/mol')
        if not energy:
            raise TypeError
        logger.info(f'{results.job.__class__.__name__}: {self.properties.name} single point '
                    f'({_name}) is successful')
    except TypeError:
        self.properties.energy.E = np.nan
        logger.error(f'{results.job.__class__.__name__}: {self.properties.name} single point '
                     f'({_name}) has failed')

    inp_name = join(my_job.path, my_job.name + '.in')
    self.properties.job_path.append(inp_name)

    # Return results
    if ret_results:
        return results
    return None


@add_to_class(Molecule)
def job_geometry_opt(self, job: Callable,
                     settings: Settings,
                     name: str = 'Geometry_optimization',
                     ret_results: bool = False) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies and final geometries.

    Paramaters
    ----------
    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = qmflows.geometry['specific'][type_to_string(job)].copy()
    s.update(settings)
    if job == AMSJob:
        s.input.ams.system.bondorders._1 = adf_connectivity(self)
        if 'uff' not in s.input:
            s.input.ams.system.charge = sum([at.properties.charge for at in self if
                                             'charge' in at.properties])

    # Run the job; extract geometries and energies
    my_job = job(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    logger.info(f'{my_job.__class__.__name__}: {self.properties.name} optimization ({_name})'
                ' has started')
    results = my_job.run()

    try:
        self.from_mol_other(results.get_main_molecule())
        self.properties.energy.E = energy = results.get_energy(unit='kcal/mol')
        if not energy:
            raise TypeError
        logger.info(f'{results.job.__class__.__name__}: {self.properties.name} optimization '
                    f'({_name}) is successful')
    except TypeError:
        self.properties.energy.E = np.nan
        logger.error(f'{results.job.__class__.__name__}: {self.properties.name} optimization '
                     f'({_name}) has failed')

    inp_name = join(my_job.path, my_job.name + '.in')
    self.properties.job_path.append(inp_name)

    # Return results
    if ret_results:
        return results
    return None


@add_to_class(Molecule)
def job_freq(self, job, settings, name='Frequency_analysis', opt=True, ret_results=False):
    """ Function for running an arbritrary Jobs

    Extracts total energies, final geometries and
    thermochemical quantities derived from vibrational frequencies.

    Paramaters
    ----------
    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    opt : bool
        Perform a geometry optimization (see :func:`.job_geometry_opt`) before calculating
        frequencies.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Preceed the frequency analysis with a geometry optimization
    if opt:
        self.job_geometry_opt(job, settings)

    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = qmflows.freq['specific'][type_to_string(job)].copy()
    s.update(settings)
    if job == AMSJob:
        s.input.ams.system.bondorders._1 = adf_connectivity(self)
        if 'uff' not in s.input:
            s.input.ams.system.charge = sum([at.properties.charge for at in self if
                                             'charge' in at.properties])

    # Run the job; extract geometries and (Gibbs free) energies
    my_job = job(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    logger.info(f'{my_job.__class__.__name__}: {self.properties.name} frequency analysis ({_name}) '
                'has started')
    results = my_job.run()

    try:
        self.properties.frequencies = freq = results.get_frequencies()
        self.properties.energy = energy = get_thermo(self, self.properties.frequencies,
                                                     results.get_energy(unit='kcal/mol'))
        if not isinstance(freq, np.ndarray) or energy:
            raise TypeError
        logger.info(f'{results.job.__class__.__name__}: {self.properties.name} frequency analysis '
                    f'({_name}) is successful')
    except TypeError:
        self.properties.frequencies = np.nan
        self.properties.energy = {'E': np.nan, 'H': np.nan, 'S': np.nan, 'G': np.nan}
        logger.error(f'{results.job.__class__.__name__}: {self.properties.name} frequency analysis '
                     f'({_name}) has failed')

    inp_name = join(my_job.path, _name + '.in')
    self.properties.job_path.append(inp_name)

    # Return results
    if ret_results:
        return results
    return None
