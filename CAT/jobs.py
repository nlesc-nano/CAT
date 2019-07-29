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

from shutil import rmtree
from typing import (Optional, Callable)
from os.path import join

import numpy as np

from scm.plams import (Molecule, Settings, Results, config, add_to_class, ResultsError)
from scm.plams.core.basejob import Job
from scm.plams.tools.units import Units

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kResults

import qmflows

from .logger import (logger, log_start, log_succes, log_fail, log_copy)
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


def pre_process_settings(mol: Molecule,
                         s: Settings,
                         job_type: type,
                         template_name: str) -> Settings:
    """Update all :class:`Settings`, **s**, with those from a QMFlows template (see **job**)."""
    ret = Settings()
    ret.input = getattr(qmflows, template_name)['specific'][type_to_string(job_type)].copy()
    ret.update(s)
    if job_type == AMSJob:
        ret.input.ams.system.bondorders._1 = adf_connectivity(mol)
        if 'uff' not in s.input:
            ret.input.ams.system.charge = sum([at.properties.charge for at in mol if
                                               'charge' in at.properties])
    return ret


def retrieve_results(results: Results,
                     job_preset: str) -> None:
    """Unpack the :class:`results` from a PLAMS-facilitated calculation.

    Performs an inplace update of **results**:

        * The Cartesian coordinates of :attr:`Results.job.molecule` are updated.
        * Energies and frequencies are added to :attr:`Results.job.molecule.properties`.
        *

    Paramaters
    ----------
    results : |plams.Results|_
        A PLAMS :class:`Results` instance.

    job_preset : str
        The name of a job preset from :mod:`CAT.jobs`.
        Accepted values are ``"single point"``, ``"geometry optimization"`` or
        ``"frequency analysis"``.

    """
    if job_preset not in ('geometry optimization', 'frequency analysis', 'single point'):
        raise ValueError(f'Invalid value for job_preset: {repr(job_preset)}')

    # Unpack arguments
    job = results.job
    name = job.name
    mol = job.molecule

    # Define more arguments
    freq = np.empty(1)
    nan_dict = {
        'frequencies': np.nan,
        'energy': {'E': np.nan, 'H': np.nan, 'S': np.nan, 'G': np.nan}
    }

    try:
        # Read all relevant results
        energy = mol.properties.energy.E = results.get_energy(unit='kcal/mol')
        if job_preset in ('geometry optimization', 'frequency analysis'):
            mol.from_mol_other(results.get_main_molecule())
        if job_preset == 'frequency analysis':
            freq = mol.properties.frequencies = results.get_frequencies()
            energy = mol.properties.energy = get_thermo(mol, freq, energy)

        # Evaluate all results
        if not (energy and isinstance(freq, np.ndarray)):
            raise ResultsError(f'Failed to retrieve results of {name}')
        log_succes(job, mol, job_preset, name)

    except Exception as ex:  # Failed to retrieve results
        mol.properties.soft_update(nan_dict)
        log_fail(job, mol, job_preset, name)
        logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)

    else:
        if job.status != 'copied':
            return None

        # results.job is a copy from a previously run job.
        # Exchange the attributes of results.job wth those of the revious job
        job_old = config.default_jobmanager.hashes[job._hash]()
        log_copy(job, mol, job_preset, name, job_old)
        rmtree(job.path)
        for key, value in vars(job_old).items():
            if key != 'molecule':
                setattr(job, key, value)
    return None


@add_to_class(Job)
def _finalize(self):
    """Modified PLAMS :meth:`Job._finalize` method.

    1. All references to the :func:`scm.plams.log` function are removed.
    2. A $JN.hash file is created which contains the hash of this instance.

    Gather the results of the job execution and organize them.
    This method collects steps 9-12 from :ref:`job-life-cycle`.
    Should not be overridden.

    """
    if config.preview is False:
        self.results.collect()
        self.results.finished.set()
        if self.status != 'crashed':
            self.status = 'finished'
            self._log_status(3)
            if self.check():
                self.results._clean(self.settings.keep)
                self.postrun()
                self.status = 'successful'
                if self.settings.pickle:
                    self.pickle()
                    filename = join(self.path, self.name + '.hash')
                    with open(filename, 'w') as f:
                        f.write(self.hash())
            else:
                self.status = 'failed'
    else:
        self.status = 'preview'
        self.results.finished.set()
    self.results.done.set()

    if self.parent and self in self.parent:
        self.parent._notify()

    self._log_status(1)


@add_to_class(Molecule)
def job_single_point(self, job_type: Callable,
                     settings: Settings,
                     name: str = 'Single_point',
                     ret_results: bool = False) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies.

    Paramaters
    ----------
    job_type : |Callable|_
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
    s = pre_process_settings(self, settings, job_type, 'singlepoint')

    # Run the job; extract energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'single point', _name)

    results = job.run()
    retrieve_results(results, 'single point')

    inp_name = join(job.path, job.name + '.in')
    self.properties.job_path.append(inp_name)

    # Return results
    if ret_results:
        return results
    return None


@add_to_class(Molecule)
def job_geometry_opt(self, job_type: Callable,
                     settings: Settings,
                     name: str = 'Geometry_optimization',
                     ret_results: bool = False) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies and final geometries.

    Paramaters
    ----------
    job_type : |Callable|_
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
    s = pre_process_settings(self, settings, job_type, 'geometry')

    # Run the job; extract geometries and energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'geometry optimization', _name)

    results = job.run()
    retrieve_results(results, 'geometry optimization')

    inp_name = join(job.path, job.name + '.in')
    self.properties.job_path.append(inp_name)

    # Return results
    if ret_results:
        return results
    return None


@add_to_class(Molecule)
def job_freq(self, job_type: Callable,
             settings: Settings,
             name: str = 'Frequency_analysis',
             opt: bool = True,
             ret_results: bool = False) -> Optional[Results]:
    """Function for running an arbritrary Jobs

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
        self.job_geometry_opt(job_type, settings)

    # Grab the default settings for a specific job and update them with user provided settings
    s = pre_process_settings(self, settings, job_type, 'freq')

    # Run the job; extract geometries and (Gibbs free) energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'frequency analysis', _name)

    results = job.run()
    retrieve_results(results, 'frequency analysis')

    inp_name = join(job.path, _name + '.in')
    self.properties.job_path.append(inp_name)

    # Return results
    if ret_results:
        return results
    return None
