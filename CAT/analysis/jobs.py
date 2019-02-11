""" A module designed for running Jobs. """

__all__ = ['job_single_point', 'job_geometry_opt', 'job_freq']

from scm.plams.mol.molecule import Molecule
from scm.plams.core.settings import Settings
from scm.plams.core.functions import add_to_class

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob

from qmflows.templates.templates import get_template

from .thermo_chem import get_thermo
from ..misc import get_time
from ..qd_functions import (adf_connectivity, from_mol_other)


def type_to_string(job):
    """ Turn a <type> object into a <str> object. """
    job_dict = {ADFJob: 'adf', AMSJob: 'ams', DiracJob: 'dirac',
                Cp2kJob: 'cp2k', GamessJob: 'gamess', ORCAJob: 'orca'}
    try:
        return job_dict[job]
    except KeyError:
        print(get_time() + 'No default settings available for ' + str(job))


@add_to_class(Molecule)
def job_single_point(self, job, settings, name='Single_point', ret_results=False):
    """ Function for running an arbritrary <Job>, extracting total energies.

    mol <plams.Molecule>: A PLAMS molecule.
    job <type>: A type object of a class derived from <Job>, e.g. AMSJob or Cp2kJob.
    settings <Settings>: The settings for *job*.
    name <str>: The name of *job*.
    """
    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = get_template('singlepoint.json')['specific'][type_to_string(job)]
    s.update(settings)
    if job == AMSJob:
        s.input.ams.system.bondorders._1 = adf_connectivity(self)

    # Run the job; extract energies
    my_job = job(molecule=self, settings=s, name=name)
    results = my_job.run()
    results.wait()
    self.properties.energy.E = results.get_energy(unit='kcal/mol')

    # Return results
    if ret_results:
        return results


@add_to_class(Molecule)
def job_geometry_opt(self, job, settings, name='Geometry_optimization', ret_results=False):
    """ Function for running an arbritrary <Job>, extracting total energies and final geometries.

    mol <plams.Molecule>: A PLAMS molecule.
    job <type>: A type object of a class derived from <Job>, e.g. AMSJob or Cp2kJob.
    settings <Settings>: The settings for *job*.
    name <str>: The name of *job*.
    """
    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = get_template('geometry.json')['specific'][type_to_string(job)]
    s.update(settings)
    if job == AMSJob:
        s.input.ams.system.bondorders._1 = adf_connectivity(self)

    # Run the job; extract geometries and energies
    my_job = job(molecule=self, settings=s, name=name)
    results = my_job.run()
    results.wait()
    self.properties.energy.E = results.get_energy(unit='kcal/mol')
    self.from_mol_other(results.get_main_molecule())

    # Return results
    if ret_results:
        return results


@add_to_class(Molecule)
def job_freq(self, job, settings, name='Frequency_analysis', opt=True, ret_results=False):
    """ Function for running an arbritrary <Job>, extracting total energies, final geometries and
    thermochemical quantities derived from vibrational frequencies.

    mol <plams.Molecule>: A PLAMS molecule.
    job <type>: A type object of a class derived from <Job>, e.g. AMSJob or Cp2kJob.
    settings <Settings>: The settings for *job*.
    name <str>: The name of *job*.
    opt <bool>: Preceed the frequency analysis with a geometry optimization.
    """
    # Preceed the frequency analysis with a geometry optimization
    if opt:
        self.job_geometry_opt(job, settings)

    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = get_template('freq.json')['specific'][type_to_string(job)]
    s.update(settings)
    if job == AMSJob:
        s.input.ams.system.bondorders._1 = adf_connectivity(self)

    # Run the job; extract geometries and (Gibbs free) energies
    my_job = job(molecule=self, settings=s, name=name)
    results = my_job.run()
    results.wait()
    self.properties.frequencies = results.get_frequencies()
    self.properties.energy = get_thermo(self,
                                        self.properties.frequencies,
                                        results.get_energy(unit='kcal/mol'))

    # Return results
    if ret_results:
        return results
