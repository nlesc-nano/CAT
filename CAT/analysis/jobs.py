""" A module designed for running Jobs. """

__all__ = ['ams_job_mopac_crs', 'ams_job_mopac_opt', 'ams_job_mopac_sp', 'ams_job_uff_opt']

import shutil
from os.path import (dirname, join)

from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units
from scm.plams.tools.kftools import KFFile
from scm.plams.core.settings import Settings
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.functions import (init, finish, add_to_class)

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob

from qmflows.templates.templates import get_template

from .crs import CRSJob
from .thermo_chem import get_thermo
from ..misc import get_time
from ..qd_functions import (adf_connectivity, fix_h, fix_carboxyl, from_mol_other)
from ..data_handling.mol_export import export_mol


def ams_job_mopac_sp(mol):
    """
    Runs a gas-phase MOPAC single point.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    # Run MOPAC
    s = Settings()
    s.input = get_template('singlepoint.json')['specific']['ams']
    s.update(get_template('qd.json')['COSMO-MOPAC'])
    job = AMSJob(molecule=mol, settings=s, name='MOPAC')
    results = job.run()
    results.wait()
    mol.properties.energy.E = KFFile(results['mopac.rkf']).read('AMSResults', 'Energy')

    return mol


def ams_job_mopac_opt(mol):
    """
    Runs a gas-phase MOPAC geometry optimization.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    # Run MOPAC
    s = get_template('qd.json')['MOPAC geometry optimization']
    job = AMSJob(molecule=mol, settings=s, name='MOPAC')
    results = job.run()
    results.wait()
    mol.properties.energy.E = KFFile(results['mopac.rkf']).read('AMSResults', 'Energy')

    return mol


def ams_job_mopac_crs(mol):
    """
    Runs a MOPAC + COSMO-RS single point.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    path = mol.properties.path
    angstrom = Units.conversion_ratio('Bohr', 'Angstrom')
    solvents = ('Acetone', 'Acetonitrile', 'DMF', 'DMSO', 'Ethanol',
                'EtOAc', 'Hexane', 'Toluene', 'Water')
    solv_path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')

    # Run MOPAC
    init(path=path, folder='ligand')
    s1 = get_template('qd.json')['COSMO-MOPAC single point']
    s1.input.ams.System.Charge = mol.properties.charge
    mopac_job = AMSJob(molecule=mol, settings=s1, name='MOPAC')
    mopac_results = mopac_job.run()

    if 'mopac.coskf' in mopac_results.files:
        # Extract properties from mopac.coskf
        coskf = mopac_results['mopac.coskf']
        mol.properties.surface = KFFile(coskf).read('COSMO', 'Area') * angstrom**2
        mol.properties.volume = KFFile(coskf).read('COSMO', 'Volume') * angstrom**3

        # Prepare COSMO-RS a list of settings; one for each solvent
        parallel = JobRunner(parallel=True)
        s2 = get_template('qd.json')['COSMO-RS activity coefficient']
        s2.input = get_template('crs.json')['MOPAC PM6']
        s2.input.Compound._h = '"' + coskf + '"'
        s2_list = []
        for solv in solvents:
            s2_tmp = s2.copy()
            s2_tmp.input.compound._h = '"' + join(solv_path, solv + '.coskf') + '"'
            s2_list.append(s2_tmp)

        # Run COSMO-RS in parallel
        crs_jobs = [CRSJob(settings=s2, name='COSMO-RS_'+solv) for s2 in s2_list]
        crs_result = [job.run(jobrunner=parallel) for job in crs_jobs]

        # Extract properties from COSMO-RS_solv.crskf
        solvation_energy = Settings()
        for results, solv in zip(crs_result, solvents):
            results.wait()
            if 'COSMO-RS_' + solv + '.crskf' in results.files:
                solvation_energy[solv] = results.get_solvation_energy('$JN.crskf')
            else:
                solvation_energy[solv] = None
        mol.properties.energy = solvation_energy
    finish()

    if len(mol.properties.solvation) == len(solvents):
        shutil.rmtree(mopac_job.path.rsplit('/', 1)[0])

    return mol


def ams_job_uff_opt(mol, maxiter=1000, get_freq=False, fix_angle=True):
    """
    Runs an AMS UFF constrained geometry optimization.
    mol <plams.Molecule>: The input PLAMS molecule with the 'path', 'name' and
        'mol_indices' properties.
    bonds <list>[<list>[<int>, <float>]]: A nested list of atomic indices and bond orders.
    maxiter <int>: The maximum number of iterations during the geometry optimization.
    get_freq <bool>: Perform a frequency analyses after the geometry optimization.
    return <plams.Molecule>: A PLAMS molecule.
    """
    name = mol.properties.name + '.opt'
    path = mol.properties.path
    constrains = mol.properties.indices

    # Prepare the job settings
    init(path=path, folder='Quantum_dot')
    s = Settings()
    s.input = get_template('geometry.json')['specific']['ams']
    s.update(get_template('qd.json')['UFF'])
    s.input.ams.constraints.atom = constrains
    s.input.ams.system.bondorders._1 = adf_connectivity(mol)
    s.input.ams.geometryoptimization.maxiterations = int(maxiter / 2)
    if get_freq:
        s.input.ams.Properties.NormalModes = 'Yes'

    # Run the job (pre-optimization)
    job = AMSJob(molecule=mol, settings=s, name='UFF_part1')
    results = job.run()
    mol.from_mol_other(results.get_main_molecule())
    if get_freq:
        results.wait()
        mol.properties.energy = Settings(results.get_thermo())

    # Fix all O=C-O and H-C=C angles and continue the job
    if fix_angle:
        job = AMSJob(molecule=fix_carboxyl(fix_h(mol)), settings=s, name='UFF_part2')
        results = job.run()
        output_mol = results.get_main_molecule()
        mol.from_mol_other(output_mol)
    finish()

    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], join(path, name + '.ams.rkf'))
    shutil.copy2(results['uff.rkf'], join(path, name + '.uff.rkf'))
    shutil.copy2(results['$JN.out'], join(path, name + '.out'))
    shutil.rmtree(job.path.rsplit('/', 1)[0])

    # Write the reuslts to an .xyz and .pdb file
    mol.properties.name += '.opt'
    export_mol(mol, message='Optimized core + ligands:\t\t')
    mol.properties.name = mol.properties.name.split('.opt')[0]

    return mol


def type_to_string(job):
    """ Turn a <type> object into a <str> object. """
    job_dict = {ADFJob: 'adf', AMSJob: 'ams', DiracJob: 'dirac',
                Cp2kJob: 'cp2k', GamessJob: 'gamess', ORCAJob: 'orca'}
    try:
        return job_dict[job]
    except KeyError:
        print(get_time() + 'No default settings available for ' + str(job))


@add_to_class(Molecule)
def job_crs(self, job, settings, name='Single_point'):
    """ Function for running an COSMO-RS with <CRSJob>, extracting solvation energies.

    mol <plams.Molecule>: A PLAMS molecule.
    job <type>: A type object of a class derived from <Job>, e.g. AMSJob or Cp2kJob.
    settings <Settings>: The settings for *job*.
    name <str>: The name of *job*.
    """
    # Grab the default settings for a specific job and update them with user provided settings
    s = Settings()
    s.input = get_template('qd.json')['MOPAC PM6']
    s.ignore_molecule = True
    s.update(settings)

    # Run the job; extract energies
    my_job = job(settings=s, name=name)
    results = my_job.run()
    results.wait()
    self.properties.energy.E = results.get_energy(unit='kcal/mol')


@add_to_class(Molecule)
def job_single_point(self, job, settings, name='Single_point'):
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


@add_to_class(Molecule)
def job_geometry_opt(self, job, settings, name='Geometry_optimization'):
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


@add_to_class(Molecule)
def job_freq(self, job, settings, name='Frequency_analysis', opt=True):
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
