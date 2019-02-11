""" A module designed for optimizing the combined ligand & core. """

__all__ = ['init_qd_opt']

from scm.plams.core.settings import Settings
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob
import scm.plams.interfaces.molecule.rdkit as molkit

from qmflows.templates.templates import get_template

from ..qd_functions import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt
from ..data_handling.mol_export import export_mol


def qd_opt(mol, job1=None, job2=None, s1=None, s2=None):
    """ """
    if job1 is None and s1 is None:
        job1 = AMSJob
        s1 = get_template('qd.json')['UFF']
    elif job1 is None or s1 is None:
        finish()
        raise TypeError('job1 & s1 should neither or both be None')

    if job2 is None and s2 is None:
        job2 = AMSJob
        s2 = get_template('qd.json')['UFF']
    elif job2 is None or s2 is None:
        finish()
        raise TypeError('job2 & s2 should neither or both be None')

    # Prepare the job settings
    init(path=mol.properties.path, folder='QD_opt')
    s = Settings()
    s.input = get_template('geometry.json')['specific']['ams']
    s.update(get_template('qd.json')['UFF'])
    s.input.ams.constraints.atom = mol.properties.indices
    mol.job_geometry_opt(job1, s1, name='QD_opt')

    # Fix broken angles
    mol = fix_carboxyl(fix_h(mol))
    mol.job_geometry_opt(job2, s2, name='QD_opt')

    # Write the reuslts to an .xyz and .pdb file
    mol.properties.name += '.opt'
    export_mol(mol, message='Optimized core + ligands:\t\t')
    mol.properties.name = mol.properties.name.split('.opt')[0]

    return mol


def init_qd_opt(mol, database, arg):
    """
    Check if the to be optimized quantom dot has previously been optimized.
    Pull if the structure from the database if it has, otherwise perform a geometry optimization.
    mol <plams.Molecule> The input quantom dot with the 'name' property.
    database <pd.DataFrame>: A database of previous calculations.
    return <plams.Molecule>: An optimized quantom dot.
    """
    name = mol.properties.name.rsplit('.', 1)[0]
    if database is None:
        mol = qd_opt(mol, arg['maxiter'])
        mol.properties.entry = False
    elif database.empty or name not in list(database['Quantum_dot_name']):
        mol = qd_opt(mol, arg['maxiter'])
        mol.properties.entry = True
    else:
        index = list(database['Quantum_dot_name']).index(name)
        try:
            mol_new = molkit.readpdb(database['Quantum_dot_opt_pdb'][index])
            mol_new.properties = mol.properties
            mol = mol_new
        except FileNotFoundError:
            mol = qd_opt(mol, arg['maxiter'])
            mol.properties.entry = True

    return mol
