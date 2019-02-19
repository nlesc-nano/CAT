""" A module designed for optimizing the combined ligand & core. """

__all__ = ['init_qd_opt']

from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob
import scm.plams.interfaces.molecule.rdkit as molkit

from ..mol_utils import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt
from ..data_handling.mol_export import export_mol


def qd_opt(mol, job_recipe):
    """ """
    if job_recipe.job1 is AMSJob:
        job_recipe.s1.input.ams.constraints.atom = mol.properties.indices
    if job_recipe.job2 is AMSJob:
        job_recipe.s2.input.ams.constraints.atom = mol.properties.indices

    # Prepare the job settings
    init(path=mol.properties.path, folder='QD_opt')
    mol.job_geometry_opt(job_recipe.job1, job_recipe.s1, name='QD_opt_part1')

    # Fix broken angles
    mol = fix_carboxyl(fix_h(mol))
    mol.job_geometry_opt(job_recipe.job2, job_recipe.s2, name='QD_opt_part2')

    # Write the reuslts to an .xyz and .pdb file
    mol.properties.name += '.opt'
    export_mol(mol, message='Optimized core + ligands:\t\t')
    mol.properties.name = mol.properties.name.split('.opt')[0]
    finish()

    return mol


def init_qd_opt(mol, database, job_recipe):
    """
    Check if the to be optimized quantom dot has previously been optimized.
    Pull if the structure from the database if it has, otherwise perform a geometry optimization.
    mol <plams.Molecule> The input quantom dot with the 'name' property.
    database <pd.DataFrame>: A database of previous calculations.
    return <plams.Molecule>: An optimized quantom dot.
    """
    name = mol.properties.name.rsplit('.', 1)[0]
    if database is None:
        mol = qd_opt(mol, job_recipe)
        mol.properties.entry = False
    elif database.empty or name not in list(database['Quantum_dot_name']):
        mol = qd_opt(mol, job_recipe)
        mol.properties.entry = True
    else:
        index = list(database['Quantum_dot_name']).index(name)
        try:
            mol_new = molkit.readpdb(database['Quantum_dot_opt_pdb'][index])
            mol_new.properties = mol.properties
            mol = mol_new
        except FileNotFoundError:
            mol = qd_opt(mol, job_recipe)
            mol.properties.entry = True

    return mol
