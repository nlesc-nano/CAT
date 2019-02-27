""" A module designed for optimizing the combined ligand & core. """

__all__ = ['init_qd_opt']

from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob

from ..mol_utils import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt
from ..data_handling.database import qd_to_database
from ..data_handling.mol_export import export_mol


def init_qd_opt(mol_list, arg):
    """
    Check if the to be optimized quantom dot has previously been optimized.
    Pull if the structure from the database if it has, otherwise perform a geometry optimization.
    mol_list <list> [<plams.Molecule>]: The list of input quantom dots with the 'name' property.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    # Optimize all geometries
    job_recipe = arg.optional.qd.optimize
    for mol in mol_list:
        qd_opt(mol, job_recipe)

    # Export the geometries to the database
    if 'qd' in arg.optional.database.write:
        qd_to_database(mol_list, arg)


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
