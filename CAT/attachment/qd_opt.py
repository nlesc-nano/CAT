""" A module designed for optimizing the combined ligand & core. """

__all__ = ['init_qd_opt']

from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob

from ..utils import get_time
from ..mol_utils import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt
from ..data_handling.database import qd_to_database


def init_qd_opt(mol_list, arg):
    """
    Check if the to be optimized quantom dot has previously been optimized.
    Pull if the structure from the database if it has, otherwise perform a geometry optimization.

    mol_list <list> [<plams.Molecule>]: The list of input quantom dots with the 'name' property.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    # Optimize all geometries
    job_recipe = arg.optional.qd.optimize
    overwrite = 'qd' in arg.optional.database.overwrite
    init(path=arg.optional.qd.dirname, folder='QD_optimize')
    for mol in mol_list:
        if overwrite or not mol.properties.read:
            qd_opt(mol, job_recipe)
            if mol.properties.read:
                print(get_time() + mol.properties.name + '\t has been reoptimized')
            else:
                print(get_time() + mol.properties.name + '\t has been optimized')
    finish()
    print('')

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
    mol.job_geometry_opt(job_recipe.job1, job_recipe.s1, name='QD_opt_part1')

    # Fix broken angles
    mol = fix_carboxyl(fix_h(mol))
    mol.job_geometry_opt(job_recipe.job2, job_recipe.s2, name='QD_opt_part2')

    return mol
