""" A module designed for optimizing the combined ligand & core. """

__all__ = ['init_qd_opt']

from scm.plams.core.settings import Settings
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob

from ..utils import get_time
from ..mol_utils import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt
from ..data_handling.CAT_database import (Database, mol_to_file)


def init_qd_opt(qd_df, arg):
    """ Initialized the quantum dot (constrained) geometry optimization.
    performs an inplace update of the *mol* column in **qd_df**.

    :parameter qd_df: A dataframe of quantum dots.
    :type qd_df: |pd.DataFrame|_ (columns: |str|_, index: |str|_, values: |plams.Molecule|_)
    :parameter arg: A settings object containing all (optional) arguments.
    :type arg: |plams.Settings|_ (superclass: |dict|_).
    """
    # Prepare slices
    job_recipe = arg.optional.qd.optimize
    overwrite = 'qd' in arg.optional.database.overwrite
    if overwrite:
        idx = qd_df['hdf5 index']
        message = '\t has been (re-)optimized'
    else:
        idx = qd_df['hdf5 index'] < 0
        message = '\t has been optimized'

    # Optimize the geometries
    if idx.any():
        init(path=arg.optional.qd.dirname, folder='QD_optimize')
        for mol in qd_df['mol'][idx]:
            qd_opt(mol, job_recipe)
            print(get_time() + mol.properties.name + message)
        finish()
        print('')

    # Export the geometries to the database
    if 'qd' in arg.optional.database.write:
        recipe = Settings()
        recipe.settings1 = {
            'name': '1',
            'key': job_recipe.job1,
            'value': job_recipe.s1,
            'template': 'geometry.json'
        }
        recipe.settings2 = {
            'name': '2',
            'key': job_recipe.job2,
            'value': job_recipe.s2,
            'template': 'geometry.json'
        }
        path = arg.optional.database.dirname
        database = Database(path=path)
        database.update_csv(qd_df, columns=['hdf5 index'], job_recipe=recipe, database='QD')
        mol_to_file(qd_df['mol'], path, overwrite, arg.optional.database.mol_format)


def qd_opt(mol, job_recipe):
    """ """
    if job_recipe.job1 is AMSJob:
        job_recipe.s1.input.ams.constraints.atom = mol.properties.indices
    if job_recipe.job2 is AMSJob:
        job_recipe.s2.input.ams.constraints.atom = mol.properties.indices

    # Prepare the job settings
    mol.job_geometry_opt(job_recipe.job1, job_recipe.s1, name='QD_opt_part1')

    # Fix broken angles
    fix_carboxyl(mol)
    fix_h(mol)
    mol.job_geometry_opt(job_recipe.job2, job_recipe.s2, name='QD_opt_part2')
