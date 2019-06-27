""" A module designed for optimizing the combined ligand & core. """

__all__ = ['init_qd_opt']

import pandas as pd

from scm.plams.core.settings import Settings
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob

import qmflows

from ..utils import (get_time, type_to_string)
from ..mol_utils import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt
from ..data_handling.database import Database
from ..data_handling.database_functions import mol_to_file


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
        idx = pd.Series(True, index=qd_df.index, name='mol')
        message = '\t has been (re-)optimized'
    else:
        idx = qd_df['opt'] == False  # noqa
        message = '\t has been optimized'

    # Optimize the geometries
    if idx.any():
        init(path=arg.optional.qd.dirname, folder='QD_optimize')
        for mol in qd_df['mol'][idx]:
            qd_opt(mol, job_recipe)
            print(get_time() + mol.properties.name + message)
        finish()

    qd_df['job_settings_QD_opt'] = [mol.properties.pop('job_path') for mol in qd_df['mol']]
    for mol in qd_df['mol']:
        mol.properties.job_path = []

    # Export the geometries to the database
    if 'qd' in arg.optional.database.write:
        _qd_to_db(qd_df, arg)


def _qd_to_db(qd_df, arg):
    """Export quantum dot optimziation results to the database."""
    job_recipe = arg.optional.qd.optimize
    overwrite = 'qd' in arg.optional.database.overwrite

    v1 = qmflows.geometry['specific'][type_to_string(job_recipe.job1)].copy()
    v1.update(job_recipe.s1)
    v2 = qmflows.geometry['specific'][type_to_string(job_recipe.job2)].copy()
    v2.update(job_recipe.s2)
    recipe = Settings({
        '1': {'key': job_recipe.job1, 'value': v1},
        '2': {'key': job_recipe.job2, 'value': v2}
    })

    columns = [('hdf5 index', ''), ('settings', '1'), ('settings', '2')]
    database = Database(path=arg.optional.database.dirname)
    database.update_csv(qd_df, columns=columns, job_recipe=recipe, database='QD', opt=True)
    path = arg.optional.qd.dirname

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
