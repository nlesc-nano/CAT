"""A module designed for optimizing the combined ligand & core."""

import pandas as pd

from scm.plams import (Molecule, Settings)
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob

import qmflows
from data_CAT import (Database, mol_to_file)

from ..properties_dataframe import PropertiesDataFrame
from ..utils import (get_time, type_to_string)
from ..mol_utils import (fix_carboxyl, fix_h)
from ..analysis.jobs import job_geometry_opt

__all__ = ['init_qd_opt']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
OPT = ('opt', '')
HDF5_INDEX = ('hdf5 index', '')
JOB_SETTINGS_QD_OPT = ('job_settings_QD_opt', '')
SETTINGS1 = ('settings', '1')
SETTINGS2 = ('settings', '2')


def init_qd_opt(qd_df: PropertiesDataFrame) -> None:
    """Initialized the quantum dot (constrained) geometry optimization.

    performs an inplace update of the *mol* column in **qd_df**.

    Parameters
    ----------
    qd_df : |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    """
    # Extract arguments
    properties = qd_df.properties
    write = 'qd' in properties.optional.database.write
    job_recipe = properties.optional.qd.optimize
    overwrite = 'qd' in properties.optional.database.overwrite
    path = properties.optional.qd.dirname

    # Prepare slices
    if overwrite:
        idx = pd.Series(True, index=qd_df.index, name='mol')
        message = '\t has been (re-)optimized'
    else:
        idx = qd_df[OPT] == False  # noqa
        message = '\t has been optimized'

    # Optimize the geometries
    if idx.any():
        init(path=path, folder='QD_optimize')
        for mol in qd_df[MOL][idx]:
            mol.properties.job_path = []
            qd_opt(mol, job_recipe)
            print(get_time() + mol.properties.name + message)
        finish()

        job_settings = []
        for mol in qd_df[MOL]:
            try:
                job_settings.append(mol.properties.pop('job_path'))
            except KeyError:
                job_settings.append([])
        qd_df[JOB_SETTINGS_QD_OPT] = job_settings
        print()
    else:  # No new molecules, move along
        return None

    # Export the geometries to the database
    if write:
        with pd.option_context('mode.chained_assignment', None):
            _qd_to_db(qd_df, idx)
    return None


def _qd_to_db(qd_df: PropertiesDataFrame,
              idx: pd.Series) -> None:
    """Export quantum dot optimziation results to the database.

    Parameters
    ----------
    qd_df : |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    idx : |pd.Series|_
        A Series for slicing **qd_df**.

    """
    # Extract arguments
    properties = qd_df.properties
    job_recipe = properties.optional.qd.optimize
    overwrite = 'qd' in properties.optional.database.overwrite
    mol_format = properties.optional.database.mol_format
    path = properties.optional.qd.dirname

    # Preapre the job recipe
    v1 = qmflows.geometry['specific'][type_to_string(job_recipe.job1)].copy()
    v1.update(job_recipe.s1)
    v2 = qmflows.geometry['specific'][type_to_string(job_recipe.job2)].copy()
    v2.update(job_recipe.s2)
    recipe = Settings({
        '1': {'key': job_recipe.job1, 'value': v1},
        '2': {'key': job_recipe.job2, 'value': v2}
    })

    # Update the database
    columns = [HDF5_INDEX, JOB_SETTINGS_QD_OPT, SETTINGS1, SETTINGS2]
    database = Database(path=path)
    database.update_csv(
        qd_df[idx],
        columns=columns,
        job_recipe=recipe,
        database='QD',
        opt=True
    )

    # Export xyz/pdb files
    mol_to_file(qd_df[MOL], path, overwrite, mol_format)


def qd_opt(mol: Molecule,
           job_recipe: Settings) -> None:
    """Perform an optimization of the quantum dot.

    Performs an inplace update of **mol**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        The to-be optimized molecule.

    job_recipe : |plams.Settings|_
        A Settings instance containing all jon settings.
        Expects 4 keys: ``"job1"``, ``"job2"``, ``"s1"``, ``"s2"``.

    """
    if job_recipe.job1 is AMSJob:
        job_recipe.s1.input.ams.constraints.atom = mol.properties.indices
    if job_recipe.job2 is AMSJob:
        job_recipe.s2.input.ams.constraints.atom = mol.properties.indices

    # Prepare the job settings
    job1, s1 = job_recipe.job1, job_recipe.s1
    mol.job_geometry_opt(job1, s1, name='QD_opt_part1')

    # Fix broken angles
    fix_carboxyl(mol)
    fix_h(mol)
    job2, s2 = job_recipe.job2, job_recipe.s2
    mol.job_geometry_opt(job2, s2, name='QD_opt_part2')
