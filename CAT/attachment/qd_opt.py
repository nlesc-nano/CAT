"""A module designed for optimizing the combined ligand & core."""

from typing import List

import pandas as pd

from scm.plams import (Molecule, Settings)
from scm.plams.core.functions import (init, finish)
from scm.plams.interfaces.adfsuite.ams import AMSJob

import qmflows

from ..jobs import job_geometry_opt
from ..utils import (get_time, type_to_string)
from ..mol_utils import (fix_carboxyl, fix_h)
from ..properties_dataframe import PropertiesDataFrame

try:
    from dataCAT import (Database, mol_to_file)
    DATA_CAT = True
except ImportError:
    DATA_CAT = False

__all__ = ['init_qd_opt']

# Aliases for pd.MultiIndex columns
MOL = ('mol', '')
OPT = ('opt', '')
HDF5_INDEX = ('hdf5 index', '')
JOB_SETTINGS_QD_OPT = ('job_settings_QD_opt', '')
SETTINGS1 = ('settings', '1')
SETTINGS2 = ('settings', '2')


def init_qd_opt(qd_df: PropertiesDataFrame) -> None:
    """Initialize the quantum dot (constrained) geometry optimization.

    performs an inplace update of the *mol* column in **qd_df**.

    Parameters
    ----------
    qd_df : |CAT.PropertiesDataFrame|_
        A dataframe of quantum dots.

    """
    # Extract arguments
    properties = qd_df.properties
    write = DATA_CAT and 'qd' in properties.optional.database.write
    overwrite = DATA_CAT and 'qd' in properties.optional.database.overwrite

    # Prepare slices
    if overwrite and DATA_CAT:
        idx = pd.Series(True, index=qd_df.index, name='mol')
        message = '\t has been (re-)optimized'
    else:
        idx = qd_df[OPT] == False  # noqa
        message = '\t has been optimized'

    # Optimize the geometries
    if idx.any():
        start_qd_opt(qd_df, idx, message)
        qd_df[JOB_SETTINGS_QD_OPT] = get_job_settings(qd_df)
        print()
    else:  # No new molecules, move along
        return None

    # Export the geometries to the database
    if write and DATA_CAT:
        with pd.option_context('mode.chained_assignment', None):
            _qd_to_db(qd_df, idx)
    return None


def start_qd_opt(qd_df: PropertiesDataFrame,
                 idx: pd.Series,
                 message: str) -> None:
    """Loop over all molecules in ``qd_df.loc[idx]`` and perform geometry optimizations."""
    # Extract arguments
    path = qd_df.properties.optional.qd.dirname
    job_recipe = qd_df.properties.optional.qd.optimize

    # Perform the main optimization loop
    init(path=path, folder='QD_optimize')
    for mol in qd_df[MOL][idx]:
        mol.properties.job_path = []
        qd_opt(mol, job_recipe)
        print(get_time() + mol.properties.name + message)
    finish()


def get_job_settings(qd_df: PropertiesDataFrame) -> List[str]:
    """Create a nested list of input files for each molecule in **ligand_df**."""
    job_settings = []
    for mol in qd_df[MOL]:
        try:
            job_settings.append(mol.properties.pop('job_path'))
        except KeyError:
            job_settings.append([])
    return job_settings


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
    overwrite = DATA_CAT and 'qd' in properties.optional.database.overwrite
    mol_format = properties.optional.database.mol_format
    qd_path = properties.optional.qd.dirname
    db_path = properties.optional.database.dirname

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
    database = Database(path=db_path)
    database.update_csv(
        qd_df[idx],
        columns=columns,
        job_recipe=recipe,
        database='QD',
        opt=True
    )

    # Export xyz/pdb files
    mol_to_file(qd_df[MOL], qd_path, overwrite, mol_format)


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
