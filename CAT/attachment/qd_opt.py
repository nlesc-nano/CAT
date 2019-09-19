"""
CAT.attachment.qd_opt
=====================

A module designed for optimizing the combined ligand & core.

Index
-----
.. currentmodule:: CAT.attachment.qd_opt
.. autosummary::
    init_qd_opt
    start_qd_opt
    get_job_settings
    _qd_to_db
    qd_opt

API
---
.. autofunction:: init_qd_opt
.. autofunction:: start_qd_opt
.. autofunction:: get_job_settings
.. autofunction:: _qd_to_db
.. autofunction:: qd_opt

"""

from typing import List, Tuple

import pandas as pd

from scm.plams import (Molecule, Settings)
from scm.plams.core.functions import finish
from scm.plams.interfaces.adfsuite.ams import AMSJob

import qmflows

from .qd_opt_ff import qd_opt_ff
from ..logger import logger
from ..jobs import job_geometry_opt
from ..utils import (restart_init, type_to_string)
from ..mol_utils import (fix_carboxyl, fix_h, round_coords)
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.mol_to_file import mol_to_file

__all__ = ['init_qd_opt']

# Aliases for pd.MultiIndex columns
MOL: Tuple[str, str] = ('mol', '')
OPT: Tuple[str, str] = ('opt', '')
HDF5_INDEX: Tuple[str, str] = ('hdf5 index', '')
JOB_SETTINGS_QD_OPT: Tuple[str, str] = ('job_settings_QD_opt', '')
SETTINGS1: Tuple[str, str] = ('settings', '1')
SETTINGS2: Tuple[str, str] = ('settings', '2')


def init_qd_opt(qd_df: SettingsDataFrame) -> None:
    """Initialize the quantum dot (constrained) geometry optimization.

    performs an inplace update of the *mol* column in **qd_df**.

    Parameters
    ----------
    qd_df : |CAT.SettingsDataFrame|_
        A dataframe of quantum dots.

    """
    # Extract arguments
    settings = qd_df.settings.optional
    db = settings.database.db
    write = db and 'qd' in settings.database.write
    overwrite = db and 'qd' in settings.database.overwrite
    mol_format = settings.database.mol_format
    qd_path = settings.qd.dirname

    # Prepare slices
    if overwrite and db:
        idx = pd.Series(True, index=qd_df.index, name='mol')
    else:
        idx = qd_df[OPT] == False  # noqa

    # Optimize the geometries
    if idx.any():
        logger.info('Starting quantum dot optimization')
        start_qd_opt(qd_df, idx)
        qd_df[JOB_SETTINGS_QD_OPT] = get_job_settings(qd_df)
        logger.info('Finishing quantum dot optimization\n')
    else:  # No new molecules, move along
        logger.info('No new to-be optimized quantum dots found\n')
        return None

    # Export the geometries to the database
    if write and db:
        with pd.option_context('mode.chained_assignment', None):
            _qd_to_db(qd_df, idx)

    # Export xyz/pdb files
    if 'qd' in settings.database.write and mol_format:
        mol_to_file(qd_df[MOL], qd_path, mol_format=mol_format)

    return None


def start_qd_opt(qd_df: SettingsDataFrame, idx: pd.Series) -> None:
    """Loop over all molecules in ``qd_df.loc[idx]`` and perform geometry optimizations."""
    # Extract arguments
    path = qd_df.settings.optional.qd.dirname
    job_recipe = qd_df.settings.optional.qd.optimize
    forcefield = bool(qd_df.settings.optional.forcefield)

    # Perform the main optimization loop
    restart_init(path=path, folder='QD_optimize')
    for mol in qd_df[MOL][idx]:
        mol.properties.job_path = []
        qd_opt(mol, job_recipe, forcefield=forcefield)
    finish()


def get_job_settings(qd_df: SettingsDataFrame) -> List[str]:
    """Create a nested list of input files for each molecule in **ligand_df**."""
    job_settings = []
    for mol in qd_df[MOL]:
        try:
            job_settings.append(mol.properties.pop('job_path'))
        except KeyError:
            job_settings.append([])
    return job_settings


def _qd_to_db(qd_df: SettingsDataFrame, idx: pd.Series) -> None:
    """Export quantum dot optimziation results to the database.

    Parameters
    ----------
    qd_df : |CAT.SettingsDataFrame|_
        A dataframe of quantum dots.

    idx : |pd.Series|_
        A Series for slicing **qd_df**.

    """
    # Extract arguments
    settings = qd_df.settings.optional
    job_recipe = settings.qd.optimize
    overwrite = 'qd' in settings.database.overwrite
    db = settings.database.db

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
    db.update_csv(
        qd_df[idx],
        overwrite=overwrite,
        columns=columns,
        job_recipe=recipe,
        database='QD',
        opt=True
    )


def qd_opt(mol: Molecule, job_recipe: Settings, forcefield: bool = False) -> None:
    """Perform an optimization of the quantum dot.

    Performs an inplace update of **mol**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        The to-be optimized molecule.

    job_recipe : |plams.Settings|_
        A Settings instance containing all jon settings.
        Expects 4 keys: ``"job1"``, ``"job2"``, ``"s1"``, ``"s2"``.

    forcefield : bool
        If ``True``, perform the job with CP2K with a user-specified forcefield.

    """
    if job_recipe.job1 is AMSJob:
        job_recipe.s1.input.ams.constraints.atom = mol.properties.indices
    if job_recipe.job2 is AMSJob:
        job_recipe.s2.input.ams.constraints.atom = mol.properties.indices

    # Prepare the job settings
    if forcefield:
        qd_opt_ff(mol, job_recipe)
        return None

    job1, s1 = job_recipe.job1, job_recipe.s1
    job2, s2 = job_recipe.job2, job_recipe.s2

    # Run the first job and fix broken angles
    mol.job_geometry_opt(job1, s1, name='QD_opt_part1')
    fix_carboxyl(mol)
    fix_h(mol)
    mol.round_coords()

    # Run the second job
    mol.job_geometry_opt(job2, s2, name='QD_opt_part2')
    mol.round_coords()
    return None
