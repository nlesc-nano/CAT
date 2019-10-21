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

from typing import List, Tuple, Iterable, Optional, Type

from scm.plams import Molecule, Settings, AMSJob
from scm.plams.core.basejob import Job

from .qd_opt_ff import qd_opt_ff
from ..jobs import job_geometry_opt
from ..mol_utils import (fix_carboxyl, fix_h, round_coords)
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.mol_to_file import mol_to_file
from ..workflows.workflow import WorkFlow

__all__ = ['init_qd_opt']

# Aliases for pd.MultiIndex columns
MOL: Tuple[str, str] = ('mol', '')
OPT: Tuple[str, str] = ('opt', '')
HDF5_INDEX: Tuple[str, str] = ('hdf5 index', '')
JOB_SETTINGS_QD_OPT: Tuple[str, str] = ('job_settings_QD_opt', '')
SETTINGS1: Tuple[str, str] = ('settings', '1')
SETTINGS2: Tuple[str, str] = ('settings', '2')


def init_qd_opt(qd_df: SettingsDataFrame) -> None:
    """Initialize the ligand optimization procedure."""
    workflow = WorkFlow.from_template(qd_df, name='qd_opt')

    # Pull from the database; push unoptimized structures
    idx = workflow.from_db(qd_df)
    workflow(start_qd_opt, qd_df, index=idx)
    qd_df.loc[idx, OPT] = True
    qd_df.loc[idx, JOB_SETTINGS_QD_OPT] = pop_job_settings(qd_df)

    # Push the optimized structures to the database
    job_recipe = workflow.get_recipe()
    workflow.to_db(qd_df, status='optimized', idx_slice=idx, job_recipe=job_recipe)

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    mol_format = qd_df.settings.database.mol_format
    if mol_format:
        path = workflow.dirname
        mol_to_file(qd_df.loc[idx, MOL], path, mol_format=mol_format)


def start_qd_opt(mol_list: Iterable[Molecule],
                 jobs: Tuple[Optional[Type[Job]], ...], settings: Tuple[Optional[Settings], ...],
                 forcefield=None, **kwargs) -> None:
    """Loop over all molecules in ``qd_df.loc[idx]`` and perform geometry optimizations."""
    for mol in mol_list:
        mol.properties.job_path = []
        qd_opt(mol, jobs, settings, forcefield=forcefield)


def pop_job_settings(qd_df: SettingsDataFrame) -> List[str]:
    """Create a nested list of input files for each molecule in **ligand_df**."""
    job_settings = []
    for mol in qd_df[MOL]:
        try:
            job_settings.append(mol.properties.pop('job_path'))
        except KeyError:
            job_settings.append([])
    return job_settings


def qd_opt(mol: Molecule, jobs: Tuple[Optional[Type[Job]], ...],
           settings: Tuple[Optional[Settings], ...], forcefield: bool = False) -> None:
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
    # Prepare the job settings
    if forcefield:
        qd_opt_ff(mol, jobs, settings)
        return None

    # Expand arguments
    job1, job2 = jobs
    s1, s2 = settings

    # Extra options for AMSJob
    if job1 is AMSJob:
        s1 = s1.copy()
        s1.input.ams.constraints.atom = mol.properties.indices
    if job2 is AMSJob:
        s2 = s2.copy()
        s2.input.ams.constraints.atom = mol.properties.indices

    # Run the first job and fix broken angles
    mol.job_geometry_opt(job1, s1, name='QD_opt_part1')
    fix_carboxyl(mol)
    fix_h(mol)
    mol.round_coords()

    # Run the second job
    mol.job_geometry_opt(job2, s2, name='QD_opt_part2')
    mol.round_coords()
    return None
