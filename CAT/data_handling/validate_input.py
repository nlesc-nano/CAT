"""A module designed for sanitizing and interpreting the input file.

Index
-----
.. currentmodule:: CAT.data_handling.validate_input
.. autosummary::
    validate_input

API
---
.. autofunction:: validate_input

"""

import sys
from os import mkdir
from os.path import (join, isdir)

from rdkit.Chem import Mol
from scm.plams import Settings

from .validation_schemas import (
    core_schema,
    ligand_schema,
    qd_schema,
    database_schema,
    mongodb_schema,
    bde_schema,
    qd_opt_schema,
    crs_schema,
    asa_schema,
    ligand_opt_schema,
    subset_schema,
    multi_ligand_schema,
    cdft_schema,
    bulkiness_schema,
)

from .validate_ff import validate_ff, update_ff_jobs
from .validate_mol import validate_mol
from ..utils import validate_path, SetEnviron
from ..logger import logger
from ..attachment.ligand_anchoring import get_functional_groups

try:
    from dataCAT import Database
    DATA_CAT = True
except ImportError:
    DATA_CAT = False

if sys.version_info >= (3, 7):
    from contextlib import nullcontext
else:
    from contextlib2 import nullcontext

__all__ = ['validate_input']


def _validate_multi_lig(s: Settings) -> None:
    """Check that one (and only one!) of ``'f'`` and ``'dummy'`` is specified."""
    f = s.optional.qd.multi_ligand.f

    anchor = s.optional.qd.multi_ligand.anchor
    if anchor is None:
        anchor = s.optional.qd.multi_ligand.dummy
        s.optional.qd.multi_ligand.anchor = anchor
    del s.optional.qd.multi_ligand.dummy

    if f is anchor is None:
        raise ValueError("'.multi_ligand.f' and '.multi_ligand.anchor' cannot be "
                         "both unspecified or set to 'None'")
    elif None not in (f, anchor):
        raise ValueError("Only one of '.multi_ligand.f' and '.multi_ligand.anchor' "
                         "should be specified")

    if anchor is not None:
        assert len(anchor) == len(s.optional.qd.multi_ligand.ligands)
    else:
        assert len(f) == len(s.optional.qd.multi_ligand.ligands) - 1


def validate_input(s: Settings, validate_only: bool = True) -> None:
    """Initialize the input-validation procedure.

    performs an inplace update of **s**.

    Parameters
    ----------
    s : |plams.Settings|_
        A Settings instance with to-be validated CAT input settings.
    validate_only : bool
        Perform only validation.

    """
    dirnames = ('database', 'ligand', 'core', 'qd')
    if not validate_only:
        # Validate the path
        s.path = validate_path(s.path)

        # Set the various working directories
        for key in dirnames:
            value = join(s.path, key)
            s.optional[key].dirname = value
            if not isdir(value):
                mkdir(value)
    else:
        if s.path is None:
            s.path = '.'
        for key in dirnames:
            s.optional[key].dirname = join(s.path, key)

    # Validate optional argument
    if not validate_only:
        context = nullcontext()
    else:
        context = SetEnviron(ADFBIN='a', ADFHOME='2019', ADFRESOURCES='b', SCMLICENSE='c')

    s.optional.database = database_schema.validate(s.optional.database)
    with context:
        s.optional.ligand = ligand_schema.validate(s.optional.ligand)
        s.optional.core = core_schema.validate(s.optional.core)
        s.optional.qd = qd_schema.validate(s.optional.qd)

        # Validate some of the more complex optionala rguments
        if s.optional.database.mongodb:
            s.optional.database.mongodb = mongodb_schema.validate(s.optional.database.mongodb)

        if s.optional.core.subset:
            s.optional.core.subset = subset_schema.validate(s.optional.core.subset)
            if 'p' in s.optional.core.subset:
                if 'weight' in s.optional.core.subset:
                    raise KeyError("'p' and 'weight' cannot be simultaneously specified")
                logger.warn("The 'subset.p' parameter is deprecated; see 'subset.weight'")
                p = s.optional.core.subset.pop('p')
                s.optional.core.subset.weight = lambda x: -(x**p)
        if s.optional.core.anchor is not None:
            pass
        elif s.optional.core.dummy is not None:
            s.optional.core.anchor = s.optional.core.dummy
        else:
            s.optional.core.anchor = 17
        del s.optional.core.dummy

        if s.optional.ligand.optimize:
            s.optional.ligand.optimize = ligand_opt_schema.validate(s.optional.ligand.optimize)
        if s.optional.ligand.cdft:
            s.optional.ligand.cdft = cdft_schema.validate(s.optional.ligand.cdft)

        if s.optional.get('crs') is not None:
            s.optional.ligand.crs = crs_schema.validate(s.optional.crs)
        elif s.optional.ligand['cosmo-rs']:
            crs = s.optional.ligand.pop('cosmo-rs')
            s.optional.ligand.crs = crs_schema.validate(crs)

        if s.optional.qd.optimize:
            s.optional.qd.optimize = qd_opt_schema.validate(s.optional.qd.optimize)
        if s.optional.qd.dissociate:
            s.optional.qd.dissociate = bde_schema.validate(s.optional.qd.dissociate)
        if s.optional.qd.activation_strain:
            s.optional.qd.activation_strain = asa_schema.validate(s.optional.qd.activation_strain)
        if s.optional.qd.multi_ligand:
            s.optional.qd.multi_ligand = multi_ligand_schema.validate(s.optional.qd.multi_ligand)
            _validate_multi_lig(s)
        if s.optional.qd.bulkiness:
            s.optional.qd.bulkiness = bulkiness_schema.validate(s.optional.qd.bulkiness)

        # Create forcefield Job Settings
        if s.optional.forcefield:
            s.optional.forcefield = validate_ff(s.optional.forcefield)
            update_ff_jobs(s)

        # Validate the input cores and ligands
        validate_mol(s.input_cores, 'input_cores', join(s.path, 'core'), not validate_only)
        validate_mol(s.input_ligands, 'input_ligands', join(s.path, 'ligand'), not validate_only)
        if 'input_qd' in s:
            validate_mol(s.input_qd, 'input_qd', join(s.path, 'qd'), not validate_only)

        # Create a dataCAT.Database instance
        if s.optional.database.get('db') is None and not validate_only:
            if DATA_CAT:
                db_path = s.optional.database.dirname
                s.optional.database.db = Database(path=db_path, **s.optional.database.mongodb)
            else:
                s.optional.database.db = False

        # Create RDKit molecules representing functional groups
        if s.optional.ligand.anchor is not None:
            func_groups = s.optional.ligand.anchor
        else:
            func_groups = s.optional.ligand.functional_groups
        del s.optional.ligand.functional_groups

        split = s.optional.ligand.split
        if func_groups is None:
            s.optional.ligand.anchor = get_functional_groups(None, split)
        elif not isinstance(func_groups[0], Mol):
            s.optional.ligand.anchor = get_functional_groups(func_groups)
