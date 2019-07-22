"""
CAT.data_handling.validate_input
================================

A module designed for sanitizing and interpreting the input file.

Index
-----
.. currentmodule:: CAT.data_handling.validate_input
.. autosummary::
    validate_input

API
---
.. autofunction:: validate_input

"""

from os import mkdir
from os.path import (join, isdir)

from scm.plams import Settings

from CAT.data_handling.validation_schemas import (
    core_schema, ligand_schema, qd_schema, database_schema,
    mongodb_schema, bde_schema, qd_opt_schema, crs_schema
)
from .validate_mol import validate_mol
from ..utils import validate_path

__all__ = ['validate_input']


def validate_input(s: Settings) -> None:
    """Initialize the input-validation procedure.

    performs an inplace update of **s**.

    Parameters
    ----------
    s : |plams.Settings|_
        A Settings instance with to-be validated CAT input settings.

    """
    # Validate the path
    s.path = path = validate_path(s.path)

    try:
        # Set the various working directories
        dirnames = ('database', 'ligand', 'core', 'qd')
        for key in dirnames:
            value = join(path, key)
            s.optional[key].dirname = value
            if not isdir(value):
                mkdir(value)

        # Validate optional argument
        s.optional.database = database_schema.validate(s.optional.database)
        s.optional.ligand = ligand_schema.validate(s.optional.ligand)
        s.optional.core = core_schema.validate(s.optional.core)
        s.optional.qd = qd_schema.validate(s.optional.qd)

        # Validate some of the more complex optionala rguments
        if s.optional.database.mongodb:
            s.optional.database.mongodb = mongodb_schema.validate(s.optional.database.mongodb)
        if s.optional.qd.optimize:
            s.optional.qd.optimize = qd_opt_schema.validate(s.optional.qd.optimize)
        if s.optional.qd.dissociate:
            s.optional.qd.dissociate = bde_schema.validate(s.optional.qd.dissociate)
        if s.optional.ligand['cosmo-rs']:
            crs = s.optional.ligand.pop('cosmo-rs')
            s.optional.ligand.crs = crs_schema.validate(crs)

        # Validate the input cores and ligands
        validate_mol(s.input_cores, 'input_cores', join(path, 'core'))
        validate_mol(s.input_ligands, 'input_ligands', join(path, 'ligand'))
    except Exception as ex:
        raise ex.__class__(ex)
