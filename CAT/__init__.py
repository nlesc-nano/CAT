"""
CAT
===

A collection of tools designed for the construction of various chemical compounds.

"""

from .__version__ import __version__

from .jobs import (
    job_single_point, job_geometry_opt, job_freq
)

from .thermo_chem import (
    get_thermo, get_entropy
)

from .attachment import (
    init_qd_opt, init_ligand_opt, init_qd_construction, init_ligand_anchoring
)

from .data_handling import (
    read_mol, set_mol_prop,
    sanitize_optional, sanitize_input_mol, sanitize_path
)

from .base import prep

from .utils import get_template

__version__ = __version__
__author__ = "Bas van Beek"
__email__ = 'b.f.van.beek@vu.nl'

__all__ = [
    'job_single_point', 'job_geometry_opt', 'job_freq',

    'get_thermo', 'get_entropy',

    'init_qd_opt', 'init_ligand_opt', 'init_qd_construction', 'init_ligand_anchoring',

    'read_mol', 'set_mol_prop',
    'sanitize_optional', 'sanitize_input_mol', 'sanitize_path',

    'prep',

    'get_template'
]
