""" A collection of tools designed for the construction,
and subsequent analysis, of various chemical compounds. """

__author__ = "Bas van Beek"
__email__ = 'b.f.van.beek@vu.nl'

from .__version__ import __version__

from .analysis import (
    init_asa, CRSJob, CRSResults, job_single_point, job_geometry_opt, job_freq,
    init_bde, get_thermo, get_entropy, dissociate_ligand, init_solv
)

from .attachment import (
    qd_opt, optimize_ligand, ligand_to_qd, bob_ligand, bob_core, substitution, multi_substitution
)

from .data_handling import (
    read_database, compare_database, write_database, export_mol, read_mol,
    set_prop, sanitize_arg_dict
)

from .base import prep

from .misc import get_template

__all__ = [
    'init_asa', 'CRSJob', 'CRSResults', 'job_single_point', 'job_geometry_opt', 'job_freq',
    'init_bde', 'get_thermo', 'get_entropy', 'dissociate_ligand', 'init_solv',

    'qd_opt', 'optimize_ligand', 'ligand_to_qd', 'substitution', 'multi_substitution',
    'bob_ligand', 'bob_core',

    'read_database', 'compare_database', 'write_database', 'export_mol', 'read_mol',
    'set_prop', 'sanitize_arg_dict',

    'prep',

    'get_template',
]
