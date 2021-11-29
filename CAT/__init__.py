"""**CAT**: A collection of tools designed for the construction of various chemical compounds."""  # noqa: E501

# flake8: noqa: E402

from .__version__ import __version__ as __version__
__author__ = "Bas van Beek"
__email__ = 'b.f.van.beek@vu.nl'

from nanoutils import VersionInfo
version_info = VersionInfo.from_str(__version__)
del VersionInfo

from .logger import logger
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
    read_mol, set_mol_prop
)
from .utils import get_template

from . import distribution

__all__ = [
    'logger',

    'job_single_point', 'job_geometry_opt', 'job_freq',

    'get_thermo', 'get_entropy',

    'init_qd_opt', 'init_ligand_opt', 'init_qd_construction', 'init_ligand_anchoring',

    'read_mol', 'set_mol_prop',

    'get_template',
]
