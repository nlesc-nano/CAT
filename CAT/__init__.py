"""
CAT
===

A collection of tools designed for the construction of various chemical compounds.

"""

from scm.plams import (
    Settings as _Settings,
    add_to_class
)

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
    read_mol, set_mol_prop
)

from .utils import get_template

if hasattr(_Settings, 'suppress_missing'):
    _Settings.supress_missing = _Settings.suppress_missing

    @add_to_class(_Settings)
    def find_case(self, key):
        """Check if this instance contains a key consisting of the same letters as *key*, but possibly with different case.

        If found, return such a key. If not, return *key*.
        """  # noqa: E501
        if not isinstance(key, str):
            return key
        lowkey = key.lower()
        for k in self:
            try:
                if k.lower() == lowkey:
                    return k
            except (AttributeError, TypeError):
                pass
        return key

    del add_to_class


__version__ = __version__
__author__ = "Bas van Beek"
__email__ = 'b.f.van.beek@vu.nl'

__all__ = [
    'job_single_point', 'job_geometry_opt', 'job_freq',

    'get_thermo', 'get_entropy',

    'init_qd_opt', 'init_ligand_opt', 'init_qd_construction', 'init_ligand_anchoring',

    'read_mol', 'set_mol_prop',

    'get_template',
]
