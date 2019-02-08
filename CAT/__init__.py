""" A collection of tools designed for the construction,
and subsequent analysis, of various chemical compounds. """

__author__ = "Bas van Beek"
__email__ = 'b.f.van.beek@vu.nl'

from .__version__ import __version__

from .qd_functions import (find_substructure, find_substructure_split,
                           merge_mol, adf_connectivity, fix_h, fix_carboxyl)

__all__ = [
    'find_substructure', 'find_substructure_split', 'merge_mol',
    'adf_connectivity', 'fix_h', 'fix_carboxyl',
]
