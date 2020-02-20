"""
CAT.recipes.mark_surface
========================

A recipe for identifying surface-atom subsets.

Index
-----
.. currentmodule:: CAT.recipes.mark_surface
.. autosummary::
    replace_surface

API
---
.. autofunction:: replace_surface

"""

from typing import Any, Union

import numpy as np

from scm.plams import Molecule

from CAT.mol_utils import to_atnum
from CAT.attachment.distribution import distribute_idx
try:
    from nanoCAT.bde.identify_surface import identify_surface
except ImportError as ex:
    tb = ex.__traceback__
    raise ImportError("Executing the content of '{__file__}' requires the Nano-CAT package: "
                      "'https://github.com/nlesc-nano/nano-CAT'").with_traceback(tb)

__all__ = ['replace_surface']


def replace_surface(mol: Molecule,
                    symbol: Union[str, int],
                    symbol_new: Union[str, int] = 'Cl',
                    f: float = 0.5,
                    mode: str = 'uniform',
                    **kwargs: Any) -> Molecule:
    r"""A workflow for identifying all surface atoms in **mol** and replacing a subset of them.

    Consists of three distinct steps:

    1. Identifying which atoms, with a user-specified atomic **symbol**,
       are located on the surface of **mol** rather than in the bulk.
    2. Define a subset of the newly identified surface atoms using one of CAT's
       distribution algorithms.
    3. Create and return a molecule where the atom subset defined in step 2
       has its atomic symbols replaced with **symbol_new**.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule
        >>> from CAT.recipes import replace_surface

        >>> mol = Molecule(...)  # Read an .xyz file
        >>> mol_new = replace_surface(mol, symbol='Cl', symbol_new='I', f=0.75)
        >>> mol_new.write(...)  # Write an .xyz file


    Parameters
    ----------
    mol : :class:`Molecule<scm.plams.mol.molecule.Molecule>`
        The input molecule.

    symbol : :class:`str` or :class:`int`
        An atomic symbol or number defining the super-set of the surface atoms.

    symbol_new : :class:`str` or :class:`int`
        An atomic symbol or number which will be assigned to the new surface-atom subset.

    f : :class:`float`
        The fraction of surface atoms whose atom types will be replaced with **symbol_new**.
        Must obey the following condition: :math:`0 < f \le 1`.

    mode : :class:`str`
        How the subset of surface atoms will be generated.
        Accepts one of the following values:

        * ``"random"``: A random distribution.
        * ``"uniform"``: A uniform distribution; maximizes the nearest-neighbor distance.
        * ``"cluster"``: A clustered distribution; minimizes the nearest-neighbor distance.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for
        :func:`distribute_idx()<CAT.attachment.distribution.distribute_idx>`.

    Returns
    -------
    :class:`Molecule<scm.plams.mol.molecule.Molecule>`
        A new Molecule with a subset of its surface atoms replaced with **symbol_new**.

    See Also
    --------
    :func:`distribute_idx()<CAT.attachment.distribution.distribute_idx>`
        Create a new distribution of atomic indices from **idx** of length :code:`f * len(idx)`.

    """
    # Parse input arguments
    xyz = np.asarray(mol)
    atnum = to_atnum(symbol)
    atnum_new = to_atnum(symbol_new)

    # Define the surface-atom subset
    idx = np.fromiter((i for i, at in enumerate(mol) if at.atnum == atnum), dtype=int)
    idx_surface = idx[identify_surface(xyz[idx])]
    idx_surface_subset = distribute_idx(xyz, idx_surface, f=f, mode=mode, **kwargs)
    idx_surface_subset += 1

    # Update atomic symbols and return
    ret = mol.copy()
    for i in idx_surface_subset:
        ret[i].atnum = atnum_new
    return ret
