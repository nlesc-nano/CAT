"""
CAT.recipes.surface_dissociation
================================

A recipe for dissociation specific sets of surface atoms.

Index
-----
.. currentmodule:: CAT.recipes.surface_dissociation
.. autosummary::
    dissociate_surface
    row_accumulator

API
---
.. autofunction:: dissociate_surface
.. autofunction:: row_accumulator

"""

from typing import Iterable, Optional, Generator, Any

import numpy as np

from scm.plams import Molecule

from CAT.utils import get_nearest_neighbors
from CAT.mol_utils import to_atnum
from CAT.attachment.distribution_brute import brute_uniform_idx
try:
    from nanoCAT.bde.dissociate_xyn import dissociate_ligand
    from nanoCAT.bde.identify_surface import identify_surface
except ImportError as ex:
    tb = ex.__traceback__
    raise ImportError("Executing the content of '{__file__}' requires the Nano-CAT package: "
                      "'https://github.com/nlesc-nano/nano-CAT'").with_traceback(tb)

__all__ = ['dissociate_surface', 'row_accumulator']


def dissociate_surface(mol: Molecule,
                       idx: np.ndarray,
                       symbol: str = 'Cl',
                       lig_count: int = 1,
                       k: int = 4, **kwargs: Any) -> Generator[Molecule, None, None]:
    r"""A workflow for dissociating :math:`(XY_{n})_{\le m}` compounds from the surface of **mol**.

    The workflow consists of four distinct steps:

    1. Identify which atoms :math:`Y`, as specified by **symbol**,
       are located on the surface of **mol**.
    2. Identify which surface atoms are neighbors of :math:`X`, the latter being defined by **idx**.
    3. Identify which pairs of :math:`n*m` neighboring surface atoms are furthest removed from
       each other.
       :math:`n` is defined by **lig_count** and :math:`m`, if applicable, by the index along axis 1
       of **idx**.
    4. Yield :math:`(XY_{n})_{\le m}` molecules constructed from **mol**.

    Examples
    --------
    .. code:: python

        >>> from pathlib import path

        >>> import numpy as np

        >>> from scm.plams import Molecule
        >>> from CAT.recipes import dissociate_surface, row_accumulator

        >>> base_path = Path(...)
        >>> mol = Molecule(base_path / 'mol.xyz')

        # The indices of, e.g., Cs-pairs
        >>> idx = np.array([
        ...     [1, 3],
        ...     [4, 5],
        ...     [6, 10],
        ...     [15, 12],
        ...     [99, 105],
        ...     [20, 4]
        ... ])

        # Convert 1- to 0-based indices by substracting 1 from idx
        >>> mol_generator = dissociate_surface(mol, idx-1, symbol='Cl', lig_count=1)

        >>> iterator = zip(row_accumulator(idx), mol_generator)
        >>> for i, mol in iterator:
        ...     mol.write(base_path / f'output{i}.xyz')


    Parameters
    ----------
    mol : :class:`Molecule<scm.plams.mol.molecule.Molecule>`
        The input molecule.

    idx : array-like, dimensions: :math:`\le 2`
        An array of indices denoting to-be dissociated atoms (*i.e.* :math:`X`).
        If a 2D array is provided then all elements along axis 1 will be dissociated
        in a cumulative manner.
        :math:`m` is herein defined as the index along axis 1.

    symbol : :class:`str` or :class:`int`
        An atomic symbol or number defining the super-set of the atoms to-be dissociated in
        combination with **idx** (*i.e.* :math:`Y`).

    lig_count : :class:`int`
        The number of atoms specified in **symbol** to-be dissociated in combination
        with a single atom from **idx** (*i.e.* :math:`n`).

    k : :class:`int`
        The number of atoms specified in **symbol** which are surrounding a single atom in **idx**.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for
        :func:`brute_uniform_idx()<CAT.attachment.distribution_brute.brute_uniform_idx>`.

    Yields
    ------
    :class:`Molecule<scm.plams.mol.molecule.Molecule>`
        Yields new :math:`(XY_{n})_{m}`-dissociated molecules.

    See Also
    --------
    :func:`brute_uniform_idx()<CAT.attachment.distribution_brute.brute_uniform_idx>`
        Brute force approach to creating uniform or clustered distributions.

    :func:`identify_surface()<nanoCAT.bde.identify_surface.identify_surface>`
        Take a molecule and identify which atoms are located on the surface,
        rather than in the bulk.

    :func:`dissociate_ligand()<nanoCAT.bde.dissociate_xyn.dissociate_ligand>`
        Remove :math:`XY_{n}` from **mol** with the help of the
        :class:`MolDissociater<nanoCAT.bde.dissociate_xyn.MolDissociater>` class.

    """
    idx = np.array(idx, ndmin=2, copy=True)
    if idx.ndim > 2:
        raise ValueError("'idx' expected a 2D array-like object; "
                         f"observed dimensionality: {idx.ndim}D")
    idx.sort(axis=1)
    idx = idx[:, ::-1]

    # Identify all atoms in **idx** located on the surface
    idx_surface_superset = _get_surface(mol, symbol=symbol)

    # Construct an array with the indices of opposing surface-atoms
    n = lig_count * idx.shape[1]
    idx_surface = _get_opposite_neighbor(mol, idx, idx_surface_superset, n=n, k=k, **kwargs)

    # Dissociate and yield new molecules
    idx += 1
    idx_surface += 1
    for idx_pair, idx_pair_surface in zip(idx, idx_surface):
        mol_tmp = mol.copy()
        _mark_atoms(mol_tmp, idx_pair_surface)

        for i in idx_pair:
            mol_tmp = next(dissociate_ligand(mol_tmp, lig_count=lig_count,
                                             core_index=i, lig_core_pairs=1,
                                             **kwargs))
            yield mol_tmp


def row_accumulator(iterable: Iterable[Iterable[Any]]) -> Generator[str, None, None]:
    """Return a generator which accumulates elements along the nested elements of **iterable**.

    Examples
    --------
    .. code:: python

        >>> iterable = [[1, 3],
        ...             [4, 5],
        ...             [6, 10]]

        >>> for i in row_accumulator(iterable):
        ...     print(repr(i))
        '_1'
        '_1_3'
        '_4'
        '_4_5'
        '_6'
        '_6_10'

    Parameters
    ----------
    iterable : :class:`Iterable<collections.abc.Iterable>` [:class:`Iterable<collections.abc.Iterable>` [:data:`Any<typing.Any>`]]
        A nested iterable.

    Yields
    ------
    :class:`str`
        The accumulated nested elements of **iterable** as strings.

    """  # noqa
    for i in iterable:
        ret = ''
        for j in i:
            ret += f'_{j}'
            yield ret


def _get_opposite_neighbor(mol: Molecule,
                           idx_center: np.ndarray,
                           idx_neighbor: np.ndarray,
                           k: int = 4, n: int = 2,
                           **kwargs) -> np.ndarray:
    """Identify the **k** nearest neighbors of **idx_center** and return those furthest removed from each other."""  # noqa
    # Sanitize arguments
    xyz = np.asarray(mol)
    idx_center = np.array(idx_center, ndmin=2, copy=False)
    idx_neighbor = np.asarray(idx_neighbor)

    # Indices of the **k** nearest neighbors in **neighbor** with respect to **center**
    xyz1 = xyz[idx_neighbor]
    xyz2 = xyz[idx_center.ravel()]
    idx_nn = idx_neighbor[get_nearest_neighbors(xyz2, xyz1, k=k)]
    idx_nn.shape = -1, idx_nn.shape[1] * idx_center.shape[1]

    # Find the **n** atoms in **idx_nn** furthest removed from each other
    return brute_uniform_idx(xyz, idx_nn, n=n, **kwargs)


def _get_surface(mol: Molecule, symbol: str, max_dist: Optional[float] = None) -> np.ndarray:
    """Return the indices of all atoms, whose atomic symbol is equal to **atom_symbol**, located on the surface."""  # noqa
    # Identify all atom with atomic symbol **atom_symbol**
    atnum = to_atnum(symbol)
    idx = np.array([i for i, atom in enumerate(mol) if atom.atnum == atnum])
    xyz = np.asarray(mol)

    # Identify all atoms on the surface
    return idx[identify_surface(xyz[idx], max_dist=max_dist)]


def _mark_atoms(mol: Molecule, idx: Iterable[int]) -> None:
    """Mark all atoms in **mol** whose index is in **idx**; indices should be 1-based."""
    for i in idx:
        mol[i].properties.anchor = True
