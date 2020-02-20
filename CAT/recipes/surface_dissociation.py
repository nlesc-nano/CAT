"""
CAT.recipes.surface_dissociation
================================

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

from typing import Tuple, Dict, Sequence, Iterable, Optional
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from scm.plams import Molecule

try:
    from nanoCAT.bde.dissociate_xyn import dissociate_ligand
    from nanoCAT.bde.identify_surface import identify_surface
except ImportError as ex:
    tb = ex.__traceback__
    raise ImportError("Executing the content of '{__file__}' requires the Nano-CAT package: "
                      "'https://github.com/nlesc-nano/nano-CAT'").with_traceback(tb)

__all__ = ['dissociate_surface']

ArrayDict = Dict[Tuple[int, int], np.ndarray]


def dissociate_surface():
    # Read the .xyz file
    base_path = Path('path/to/my/qd/directory')
    quantum_dot = Molecule(base_path / 'fancy_qd.xyz')

    # Put the indices of two opposing Cs atoms right here
    Cs_array = np.array([
        [1, 3],
        [4, 5],
        [6, 10],
        [15, 12],
        [99, 105],
        [20, 4]
    ])
    Cs_array.sort(axis=1)
    Cs_array = Cs_array[:, ::-1]  # Sort in reverse order, this is important

    # Identify all Cl atoms on the surface
    Cl_array = get_surface(quantum_dot, atom_symbol='Cl')

    # Construct a dictionary with all 2*Cs & 2*Cl pairs
    Cs_Cl_dict = get_neighbor_dict(quantum_dot, Cs_array-1, Cl_array)

    for Cs_pair, Cl_pair in Cs_Cl_dict.items():
        filename = 'output'

        mol = quantum_dot.copy()
        mark_atoms(mol, Cl_pair)

        for Cs_index in Cs_pair:
            mol = next(dissociate_ligand(mol, lig_count=1,
                                        core_index=Cs_index,
                                        lig_core_pairs=1))
            filename += f'_{Cs_index}'
            mol.write(base_path / f'{filename}.xyz')


def get_neighbor_dict(mol: Molecule,
                      Cs_idx: Sequence[Tuple[int, int]],
                      Cl_idx: Sequence[int],
                      k: int = 4) -> ArrayDict:
    """Create a dictionary with Cs-pairs as keys and Cl-pairs as values."""
    # Sanitize arguments
    xyz = np.asarray(mol)
    Cs_idx = np.array(Cs_idx).ravel()
    Cl_idx = np.asarray(Cl_idx)
    k2 = 2 * k

    # Find the k Cl atoms closest to each Cs atom
    tree = cKDTree(xyz[Cl_idx])
    _, Cl_Cl_idx = tree.query(xyz[Cs_idx], k=k)
    Cl_Cl_idx.shape = -1, k2

    # Evaluate all distances within the Cl-neighbor subset
    xyz_tensor = xyz[Cl_idx[Cl_Cl_idx]]
    dist = np.linalg.norm(xyz_tensor[..., None, :] - xyz_tensor[:, None, ...], axis=-1)
    dist.shape = -1, k2**2

    # Find the index-pair yielding the maximum Cl-Cl dist per Cs atom
    i = np.repeat(np.arange(len(dist)), 2)
    j = np.ravel(np.unravel_index(dist.argmax(axis=1), shape=(k2, k2)))

    # Parse and return the indices
    idx_ret = Cl_idx[Cl_Cl_idx[i, j].T]
    idx_ret.shape = -1, 2
    Cs_idx.shape = -1, 2
    idx_ret += 1  # Switch the 1-based indices
    Cs_idx += 1
    return {tuple(pair1): pair2 for pair1, pair2 in zip(Cs_idx, idx_ret)}


def get_surface(mol: Molecule, atom_symbol: str, max_dist: Optional[float] = None) -> np.ndarray:
    """Return the (0-based) indices of all atoms, whose atomic symbol is equal to **atom_symbol**, located on the surface."""  # noqa
    # Identify all atom with atomic symbol **atom_symbol**
    idx_superset = np.array([i for i, atom in enumerate(mol) if atom.symbol == atom_symbol])
    xyz = np.asarray(mol)[idx_superset]

    # Identify all atoms on the surface
    return idx_superset[identify_surface(xyz, max_dist=max_dist)]


def mark_atoms(mol: Molecule, idx: Iterable[int]) -> None:
    """Mark all atoms in **mol** whose index is in **idx**."""
    for i in idx:
        mol[i].properties.anchor = True
