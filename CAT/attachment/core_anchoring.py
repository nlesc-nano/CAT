"""A module designed for finding core functional groups.

Index
-----
.. currentmodule:: CAT.attachment.ligand_anchoring
.. autosummary::
    set_core_anchors
    find_core_substructure

API
---
.. autofunction:: set_core_anchors
.. autofunction:: find_core_substructure

"""

from typing import Tuple, Any, Mapping, Iterable, TYPE_CHECKING

import numpy as np
from scm.plams import Molecule, Atom, MoleculeError, to_rdmol

from .perp_surface import get_surface_vec
from .distribution import distribute_idx
from ..utils import AllignmentEnum, AllignmentTup, AnchorTup

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from numpy import int64 as i8, float64 as f8

__all__ = ["set_core_anchors", "find_core_substructure"]


def set_core_anchors(
    mol: Molecule,
    anchor_tup: AnchorTup,
    allignment_tup: AllignmentTup,
    subset_kwargs: "None | Mapping[str, Any]" = None,
) -> Tuple[str, str]:
    """Identify and parse the core anchors within the passed molecule.

    Returns two strings: The (parsed) molecular formula and the anchor indices.

    """
    # Checks the if the anchor is a string (atomic symbol) or integer (atomic number)
    formula = mol.get_formula()

    # Get the indices of all anchor atom ligand placeholders in the core
    anchors = mol.properties.dummies
    if anchors is None:
        anchor_idx, remove_idx = find_core_substructure(mol, anchor_tup)
    else:
        anchor_idx = np.fromiter(anchors, count=len(anchors), dtype=np.int64)
        anchor_idx -= 1
        remove_idx = anchor_idx.copy()
    if subset_kwargs:
        anchor_idx = distribute_idx(mol, anchor_idx, **subset_kwargs)
    if not len(anchor_idx):
        raise MoleculeError(f"No valid anchoring groups found in the core {formula!r}")

    # Convert atomic indices into Atoms
    anchor_idx += 1
    anchor_idx.sort()
    mol.properties.dummies = [mol[i] for i in anchor_idx]

    # Returns an error if no anchor atoms were found
    if (len(mol) - len(anchor_idx)) < 4 and allignment_tup.kind == AllignmentEnum.SURFACE:
        raise ValueError(
            '`optional.core.allignment = "surface"` is not supported for cores with less '
            f'than 4 (non-anchor) atoms ({mol.get_formula()}); consider using '
            '`optional.core.allignment = "sphere"`'
        )
    elif len(anchor_tup.mol.GetAtoms()) < 2 and allignment_tup.kind == AllignmentEnum.ANCHOR:
        raise ValueError(
            '`optional.core.allignment = "anchor"` is not supported for mono-atomic core anchors'
        )

    # Define all core vectors
    mol.properties.core_vec = _get_core_vectors(mol, mol.properties.dummies, allignment_tup)

    # Delete all core anchor atoms
    if remove_idx is not None:
        remove_idx += 1
        remove_idx.sort()
        for i in reversed(remove_idx):
            mol.delete_atom(mol[i])
    return formula, ' '.join(anchor_idx.astype(str))


def _get_core_vectors(
    core: Molecule,
    dummies: Iterable[Atom],
    allignment: AllignmentTup,
) -> "NDArray[f8]":
    """Return a 2D array with all core (unit) vectors."""
    anchor = Molecule.as_array(None, atom_subset=dummies)

    if allignment.kind == AllignmentEnum.SPHERE:
        vec = np.array(core.get_center_of_mass()) - anchor
        vec /= np.linalg.norm(vec, axis=1)[..., None]
    elif allignment.kind == AllignmentEnum.SURFACE:
        vec = -get_surface_vec(np.array(core), anchor)
    elif allignment.kind == AllignmentEnum.ANCHOR:
        raise NotImplementedError
    else:
        raise ValueError(f"Unknown allignment kind: {allignment.kind}")

    if allignment.invert:
        vec *= -1
    return vec


def find_core_substructure(
    mol: Molecule,
    anchor_tup: AnchorTup,
) -> Tuple["NDArray[i8]", "None | NDArray[i8]"]:
    """Identify substructures within the passed core based on **anchor_tup**.

    Returns two indice-arrays, respectivelly containing the indices of the anchor
    atoms and all to-be removed atoms.

    """
    rdmol = to_rdmol(mol)
    matches = rdmol.GetSubstructMatches(anchor_tup.mol, useChirality=True)
    remove = anchor_tup.remove

    # Remove all duplicate matches, each heteroatom (match[0]) should have <= 1 entry
    ref_set = set()
    anchor_list = []
    remove_list = []
    for idx_tup in matches:
        anchor_idx_tup = tuple(idx_tup[i] for i in anchor_tup.group_idx)
        if anchor_idx_tup in ref_set:
            continue  # Skip duplicates
        else:
            ref_set.add(anchor_idx_tup)

        if remove is not None:
            remove_list += [idx_tup[i] for i in remove]
        anchor_list.append(anchor_idx_tup[0])

    anchor_array = np.fromiter(anchor_list, dtype=np.int64, count=len(anchor_list))
    if remove is not None:
        remove_array = np.fromiter(remove_list, dtype=np.int64, count=len(remove_list))
        return anchor_array, remove_array
    else:
        return anchor_array, None
