"""
CAT.recipes.bulk
================

A short recipe for accessing the ligand-bulkiness workflow.

Index
-----
.. currentmodule:: CAT.recipes.bulk
.. autosummary::
    bulk_workflow

API
---
.. autofunction:: bulk_workflow

"""

from typing import Iterable, List, Iterator, Callable, Generator, Tuple
from itertools import chain

import numpy as np

from scm.plams import Molecule

from CAT.data_handling.mol_import import read_mol
from CAT.data_handling.validate_mol import validate_mol
from CAT.attachment.ligand_anchoring import _smiles_to_rdmol, find_substructure
from CAT.attachment.ligand_opt import optimize_ligand, allign_axis

try:
    from nanoCAT.mol_bulk import get_V, get_lig_radius
except ImportError as ex:
    tb = ex.__traceback__
    raise ImportError("Executing the content of '{__file__}' requires the Nano-CAT package: "
                      "'https://github.com/nlesc-nano/nano-CAT'").with_traceback(tb)

__all__ = ['bulk_workflow']


def bulk_workflow(smiles_list: Iterable[str],
                  optimize: bool = True) -> Tuple[List[Molecule], np.ndarray]:
    """Start the CAT ligand bulkiness workflow with an iterable of smiles strings.

    Examples
    --------
    .. code:: python

        >>> from CAT.recipes import bulk_workflow

        >>> smiles_list = [...]
        >>> mol_list, bulk_array = bulk_workflow(smiles_list, optimize=True)


    Parameters
    ----------
    smiles_list : :class:`Iterable<collections.abc.Iterable>` [:class:`str`]
        An iterable of SMILES strings.

    optimize : :class:`bool`
        Enable or disable the ligand geometry optimization.

    Returns
    -------
    :class:`list` [:class:`Molecule`] & :class:`numpy.ndarray`
        A list of plams Molecules and a matching array of :math:`V_{bulk}` values.

    """
    _mol_list = read_smiles(smiles_list)  # smiles to molecule
    mol_list = list(filter_mol(_mol_list))  # filter based on functional groups
    opt_and_allign(mol_list, opt=optimize)  # optimize and allign
    V_bulk = bulkiness(mol_list)  # calculate bulkiness
    return mol_list, V_bulk


def read_smiles(smiles_list: Iterable[str]) -> List[Molecule]:
    """Convert smiles strings into CAT-compatible plams molecules."""
    input_mol = list(smiles_list)
    validate_mol(input_mol, 'input_ligands')
    return read_mol(input_mol)


def filter_mol(mol_list: Iterable[Molecule], anchor: str = 'O(C=O)[H]') -> Iterator[Molecule]:
    """Filter all input molecules based on the presence of a functional group (the "anchor")."""
    anchor_rdmols = (_smiles_to_rdmol(anchor),)
    return chain.from_iterable(find_substructure(mol, anchor_rdmols) for mol in mol_list)


def opt_and_allign(mol_list: Iterable[Molecule], opt: bool = True) -> None:
    """Optimize all molecules and allign them along the x-axis; set :code:`opt=False` to disable the optimization."""  # noqa
    def _allign_axis(mol: Molecule) -> None:
        return allign_axis(mol, mol.properties.dummies)

    process_mol: Callable[[Molecule], None] = optimize_ligand if opt else _allign_axis
    for mol in mol_list:
        process_mol(mol)


def bulkiness(mol_list: Iterable[Molecule], diameter: float = 4.5,
              height_lim: float = 10.0) -> np.ndarray:
    r"""Calculate the ligand bulkiness descriptor :math:`V_{bulk}`.

    .. math::

        V(d, h_{lim}) =
        \sum_{i=1}^{n} e^{r_{i}} (\frac{2 r_{i}}{d} - 1)^{+} (1 - \frac{h_{i}}{h_{lim}})^{+}

    """
    def _iterate(mol_list: Iterable[Molecule]) -> Generator[float, None, None]:
        for mol in mol_list:
            radius, height = get_lig_radius(mol)  # From cartesian to cylindrical
            yield get_V(radius, height, diameter, None, h_lim=height_lim)

    try:
        count = len(mol_list)
    except TypeError:
        count = -1
    return np.fromiter(_iterate(mol_list), count=count, dtype=float)
