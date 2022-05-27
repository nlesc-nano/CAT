"""A module designed for finding ligand functional groups.

Index
-----
.. currentmodule:: CAT.attachment.ligand_anchoring
.. autosummary::
    init_ligand_anchoring
    get_functional_groups
    _smiles_to_rdmol
    find_substructure
    substructure_split
    _get_df

API
---
.. autofunction:: init_ligand_anchoring
.. autofunction:: get_functional_groups
.. autofunction:: find_substructure
.. autofunction:: _smiles_to_rdmol
.. autofunction:: substructure_split
.. autofunction:: _get_df

"""

from itertools import chain, repeat
from collections import defaultdict
from typing import Sequence, List, Tuple, Iterable, Callable, Iterator

import pandas as pd

from scm.plams import Molecule, Settings, MoleculeError
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..logger import logger
from ..utils import get_template, AnchorTup, KindEnum, get_formula, FormatEnum, MultiAnchorEnum
from ..mol_utils import separate_mod   # noqa: F401
from ..workflows import MOL, FORMULA, HDF5_INDEX, OPT
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.validate_mol import santize_smiles

__all__ = ['init_ligand_anchoring']


def init_ligand_anchoring(ligand_df: SettingsDataFrame) -> SettingsDataFrame:
    """Initialize the ligand functional group searcher.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of valid ligands.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of ligands with functional groups that can serve as valid anchor points.

    """
    # Unpack arguments
    settings = ligand_df.settings.optional
    functional_groups = settings.ligand.anchor

    # Find all functional groups; return a copy of each mol for each functional group
    mol_list = []
    for lig in ligand_df[MOL]:
        # Functional group search
        dummies = lig.properties.dummies
        if dummies is None:
            mol_list += find_substructure(lig, functional_groups)
            continue
        else:
            dummies = tuple(i - 1 for i in dummies)

        # Manual specification of a functional group
        if len(dummies) == 1:  # optional.ligand.split = False
            remove = None
        elif len(dummies) == 2:  # optional.ligand.split = True
            remove = (1,)
        else:
            raise NotImplementedError
        anchor_tup = AnchorTup(None, group_idx=(0,), remove=remove)
        mol_list += [substructure_split(lig, dummies, anchor_tup)]

    # Convert the results into a dataframe
    return _get_df(mol_list, ligand_df.settings)


def _get_df(
    mol_list: Sequence[Molecule],
    settings: Settings,
) -> SettingsDataFrame:
    """Create and return a new ligand dataframe.

    Parameters
    ----------
    mol_list : |list|_ [|plams.Molecule|_]
        A list of PLAMS molecules.
    settings : |Settings|_
        A Settings instance containing all CAT parameters.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of ligands with functional groups that can serve as valid anchor points.

    """
    # Create the dataframe index and columns
    idx_tuples = [(mol.properties.smiles, mol.properties.anchor) for mol in mol_list]
    idx = pd.MultiIndex.from_tuples(idx_tuples, names=['smiles', 'anchor'])
    columns_tuples = [MOL, FORMULA, HDF5_INDEX, OPT]
    columns = pd.MultiIndex.from_tuples(columns_tuples, names=['index', 'sub index'])

    # Create, fill and return the dataframe
    df = SettingsDataFrame(-1, index=idx, columns=columns, settings=settings)
    df[MOL] = mol_list
    df[FORMULA] = [get_formula(lig) for lig in df[MOL]]
    df[OPT] = False
    return df[~df.index.duplicated(keep='first')]  # Remove duplicate indices


def get_functional_groups(
    functional_groups: "None | Iterable[str]" = None,
    split: bool = True,
) -> Tuple[Chem.Mol, ...]:
    """Construct a list of RDKit molecules representing functional groups.

    Parameters
    ----------
    functional_groups : |list|_ [|str|_]
        Optional: A list of smiles strings representing functional groups.
        Will default to templates provided by CAT if ``None``.
    split : bool
        If templates should be pulled from the ``['split']`` or ``['no_split']`` block.
        Only relevant if **functional_groups** is ``None``.

    Returns
    -------
    |tuple|_ [|Chem.Mol|_]
        A list of RDKit molecules constructed from either **functional_group** or
        the default smiles templates in CAT.

    """
    # The user has, explicitly, provided functional groups
    if functional_groups:
        return tuple(_smiles_to_rdmol(smiles) for smiles in functional_groups)

    # Read functional groups from the default CAT SMILES templates
    if split:
        func_groups = get_template('smiles.yaml').split
    else:
        func_groups = get_template('smiles.yaml').no_split

    return tuple(_smiles_to_rdmol(smiles) for smiles in func_groups)


_smiles_to_rdmol = FormatEnum.SMILES.value


def find_substructure(
    ligand: Molecule,
    func_groups: Iterable[AnchorTup],
    split: bool = True,
    condition: "None | Callable[[int], bool]" = None,
) -> List[Molecule]:
    """Identify interesting functional groups within the ligand.

    Parameters
    ----------
    ligand : |plams.Molecule|_
        The ligand molecule.
    func_groups : |tuple|_ [|Chem.Mol|_]
        A collection of RDKit molecules representing functional groups.

    Returns
    -------
    |list|_ [|plams.Molecule|_]
        A list of ligands.
        A single copy of **ligand** is created for each identified functional group,
        removing parts of the functional group if required (see **split**).
        An empty list is returned if no valid functional groups are found.

    """
    rdmol = molkit.to_rdmol(ligand)

    # Searches for functional groups (defined by functional_group_list) within the ligand
    matches: Iterator[Tuple[AnchorTup, Tuple[int, ...]]] = chain.from_iterable(
        zip(repeat(anchor_tup), rdmol.GetSubstructMatches(anchor_tup.mol, useChirality=True))
        for anchor_tup in func_groups
    )

    # Remove all duplicate matches, each heteroatom (match[0]) should have <= 1 entry
    ligand_idx_dict = defaultdict(list)
    ref_dict = defaultdict(set)
    for anchor_tup, idx_tup in matches:
        ref_set = ref_dict[anchor_tup]
        anchor_idx_tup = tuple(idx_tup[i] for i in anchor_tup.group_idx)
        if anchor_idx_tup in ref_set:
            continue  # Skip duplicates

        ligand_idx_dict[anchor_tup].append(idx_tup)
        ref_set.add(anchor_idx_tup)

    # Apply some further filtering to the ligands
    if condition is not None:
        if not condition(sum((len(j) for j in ligand_idx_dict.values()), 0)):
            err = (f"Failed to satisfy the passed condition ({condition!r}) for "
                   f"ligand: {ligand.properties.name!r}")
            logger.error(err)
            return []
    else:
        for anchor_tup, j in ligand_idx_dict.items():
            if anchor_tup.multi_anchor_filter == MultiAnchorEnum.ALL:
                pass
            elif anchor_tup.multi_anchor_filter == MultiAnchorEnum.FIRST and len(j) > 1:
                ligand_idx_dict[anchor_tup] = j[:1]
            elif anchor_tup.multi_anchor_filter == MultiAnchorEnum.RAISE and len(j) > 1:
                logger.error(
                    f"Found multiple valid functional groups for {ligand.properties.name!r}"
                )
                return []

    ret = []
    idx_dict_items = chain.from_iterable(zip(repeat(k), v) for k, v in ligand_idx_dict.items())
    for i, (anchor_tup, idx_tup) in enumerate(idx_dict_items):
        try:
            value = substructure_split(ligand, idx_tup, anchor_tup)
        except Exception as ex:
            logger.warning(
                f"Failed to parse {ligand.properties.name!r} anchoring group {i}", exc_info=ex
            )
        else:
            ret.append(value)

    if not ret:
        err = (f"No functional groups were found (optional.ligand.split = {split!r}) for "
               f"ligand: {ligand.properties.name!r}")
        logger.error(err)
    return ret


def substructure_split(
    ligand: Molecule,
    idx_tup: Tuple[int, ...],
    anchor_tup: AnchorTup,
) -> Molecule:
    """Delete the hydrogen or mono-/polyatomic counterion attached to the functional group.

    Sets the charge of the remaining heteroatom to -1 if ``split=True``.

    Parameters
    ----------
    ligand: plams.Molecule
        The ligand molecule.
    idx_tup : tuple[int, ...]
        A tuple with 2 atomic indices associated with a functional group.
    anchor_tup : AnchorTup
        Named tuple.

    Returns
    -------
    plams.Molecule
        A copy of **ligand**, with part of its functional group removed (see **split**).

    """
    lig = ligand.copy()
    anchors = [lig[1 + idx_tup[i]] for i in anchor_tup.group_idx]
    anchor = anchors[0]

    if anchor_tup.remove is not None:
        remove_iter = sorted([1 + idx_tup[i] for i in anchor_tup.remove], reverse=True)
        for i in remove_iter:
            lig.delete_atom(lig[i])

        mol_list: Tuple[Molecule, ...] = lig.separate_mod()
        for mol in mol_list:
            if anchor not in mol:
                continue
            lig = mol
            break
        else:
            raise MoleculeError("Failed to identify an anchor-containing fragment")

        # Check if the ligand heteroatom has a charge assigned, assigns a charge if not
        if not anchor.properties.charge:
            anchor.properties.charge = -1

    anchor_tup = anchor_tup._replace(anchor_idx=tuple(lig.atoms.index(at) for at in anchors))

    # Update ligand properties
    lig.properties.dummies = anchor
    if anchor_tup.kind == KindEnum.FIRST:
        i = 1 + anchor_tup.anchor_idx[0]
        lig.properties.anchor = f"{lig[i].symbol}{i}"
    else:
        lig.properties.anchor = "".join(
            f"{lig[1 + i].symbol}{1 + i}" for i in anchor_tup.anchor_idx
        )
    lig.properties.anchor_tup = anchor_tup
    lig.properties.charge = sum(atom.properties.get('charge', 0) for atom in lig)

    # Update the ligand smiles string
    rdmol = molkit.to_rdmol(lig)
    smiles = Chem.MolToSmiles(rdmol)  # NOTE: can return `None` under rare circumstances
    lig.properties.smiles = Chem.CanonSmiles(smiles)
    lig.properties.name = santize_smiles(lig.properties.smiles) + '@' + lig.properties.anchor
    lig.properties.path = ligand.properties.path
    return lig
