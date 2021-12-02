"""A module for parsing the ``ligand.anchor`` keyword."""

import re
import operator
from typing import Union, Tuple, Iterable, SupportsFloat

from rdkit.Chem import Mol
from scm.plams import Units, PT
from schema import Schema, Use, Optional
from typing_extensions import TypedDict, SupportsIndex

from ..utils import AnchorTup, KindEnum
from ..attachment.ligand_anchoring import _smiles_to_rdmol, get_functional_groups

__all__ = ["parse_anchors"]


class _UnparsedAnchorDictBase(TypedDict):
    group: str
    anchor_idx: "SupportsIndex | Iterable[SupportsIndex]"


class _UnparsedAnchorDict(_UnparsedAnchorDictBase, total=False):
    remove: "None | SupportsIndex | Iterable[SupportsIndex]"
    angle_offset: "None | SupportsFloat | SupportsIndex | bytes | str"
    dihedral: "None | SupportsFloat | SupportsIndex | bytes | str"


class _AnchorDict(TypedDict):
    group: str
    group_idx: Tuple[int, ...]
    remove: "None | Tuple[int, ...]"
    kind: KindEnum
    angle_offset: "None | float"
    dihedral: "None | float"


def _parse_group_idx(item: "SupportsIndex | Iterable[SupportsIndex]") -> Tuple[int, ...]:
    """Parse the ``group_idx`` option."""
    try:
        return (operator.index(item),)
    except TypeError:
        pass

    ret = tuple(operator.index(i) for i in item)  # type: ignore[union-attr]
    n = len(ret) - len(set(ret))
    if n:
        raise ValueError(f"Found {n} duplicate elements")
    elif not ret:
        raise ValueError("Requires at least one element")
    return ret


def _parse_remove(
    item: "None | SupportsIndex | Iterable[SupportsIndex]"
) -> "None | Tuple[int, ...]":
    """Parse the ``remove`` option."""
    if item is None:
        return None
    return _parse_group_idx(item)


def _parse_kind(typ: "None | str | KindEnum") -> KindEnum:
    """Parse the ``kind`` option."""
    if typ is None:
        return KindEnum.FIRST
    elif isinstance(typ, KindEnum):
        return typ
    else:
        return KindEnum[typ.upper()]


_UNIT_PATTERN = re.compile(r"([\.\_0-9]+)(\s+)?(\w+)?")


def _parse_angle_offset(
    offset: "None | SupportsFloat | SupportsIndex | bytes | str"
) -> "None | float":
    """Parse the ``angle_offset`` and ``dihedral`` options; convert the offset to radians."""
    if offset is None:
        return None
    elif not isinstance(offset, str):
        return Units.convert(float(offset), "deg", "rad")

    # match offsets such as `"5.5 degree"`
    match = _UNIT_PATTERN.match(offset)
    if match is None:
        raise ValueError(f"Invalid offset string: {offset!r}")

    offset, _, unit = match.groups()
    if unit is None:
        unit = "deg"
    return Units.convert(float(offset), unit, "rad")


anchor_schema = Schema({
    "group": str,
    "group_idx": Use(_parse_group_idx),
    Optional("remove", default=None): Use(_parse_remove),
    Optional("kind", default=KindEnum.FIRST): Use(_parse_kind),
    Optional("angle_offset", default=None): Use(_parse_angle_offset),
    Optional("dihedral", default=None): Use(_parse_angle_offset)
})


#: All atom types that have to be encapsulated in square brackets when parsing SMILES strings
SQUARE_BRACKET_ATOMS = frozenset(
    PT.symtonum.keys() - {'B', 'Br', 'C', 'Cl', 'F', 'I', 'N', 'O', 'P', 'S'}
)


def parse_anchors(
    patterns: Union[
        None,
        SupportsIndex,
        str,
        Mol,
        AnchorTup,
        _UnparsedAnchorDict,
        "Iterable[str | SupportsIndex | Mol | AnchorTup | _UnparsedAnchorDict]",
    ] = None,
    split: bool = True,
    is_core: bool = False,
) -> Tuple[AnchorTup, ...]:
    """Parse the user-specified anchors."""
    if patterns is None:
        if is_core:
            raise TypeError("`anchor=None` is not supported for core anchors")
        patterns = get_functional_groups(None, split)
    elif isinstance(patterns, (Mol, str, dict, AnchorTup, SupportsIndex)):
        patterns = [patterns]

    ret = []
    for p in patterns:  # type: _UnparsedAnchorDict | str | Mol | SupportsIndex | AnchorTup
        try:
            atnum = operator.index(p)  # Check for atomic symbols
        except TypeError:
            pass
        else:
            p = PT.get_symbol(atnum)

        if isinstance(p, AnchorTup):
            ret.append(p)
        elif isinstance(p, Mol):
            mol = p
            remove: "None | Tuple[int, ...]" = (
                None if not split else (list(mol.GetAtoms())[-1].GetIdx(),)
            )
            ret.append(AnchorTup(mol=mol, remove=remove))
        elif isinstance(p, str):
            group = p
            if group in SQUARE_BRACKET_ATOMS:
                group = f"[{group}]"
            mol = _smiles_to_rdmol(group)
            remove = None if not split else (list(mol.GetAtoms())[-1].GetIdx(),)
            ret.append(AnchorTup(mol=mol, group=group, remove=remove))
        else:
            kwargs: _AnchorDict = anchor_schema.validate(p)

            group_idx = kwargs["group_idx"]
            remove = kwargs["remove"]
            angle_offset = kwargs["angle_offset"]
            dihedral = kwargs["dihedral"]

            group = kwargs.pop("group")
            if group in SQUARE_BRACKET_ATOMS:
                group = f"[{group}]"
            mol = _smiles_to_rdmol(group)

            # Dihedral and angle-offset options are not supported for core anchors
            if is_core:
                if dihedral is not None:
                    raise TypeError("`dihedral != None` is not supported for core anchors")
                elif angle_offset is not None:
                    raise TypeError("`angle_offset != None` is not supported for core anchors")
                elif kwargs["kind"] != KindEnum.FIRST:
                    raise NotImplementedError('`kind != "first"` is not yet supported')
            else:
                # Check that at least 3 atoms are available for `angle_offset`
                # (so a plane can be defined)
                if angle_offset is not None and len(group_idx) < 3:
                    raise ValueError("`group_idx` must contain at least 3 atoms when "
                                     "`angle_offset` is specified")

                # Check that at least 2 atoms are available for `dihedral`
                # (so the third dihedral-defining vector can be defined)
                if dihedral is not None and len(group_idx) < 2:
                    raise ValueError("`group_idx` must contain at least 2 atoms when "
                                     "`dihedral` is specified")

                # Check that `group_idx` and `remove` are disjoint
                # TODO: Investigate if this check can be removed
                if remove is not None and not set(group_idx).isdisjoint(remove):
                    raise ValueError("`group_idx` and `remove` must be disjoint")

            # Check that the indices in `group_idx` and `remove` are not out of bounds
            atom_count = len(mol.GetAtoms())
            if atom_count <= max(group_idx):
                raise IndexError(f"`group_idx` index {max(group_idx)} is out of bounds "
                                 f"for a `group` with {atom_count} atoms")
            elif remove is not None and atom_count <= max(remove):
                raise IndexError(f"`remove` index {max(remove)} is out of bounds "
                                 f"for a `group` with {atom_count} atoms")
            ret.append(AnchorTup(**kwargs, group=group, mol=mol))
    if is_core and len(ret) > 1:
        raise NotImplementedError("Cores with multiple anchor types aren't supported yet")
    return tuple(ret)
