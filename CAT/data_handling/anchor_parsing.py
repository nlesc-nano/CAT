"""A module for parsing the ``ligand.anchor`` keyword."""

import re
import operator
from typing import Union, Tuple, Iterable, SupportsFloat

from rdkit import Chem
from scm.plams import Units, PT
from schema import Schema, Use, Optional
from typing_extensions import TypedDict, SupportsIndex

from ..utils import AnchorTup, KindEnum, FormatEnum, MultiAnchorEnum
from ..attachment.ligand_anchoring import _smiles_to_rdmol, get_functional_groups

__all__ = ["parse_anchors"]


class _UnparsedAnchorDictBase(TypedDict):
    group: "str | SupportsIndex"
    anchor_idx: "SupportsIndex | Iterable[SupportsIndex]"


class _UnparsedAnchorDict(_UnparsedAnchorDictBase, total=False):
    remove: "None | SupportsIndex | Iterable[SupportsIndex]"
    angle_offset: "None | SupportsFloat | SupportsIndex | bytes | str"
    dihedral: "None | SupportsFloat | SupportsIndex | bytes | str"
    kind: "None | str | KindEnum"
    group_format: "None | str | FormatEnum"
    multi_anchor_filter: "None | str | MultiAnchorEnum"


class _AnchorDict(TypedDict):
    group: str
    group_idx: Tuple[int, ...]
    remove: "None | Tuple[int, ...]"
    kind: KindEnum
    angle_offset: "None | float"
    dihedral: "None | float"
    group_format: FormatEnum
    multi_anchor_filter: MultiAnchorEnum


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
    elif isinstance(typ, str):
        return KindEnum[typ.upper()]
    raise TypeError("`kind` expected None or a string")


def _parse_group_format(typ: "None | str | FormatEnum") -> FormatEnum:
    """Parse the ``group_format`` option."""
    if typ is None:
        return FormatEnum.SMILES
    elif isinstance(typ, FormatEnum):
        return typ
    elif isinstance(typ, str):
        return FormatEnum[typ.upper()]
    raise TypeError("`group_format` expected None or a string")


def _parse_multi_anchor_filter(typ: "None | str | MultiAnchorEnum") -> MultiAnchorEnum:
    """Parse the ``multi_anchor_filter`` option."""
    if typ is None:
        return MultiAnchorEnum.ALL
    elif isinstance(typ, MultiAnchorEnum):
        return typ
    elif isinstance(typ, str):
        return MultiAnchorEnum[typ.upper()]
    raise TypeError("`multi_anchor_filter` expected None or a string")


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


def _parse_group(group: "str | SupportsIndex") -> str:
    """Parse the ``group`` option."""
    # Pre-process in case the group is passed as an atomic number
    try:
        atnum = operator.index(group)
    except TypeError:
        pass
    else:
        group = PT.get_symbol(atnum)

    # String parsing
    if not isinstance(group, str):
        raise TypeError("`group` expected a string or integer")
    elif group in DUMMY_SYMBOLS:
        return "*"
    else:
        return group


def _symbol_to_rdmol(symbol: str) -> Chem.Mol:
    """Helper function for converting atomic symbols to rdkit molecules."""
    atom = Chem.Atom(PT.get_atomic_number(symbol))
    mol = Chem.EditableMol(Chem.Mol())
    mol.AddAtom(atom)
    return mol.GetMol()


anchor_schema = Schema({
    "group": Use(_parse_group),
    "group_idx": Use(_parse_group_idx),
    Optional("remove", default=None): Use(_parse_remove),
    Optional("kind", default=KindEnum.FIRST): Use(_parse_kind),
    Optional("angle_offset", default=None): Use(_parse_angle_offset),
    Optional("dihedral", default=None): Use(_parse_angle_offset),
    Optional("group_format", default=FormatEnum.SMILES): Use(_parse_group_format),
    Optional("multi_anchor_filter", default=MultiAnchorEnum.ALL): Use(_parse_multi_anchor_filter),
})

#: A collection of symbols used for different kinds of dummy atoms.
DUMMY_SYMBOLS = frozenset(PT.dummysymbols)

#: All atom types that are not directly supported by the rdkit SMILES parser.
INVALID_SMILES_ATOMS = frozenset(
    PT.symtonum.keys() - DUMMY_SYMBOLS - {'B', 'Br', 'C', 'Cl', 'F', 'I', 'N', 'O', 'P', 'S'}
)


def parse_anchors(
    patterns: Union[
        None,
        SupportsIndex,
        str,
        Chem.Mol,
        AnchorTup,
        _UnparsedAnchorDict,
        "Iterable[str | SupportsIndex | Chem.Mol | AnchorTup | _UnparsedAnchorDict]",
    ] = None,
    split: bool = True,
    is_core: bool = False,
) -> Tuple[AnchorTup, ...]:
    """Parse the user-specified anchors."""
    if patterns is None:
        if is_core:
            raise TypeError("`anchor=None` is not supported for core anchors")
        patterns = get_functional_groups(None, split)
    elif isinstance(patterns, (Chem.Mol, str, dict, AnchorTup, SupportsIndex)):
        patterns = [patterns]

    ret = []
    for p in patterns:  # type: _UnparsedAnchorDict | str | Chem.Mol | SupportsIndex | AnchorTup
        if isinstance(p, AnchorTup):
            ret.append(p)
        elif isinstance(p, Chem.Mol):
            mol = p
            remove: "None | Tuple[int, ...]" = (
                None if not split else (list(mol.GetAtoms())[-1].GetIdx(),)
            )
            ret.append(AnchorTup(mol=mol, remove=remove))
        elif isinstance(p, dict):
            kwargs: _AnchorDict = anchor_schema.validate(p)

            group_idx = kwargs["group_idx"]
            remove = kwargs["remove"]
            angle_offset = kwargs["angle_offset"]
            dihedral = kwargs["dihedral"]
            group_parser = kwargs["group_format"].value
            group = kwargs["group"]
            if group in INVALID_SMILES_ATOMS:
                mol = _symbol_to_rdmol(group)
            else:
                mol = group_parser(group)

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
            ret.append(AnchorTup(**kwargs, mol=mol))
        else:
            group = _parse_group(p)
            if group in INVALID_SMILES_ATOMS:
                mol = _symbol_to_rdmol(group)
            else:
                mol = _smiles_to_rdmol(group)
            remove = None if not split else (list(mol.GetAtoms())[-1].GetIdx(),)
            ret.append(AnchorTup(mol=mol, group=group, remove=remove))

    if is_core and len(ret) > 1:
        raise NotImplementedError("Cores with multiple anchor types aren't supported yet")
    return tuple(ret)
