"""A module for parsing the ``ligand.anchor`` keyword."""

import operator
from typing import Union, Tuple, Collection, Iterable

from rdkit.Chem import Mol
from schema import Schema, Use, Optional
from typing_extensions import TypedDict, SupportsIndex

from ..utils import AnchorTup, KindEnum
from ..attachment.ligand_anchoring import _smiles_to_rdmol, get_functional_groups

__all__ = ["parse_anchors"]


class _UnparsedAnchorDictBase(TypedDict):
    group: str
    anchor_idx: "SupportsIndex | Collection[SupportsIndex]"


class _UnparsedAnchorDict(_UnparsedAnchorDictBase, total=False):
    remove: "None | SupportsIndex | Collection[SupportsIndex]"


def _parse_group_idx(item: "SupportsIndex | Iterable[SupportsIndex]") -> Tuple[int, ...]:
    """Parse the ``group_idx`` option."""
    try:
        return (operator.index(item),)
    except TypeError:
        pass

    ret = tuple(operator.index(i) for i in item)
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


anchor_schema = Schema({
    "group": str,
    "group_idx": Use(lambda i: tuple(_parse_group_idx(i))),
    Optional("remove", default=None): Use(_parse_remove),
    Optional("kind", default=KindEnum.FIRST): Use(_parse_kind),
})


def parse_anchors(
    patterns: Union[
        None,
        str,
        Mol,
        AnchorTup,
        _UnparsedAnchorDict,
        "Collection[str | Mol | AnchorTup | _UnparsedAnchorDict]",
    ] = None,
    split: bool = True,
) -> Tuple[AnchorTup, ...]:
    """Parse the user-specified anchors."""
    if patterns is None:
        patterns = get_functional_groups(None, split)
    elif isinstance(patterns, (Mol, str, dict)):
        patterns = [patterns]

    ret = []
    for p in patterns:  # type: _UnparsedAnchorDict | str | Mol
        if isinstance(p, AnchorTup):
            ret.append(p)
        elif isinstance(p, Mol):
            mol = p
            remove = None if not split else (list(mol.GetAtoms())[-1].GetIdx(),)
            ret.append(AnchorTup(mol=mol, remove=remove))
        elif isinstance(p, str):
            group = p
            mol = _smiles_to_rdmol(group)
            remove = None if not split else (list(mol.GetAtoms())[-1].GetIdx(),)
            ret.append(AnchorTup(mol=mol, group=group, remove=remove))
        else:
            kwargs = anchor_schema.validate(p)
            group_idx = kwargs["group_idx"]
            remove = kwargs["remove"]
            if remove is not None and not set(group_idx).isdisjoint(remove):
                raise ValueError("`group_idx` and `remove` must be disjoint")

            mol = _smiles_to_rdmol(kwargs["group"])
            ret.append(AnchorTup(**kwargs, mol=mol))
    return tuple(ret)
