"""A module for parsing the ``ligand.anchor`` keyword."""

import re
import enum
import operator
import pprint
from itertools import chain
from typing import Union, NamedTuple, Tuple, Collection, Iterable, Iterator

from scm import plams
from rdkit.Chem import Mol
from schema import Schema, Use, Optional
from typing_extensions import TypedDict, SupportsIndex

from ..attachment.ligand_anchoring import _smiles_to_rdmol, get_functional_groups

__all__ = ["parse_anchors"]

_ATOMS = "|".join(sorted(
    chain(plams.PT.symtonum, plams.PT.dummysymbols),
    key=lambda i: (len(i), i),
))
PATTERN = re.compile(f"({_ATOMS})")


class _UnparsedAnchorDictBase(TypedDict):
    group: str
    anchor_idx: "SupportsIndex | Collection[SupportsIndex]"


class _UnparsedAnchorDict(_UnparsedAnchorDictBase, total=False):
    remove: "None | SupportsIndex | Collection[SupportsIndex]"


class KindEnum(enum.Enum):
    FIRST = 0
    MEAN = 1


class AnchorTup(NamedTuple):
    mol: Mol
    anchor_idx: Tuple[int, ...] = (0,)
    atoms: "None | Tuple[str, ...]" = None
    group: "None | str" = None
    remove: "None | Tuple[int, ...]" = None
    kind: KindEnum = KindEnum.FIRST


class MolMatches:
    __slots__ = ("__weakref__", "anchor_tup")

    @property
    def mol(self) -> Iterator[Mol]:
        return (i.mol for i in self.anchor_tup)

    @property
    def group(self) -> Iterator["None | str"]:
        return (i.group for i in self.anchor_tup)

    @property
    def atoms(self) -> Iterator["None | Tuple[str, ...]"]:
        return (i.atoms for i in self.anchor_tup)

    @property
    def anchor_idx(self) -> Iterator["Tuple[int, ...]"]:
        return (i.anchor_idx for i in self.anchor_tup)

    @property
    def remove(self) -> Iterator["None | Tuple[int, ...]"]:
        return (i.remove for i in self.anchor_tup)

    @property
    def kind(self) -> Iterator[KindEnum]:
        return (i.kind for i in self.anchor_tup)

    def __init__(self, anchor_tups: Iterable[AnchorTup]) -> None:
        self.anchor_tup = tuple(anchor_tups)

    def __repr__(self) -> str:
        name = type(self).__name__
        width = 80 - len(name)
        indent = len(name) + 2
        values = pprint.pformat(self.anchor_tup, width=width, indent=indent)[indent:]
        return f"{name}(({values})"

    def get_matches(mol: "plams.Molecule | Mol") -> None:
        mol = plams.to_rdmol(mol)


def _parse_anchor_idx(item: "SupportsIndex | Iterable[SupportsIndex]") -> Tuple[int, ...]:
    """Parse the ``anchor_idx`` option."""
    try:
        return (operator.index(item),)
    except TypeError:
        pass

    ret = tuple(operator.index(i) for i in item)
    n = len(ret) - len(set(ret))
    if n:
        raise ValueError(f"Found {n} duplicate elements")
    return ret


def _parse_remove(
    item: "None | SupportsIndex | Iterable[SupportsIndex]"
) -> "None | Tuple[int, ...]":
    """Parse the ``remove`` option."""
    if item is None:
        return None
    else:
        return _parse_anchor_idx(item)


def _parse_kind(typ: "None | str") -> KindEnum:
    """Parse the ``kind`` option."""
    if typ is None:
        return KindEnum.FIRST
    typ = typ.upper()
    return KindEnum[typ]


anchor_schema = Schema({
    "group": str,
    "anchor_idx": Use(_parse_anchor_idx),
    Optional("remove", default=None): Use(_parse_remove),
    Optional("kind", default=KindEnum.FIRST): Use(_parse_kind),
})


def parse_anchors(
    patterns: Union[
        None,
        str,
        Mol,
        MolMatches,
        _UnparsedAnchorDict,
        Collection[str],
        Collection[Mol],
        Collection[_UnparsedAnchorDict],
    ] = None,
    split: bool = True,
) -> MolMatches:
    """Parse the user-specified anchors."""
    if isinstance(patterns, MolMatches):
        return patterns
    elif patterns is None:
        return MolMatches(AnchorTup(m) for m in get_functional_groups(None, split))
    if isinstance(patterns, (Mol, str, dict)):
        patterns = [patterns]

    ret = []
    for _p in patterns:  # type: _UnparsedAnchorDict | str | Mol
        if isinstance(_p, Mol):
            ret.append(AnchorTup(mol=_p))
            continue
        elif isinstance(_p, str):
            group = _p
            mol = _smiles_to_rdmol(group)
            atoms = PATTERN.findall(group)
            remove = None if not split else list(mol.GetAtoms())[-1].GetIdx()
            if atoms is None:
                raise ValueError(f"Failed to extract atomic symbols from {group}")
            ret.append(AnchorTup(mol=mol, group=group, atoms=tuple(atoms), remove=remove))
            continue

        p = anchor_schema.validate(_p)
        group = p["group"]
        atoms = PATTERN.findall(group)
        if atoms is None:
            raise ValueError(f"Failed to extract atomic symbols from {group}")
        p["atoms"] = tuple(atoms)
        p["mol"] = _smiles_to_rdmol(group)
        ret.append(AnchorTup(**p))
    return MolMatches(ret)
