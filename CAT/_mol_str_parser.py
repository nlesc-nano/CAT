"""A module with various RDKit-based string-to-mol parser.

Index
-----
.. currentmodule:: CAT._mol_str_parsers
.. autosummary::
    str_to_rdmol
    FormatEnum

API
---
.. autofunction:: str_to_rdmol
.. autoclass:: FormatEnum

"""

# flake8: noqa: E501

import enum
import textwrap
from typing import Any, Callable

from nanoutils import PartialPrepend
from packaging.version import Version

try:
    import rdkit
    from rdkit import Chem
except ModuleNotFoundError:
    # Precaution against sphinx mocking raising a TypeError
    FLAGS = 0
    RDKIT_VERSION = Version("0.0.0")
else:
    FLAGS = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    RDKIT_VERSION = Version(rdkit.__version__)

__all__ = ["FormatEnum", "str_to_rdmol"]

SMILES_PARAMS = Chem.SmilesParserParams()
SMILES_PARAMS.removeHs = False


def str_to_rdmol(
    string: str,
    parser: Callable[[str], Chem.Mol],
    kind: str = "string",
    **kwargs: Any,
) -> Chem.Mol:
    """Convert a SMILES string into an rdkit Mol; supports explicit hydrogens."""
    try:
        mol = parser(string, **kwargs)
        assert mol is not None
        if not kwargs.get("sanitize", True):
            Chem.SanitizeMol(mol, sanitizeOps=FLAGS)
    except Exception as ex:
        err_str = string if "\n" not in string else "\n" + textwrap.indent(string, 4 * " ")
        raise ValueError(f'Failed to parse the following {kind}: {err_str!r}') from ex
    return mol


class FormatEnum(enum.Enum):
    """An enum with various rdkit-based string-parsing options."""

    FASTA = PartialPrepend(str_to_rdmol, Chem.MolFromFASTA, "FASTA string", sanitize=False)
    HELM = PartialPrepend(str_to_rdmol, Chem.MolFromHELM, "HELM string", sanitize=False)
    INCHI = PartialPrepend(str_to_rdmol, Chem.MolFromInchi, "InChi string", sanitize=False, removeHs=False)
    MOL2 = PartialPrepend(str_to_rdmol, Chem.MolFromMol2Block, "Mol2 block", sanitize=False, removeHs=False)
    MOL2_FILE = PartialPrepend(str_to_rdmol, Chem.MolFromMol2File, "Mol2 file", sanitize=False, removeHs=False)
    MOL = PartialPrepend(str_to_rdmol, Chem.MolFromMolBlock, "Mol block", sanitize=False, removeHs=False)
    MOL_FILE = PartialPrepend(str_to_rdmol, Chem.MolFromMolFile, "Mol file", sanitize=False, removeHs=False)
    PDB = PartialPrepend(str_to_rdmol, Chem.MolFromPDBBlock, "PDB block", sanitize=False, removeHs=False)
    PDB_FILE = PartialPrepend(str_to_rdmol, Chem.MolFromPDBFile, "PDB file", sanitize=False, removeHs=False)
    SEQUENCE = PartialPrepend(str_to_rdmol, Chem.MolFromSequence, "sequence string", sanitize=False)
    SMARTS = PartialPrepend(str_to_rdmol, Chem.MolFromSmarts, "SMARTS string")
    SMILES = PartialPrepend(str_to_rdmol, Chem.MolFromSmiles, "SMILES string", params=SMILES_PARAMS)
    TPL = PartialPrepend(str_to_rdmol, Chem.MolFromTPLBlock, "TPL block", sanitize=False)
    TPL_FILE = PartialPrepend(str_to_rdmol, Chem.MolFromTPLFile, "TPL file", sanitize=False)
    if RDKIT_VERSION >= Version("2020.03.3"):
        SVG = PartialPrepend(str_to_rdmol, Chem.MolFromRDKitSVG, "SVG string", sanitize=False, removeHs=False)
    if RDKIT_VERSION >= Version("2020.09.1"):
        PNG = PartialPrepend(str_to_rdmol, Chem.MolFromPNGString, "PNG string")
        PNG_FILE = PartialPrepend(str_to_rdmol, Chem.MolFromPNGFile, "PNG file")
