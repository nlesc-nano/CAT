import sys
import textwrap

import pytest
from assertionlib import assertion
from rdkit import Chem

from CAT.utils import FormatEnum

if sys.version_info >= (3, 7):
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict

PARAM = OrderedDict(
    FASTA="MDSKGSSQ",
    HELM="RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$",
    INCHI=r"InChI=1S/H2O/h1H2",
    MOL=textwrap.dedent("""
             RDKit          2D

          1  0  0  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
    """),
    PDB=textwrap.dedent("""\
        HETATM    1  O1  UNL     1       0.000   0.000   0.000  1.00  0.00           O
        END
    """),
    SMARTS=r"[N&H2&+0&D1:6]/[N&H0&+0&D2:5]=[c&H0&+0&D3:4](/[n:3]):[n&H1&+0&D2:2]:[c:1]",
    SMILES="CCCCCO",
)


@pytest.mark.parametrize("kind,string", PARAM.items(), ids=PARAM)
def test_str_parser(kind: str, string: str) -> None:
    func = FormatEnum[kind].value
    mol = func(string)
    assertion.isinstance(mol, Chem.Mol)


@pytest.mark.parametrize("kind", FormatEnum.__members__)
def test_raise(kind: str) -> None:
    func = FormatEnum[kind].value
    with pytest.raises(ValueError):
        func(r"bo%b")
