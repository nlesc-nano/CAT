"""Tests for :mod:`CAT.attachment.ligand_anchoring`."""

# flake8: noqa: E501

import os
import sys
import math
import h5py
from shutil import rmtree
from os.path import join
from typing import Tuple, Generator, Dict, Any
from pathlib import Path

import rdkit
import pytest
import numpy as np
from unittest import mock
from rdkit import Chem
from scm.plams import from_smiles, Molecule, to_rdmol, PT
from assertionlib import assertion
from schema import SchemaError
from packaging.version import Version

from CAT.utils import get_template, KindEnum, AnchorTup, FormatEnum, MultiAnchorEnum
from CAT.base import prep_input
from CAT.attachment.ligand_anchoring import (
    get_functional_groups, _smiles_to_rdmol, find_substructure, init_ligand_anchoring
)
from CAT.attachment.ligand_opt import optimize_ligand
from CAT.data_handling.anchor_parsing import parse_anchors, INVALID_SMILES_ATOMS

if sys.version_info >= (3, 7):
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict

PATH = Path('tests') / 'test_files'


class TestGetFunctionalGroups:
    """Tests for :func:`CAT.attachment.ligand_anchoring.get_functional_groups`."""

    @staticmethod
    def get_symbols_and_bonds(mol: Chem.Mol) -> Tuple[np.ndarray, np.ndarray]:
        atoms_iter = ((at.GetFormalCharge(), at.GetSymbol()) for at in mol.GetAtoms())
        atoms_dtype = [("charge", "i8"), ("symbol", "S2")]
        atoms_ar = np.fromiter(atoms_iter, dtype=atoms_dtype)

        bonds_iter = ((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds())
        bonds_dtype = [("atom1", "i8"), ("atom2", "i8")]
        bonds_ar = np.fromiter(bonds_iter, dtype=bonds_dtype)
        return atoms_ar, bonds_ar

    @pytest.fixture(scope="class", autouse=True, name="group")
    def get_h5py_file(self) -> Generator[h5py.Group, None, None]:
        with h5py.File(PATH / "test_ligand_anchoring.hdf5", "r") as f:
            yield f["TestGetFunctionalGroups"]

    def test_func_groups1(self, group: h5py.Group) -> None:
        smiles_list = ['[O-]C', '[O-]CC', '[O-]CCC']
        func_groups = dict(zip(smiles_list, get_functional_groups(smiles_list)))
        for smiles, mol in func_groups.items():
            atoms, bonds = self.get_symbols_and_bonds(mol)
            atoms_ref = group[f'test_func_groups1/{smiles}/atoms'][...]
            bonds_ref = group[f'test_func_groups1/{smiles}/bonds'][...]

            np.testing.assert_array_equal(atoms, atoms_ref)
            np.testing.assert_array_equal(bonds, bonds_ref)

    def test_func_groups2(self, group: h5py.Group) -> None:
        smiles_list = get_template('smiles.yaml').split
        _func_groups = get_functional_groups(split=True)
        assertion.eq(len(_func_groups), len(smiles_list))

        func_groups = dict(zip(smiles_list, _func_groups))
        for smiles, mol in func_groups.items():
            atoms, bonds = self.get_symbols_and_bonds(mol)
            atoms_ref = group[f'test_func_groups2/{smiles}/atoms'][...]
            bonds_ref = group[f'test_func_groups2/{smiles}/bonds'][...]

            np.testing.assert_array_equal(atoms, atoms_ref)
            np.testing.assert_array_equal(bonds, bonds_ref)

    def test_func_groups3(self, group: h5py.Group) -> None:
        smiles_list = get_template('smiles.yaml').no_split
        _func_groups = get_functional_groups(split=False)
        assertion.eq(len(_func_groups), len(smiles_list))

        func_groups = dict(zip(smiles_list, _func_groups))
        for smiles, mol in func_groups.items():
            atoms, bonds = self.get_symbols_and_bonds(mol)
            atoms_ref = group[f'test_func_groups3/{smiles}/atoms'][...]
            bonds_ref = group[f'test_func_groups3/{smiles}/bonds'][...]

            np.testing.assert_array_equal(atoms, atoms_ref)
            np.testing.assert_array_equal(bonds, bonds_ref)


def test_smiles_to_rdmol() -> None:
    """Tests for :meth:`CAT.attachment.ligand_anchoring._smiles_to_rdmol`."""
    smiles1 = 'CO'
    rdmol1 = _smiles_to_rdmol(smiles1)
    C, O_ = rdmol1.GetAtoms()

    assertion.len_eq(C.GetBonds(), 1)
    assertion.eq(C.GetSymbol(), 'C')
    assertion.eq(C.GetFormalCharge(), 0)
    assertion.len_eq(O_.GetBonds(), 1)
    assertion.eq(O_.GetSymbol(), 'O')
    assertion.eq(O_.GetFormalCharge(), 0)

    smiles2 = 'CO[H]'
    rdmol2 = _smiles_to_rdmol(smiles2)
    C, O_, H = rdmol2.GetAtoms()

    assertion.len_eq(C.GetBonds(), 1)
    assertion.eq(C.GetSymbol(), 'C')
    assertion.eq(C.GetFormalCharge(), 0)
    assertion.len_eq(O_.GetBonds(), 2)
    assertion.eq(O_.GetSymbol(), 'O')
    assertion.eq(O_.GetFormalCharge(), 0)
    assertion.len_eq(H.GetBonds(), 1)
    assertion.eq(H.GetSymbol(), 'H')
    assertion.eq(H.GetFormalCharge(), 0)


class TestFindSubstructure:
    """Tests for :func:`CAT.attachment.ligand_anchoring.find_substructure`."""

    @staticmethod
    def get_symbols_and_bonds(mol: Molecule) -> Tuple[np.ndarray, np.ndarray]:
        atoms_iter = ((at.properties.get("charge", 0.0), at.symbol) for at in mol)
        atoms_dtype = [("charge", "f8"), ("symbol", "S2")]
        atoms_ar = np.fromiter(atoms_iter, dtype=atoms_dtype)

        try:
            mol.set_atoms_id()
            bonds_iter = ((b.atom1.id, b.atom2.id, b.order) for b in mol.bonds)
            bonds_dtype = [("atom1", "i8"), ("atom2", "i8"), ("order", "i8")]
            bonds_ar = np.fromiter(bonds_iter, dtype=bonds_dtype)
        finally:
            mol.unset_atoms_id()
        return atoms_ar, bonds_ar

    @pytest.fixture(scope="class", autouse=True, name="group")
    def get_h5py_group(self) -> Generator[h5py.Group, None, None]:
        with h5py.File(PATH / "test_ligand_anchoring.hdf5", "r") as f:
            yield f["TestFindSubstructure"]

    OPTIONS_DICT = OrderedDict({
        "remove": dict(group="OC(=O)C", group_idx=1, remove=[0, 2, 3]),
        "kind_first": dict(group="OC(=O)C", group_idx=0, kind="FIRST"),
        "kind_mean": dict(group="OC(=O)C", group_idx=[0, 2], kind="MEAN"),
        "kind_mean_translate": dict(group="OC(=O)C", group_idx=[0, 2], kind="MEAN_TRANSLATE"),
        "angle": dict(group="OC(=O)C", group_idx=[0, 1, 2], angle_offset=45),
    })

    @pytest.mark.parametrize("kwargs_id,kwargs", OPTIONS_DICT.items(), ids=OPTIONS_DICT.keys())
    @pytest.mark.xfail(
        Version(rdkit.__version__) < Version("2021.03.4"),
        reason="requires rdkit >= 2021.03.4",
    )
    def test_options(self, kwargs_id: str, kwargs: Dict[str, Any], group: h5py.Group) -> None:
        # Parse the anchoring options
        mol = from_smiles("OC(=O)CCCCC")
        anchor_tup = parse_anchors(kwargs)

        # Apply the anchoring options and optimize
        mol_list = find_substructure(mol, anchor_tup)
        assertion.len_eq(mol_list, 1)
        mol_out = mol_list[0]
        optimize_ligand(mol_out)

        # Analyze the results
        atoms, bonds = self.get_symbols_and_bonds(mol_out)
        coords = np.array(mol)
        atoms_ref = group[f'test_options/{kwargs_id}/atoms'][...]
        coords_ref = group[f'test_options/{kwargs_id}/coords'][...]
        bonds_ref = group[f'test_options/{kwargs_id}/bonds'][...]

        np.testing.assert_array_equal(atoms, atoms_ref)
        np.testing.assert_array_equal(bonds, bonds_ref)
        np.testing.assert_allclose(coords, coords_ref, atol=10e-3)

    def test_split(self) -> None:
        func_groups = parse_anchors(split=True)
        mol = from_smiles('O=C(O)C(O)N')
        out = find_substructure(mol, func_groups)

        assertion.len_eq(out, 3)
        mol1, mol2, mol3 = out
        for m in out:
            assertion.len_eq(m, len(mol) - 1)

        assertion.eq(mol1[1].properties.charge, -1)
        assertion.eq(mol1[1].symbol, 'N')
        assertion.len_eq(mol1[1].bonds, 2)
        assertion.is_(mol1.properties.dummies, mol1[1])
        assertion.eq(mol1.properties.anchor, 'N1')
        assertion.eq(mol1.properties.smiles, '[NH-]C(O)C(=O)O')
        assertion.is_not(mol1, mol)

        assertion.eq(mol2[3].properties.charge, -1)
        assertion.eq(mol2[3].symbol, 'O')
        assertion.len_eq(mol2[3].bonds, 1)
        assertion.is_(mol2.properties.dummies, mol2[3])
        assertion.eq(mol2.properties.anchor, 'O3')
        assertion.eq(mol2.properties.smiles, 'NC([O-])C(=O)O')
        assertion.is_not(mol1, mol)

        assertion.eq(mol3[6].properties.charge, -1)
        assertion.eq(mol3[6].symbol, 'O')
        assertion.len_eq(mol3[6].bonds, 1)
        assertion.is_(mol3.properties.dummies, mol3[6])
        assertion.eq(mol3.properties.anchor, 'O6')
        assertion.eq(mol3.properties.smiles, 'NC(O)C(=O)[O-]')
        assertion.is_not(mol1, mol)

        mol_invalid = from_smiles('CCC')
        out_invalid = find_substructure(mol_invalid, func_groups)
        assertion.eq(out_invalid, [])

    def test_no_split(self) -> None:
        mol = from_smiles('O=C(O)C(O)N')
        func_groups_nosplit = parse_anchors(split=False)
        out_nosplit = find_substructure(mol, func_groups_nosplit, split=False)
        for m in out_nosplit:
            assertion.len_eq(m, len(mol))
            assertion.is_not(m, mol)


@mock.patch.dict(
    os.environ,
    {'AMSBIN': '', 'AMSHOME': '', 'AMSRESOURCES': '', 'SCMLICENSE': ''},
)
def test_init_ligand_anchoring() -> None:
    """Tests for :meth:`CAT.attachment.ligand_anchoring.init_ligand_anchoring`."""
    try:
        filename = join(PATH, 'input1.yaml')
        s = get_template(filename, from_cat_data=False)
        ligand_df, *_ = prep_input(s)
        df = init_ligand_anchoring(ligand_df)

        idx = [('C[O-]', 'O2'), ('CC[O-]', 'O3')]
        assertion.eq(df.index.tolist(), idx)

        columns = [('mol', ''), ('formula', ''), ('hdf5 index', ''), ('opt', '')]
        assertion.eq(df.columns.tolist(), columns)

        formula_ref = np.array(['C1H3O1', 'C2H5O1'], dtype=object)
        formula = df[('formula', '')].values
        np.testing.assert_array_equal(formula, formula_ref)

        opt_ref = np.array([False, False])
        opt = df[('opt', '')].values
        np.testing.assert_array_equal(opt, opt_ref)

        assertion.len_eq(df[('mol', '')].values[0], 5)
        assertion.len_eq(df[('mol', '')].values[1], 8)
    finally:
        rmtree(join(PATH, 'ligand'))
        rmtree(join(PATH, 'qd'))
        rmtree(join(PATH, 'database'))


class TestInputParsing:
    PARAM_RAISE: "OrderedDict[str, tuple[Any, type[Exception]]]" = OrderedDict(
        invalid_smiles=("test", ValueError),
        idx_len_0=({"group": "OC", "group_idx": []}, SchemaError),
        idx_len_2=({"group": "OC", "group_idx": 0, "angle_offset": 45}, ValueError),
        idx_len_3=({"group": "OC", "group_idx": 0, "dihedral": 0.5}, ValueError),
        duplicate_idx=({"group": "OC", "group_idx": [0, 0]}, SchemaError),
        invalid_idx_type=({"group": "OC", "group_idx": 0.0}, SchemaError),
        idx_intersection=({"group": "OC", "group_idx": 0, "remove": 0}, ValueError),
        out_of_bounds_idx=({"group": "OC", "group_idx": 99}, IndexError),
        out_of_bounds_remove=({"group": "OC", "group_idx": 0, "remove": 99}, IndexError),
        angle_unit=({"group": "OC", "group_idx": 0, "angle_offset": "0.5 bob"}, SchemaError),
        angle_invalid=({"group": "OC", "group_idx": 0, "angle_offset": "bob"}, SchemaError),
        invalid_group=({"group": "OC", "group_idx": 0, "group": 1.0}, SchemaError),
        invalid_group_format=({"group": "OC", "group_idx": 0, "group_format": 1}, SchemaError),
        invalid_kind=({"group": "OC", "group_idx": 0, "kind": 1}, SchemaError),
        invalid_multi_anchor_filter=({"group": "OC", "group_idx": 0, "multi_anchor_filter": 1}, SchemaError),
    )

    @pytest.mark.parametrize("inp,exc_type", PARAM_RAISE.values(), ids=PARAM_RAISE.keys())
    def test_raise(self, inp: Any, exc_type: "type[Exception]") -> None:
        with pytest.raises(exc_type):
            parse_anchors(inp)

    PARAM_RAISE_CORE: "OrderedDict[str, tuple[Any, type[Exception]]]" = OrderedDict(
        none=(None, TypeError),
        angle_offset=({"group": "Cl", "group_idx": 0, "angle_offset": 90}, TypeError),
        dihedral=({"group": "Cl", "group_idx": 0, "dihedral": 90}, TypeError),
        multiple=(["OC", "OCC"], NotImplementedError),
        kind=({"group": "Cl", "group_idx": 0, "kind": "mean_translate"}, ValueError),
    )

    @pytest.mark.parametrize("inp,exc_type", PARAM_RAISE_CORE.values(), ids=PARAM_RAISE_CORE)
    def test_raise_core(self, inp: Any, exc_type: "type[Exception]") -> None:
        with pytest.raises(exc_type):
            parse_anchors(inp, is_core=True)

    _PARAM_PASS1 = OrderedDict(
        idx_scalar={"group": "OCC", "group_idx": 0},
        idx_list={"group": "OCC", "group_idx": [0]},
        list=[{"group": "OCC", "group_idx": 0}],
        str_COH=["O(C)[H]"],
        str_Cd=["Cd"],
        str_Cl=["Cl"],
        str_dummy=["Xx"],
        int_Cd=[48],
        int_Cl=[17],
        int_dummy=[0],
        group_cd={"group": "Cd", "group_idx": 0},
        group_cl={"group": "Cl", "group_idx": 0},
        angle_unit={"group": "OCC", "group_idx": range(3), "angle_offset": "1 rad"},
        angle_no_unit={"group": "OCC", "group_idx": range(3), "angle_offset": "180"},
        angle_none={"group": "OCC", "group_idx": range(3), "angle_offset": None},
        angle_float={"group": "OCC", "group_idx": range(3), "angle_offset": 180.0},
        remove_none={"group": "OCC", "group_idx": 0, "remove": None},
        kind_none={"group": "OCC", "group_idx": 0, "kind": None},
        kind_str={"group": "OCC", "group_idx": 0, "kind": "mean"},
        kind_enum={"group": "OCC", "group_idx": 0, "kind": KindEnum.MEAN_TRANSLATE},
        group_format_none={"group": "OCC", "group_idx": 0, "group_format": None},
        group_format_str={"group": "OCC", "group_idx": 0, "group_format": "SMARTS"},
        group_format_enum={"group": "OCC", "group_idx": 0, "group_format": FormatEnum.SMARTS},
        multi_anchor_filter_none={"group": "OCC", "group_idx": 0, "multi_anchor_filter": None},
        multi_anchor_filter_str={"group": "OCC", "group_idx": 0, "multi_anchor_filter": "ALL"},
        multi_anchor_filter_enum={"group": "OCC", "group_idx": 0, "multi_anchor_filter": MultiAnchorEnum.ALL},
    )
    _PARAM_PASS2 = OrderedDict(
        idx_scalar=AnchorTup(None, group="OCC", group_idx=(0,)),
        idx_list=AnchorTup(None, group="OCC", group_idx=(0,)),
        list=AnchorTup(None, group="OCC", group_idx=(0,)),
        str_COH=AnchorTup(None, group="O(C)[H]", group_idx=(0,), remove=(2,)),
        str_Cd=AnchorTup(None, group="Cd", group_idx=(0,), remove=(0,)),
        str_Cl=AnchorTup(None, group="Cl", group_idx=(0,), remove=(0,)),
        str_dummy=AnchorTup(None, group="*", group_idx=(0,), remove=(0,)),
        int_Cd=AnchorTup(None, group="Cd", group_idx=(0,), remove=(0,)),
        int_Cl=AnchorTup(None, group="Cl", group_idx=(0,), remove=(0,)),
        int_dummy=AnchorTup(None, group="*", group_idx=(0,), remove=(0,)),
        group_cd=AnchorTup(None, group="Cd", group_idx=(0,)),
        group_cl=AnchorTup(None, group="Cl", group_idx=(0,)),
        angle_unit=AnchorTup(None, group="OCC", group_idx=(0, 1, 2), angle_offset=1.0),
        angle_no_unit=AnchorTup(None, group="OCC", group_idx=(0, 1, 2), angle_offset=math.pi),
        angle_none=AnchorTup(None, group="OCC", group_idx=(0, 1, 2), angle_offset=None),
        angle_float=AnchorTup(None, group="OCC", group_idx=(0, 1, 2), angle_offset=math.pi),
        remove_none=AnchorTup(None, group="OCC", group_idx=(0,), remove=None),
        kind_none=AnchorTup(None, group="OCC", group_idx=(0,), kind=KindEnum.FIRST),
        kind_str=AnchorTup(None, group="OCC", group_idx=(0,), kind=KindEnum.MEAN),
        kind_enum=AnchorTup(None, group="OCC", group_idx=(0,), kind=KindEnum.MEAN_TRANSLATE),
        group_format_none=AnchorTup(None, group="OCC", group_idx=(0,), group_format=FormatEnum.SMILES),
        group_format_str=AnchorTup(None, group="OCC", group_idx=(0,), group_format=FormatEnum.SMARTS),
        group_format_enum=AnchorTup(None, group="OCC", group_idx=(0,), group_format=FormatEnum.SMARTS),
        multi_anchor_filter_none=AnchorTup(None, group="OCC", group_idx=(0,), multi_anchor_filter=MultiAnchorEnum.ALL),
        multi_anchor_filter_str=AnchorTup(None, group="OCC", group_idx=(0,), multi_anchor_filter=MultiAnchorEnum.ALL),
        multi_anchor_filter_enum=AnchorTup(None, group="OCC", group_idx=(0,), multi_anchor_filter=MultiAnchorEnum.ALL),
    )
    PARAM_PASS = OrderedDict({
        k: (v1, v2) for (k, v1), v2 in zip(_PARAM_PASS1.items(), _PARAM_PASS2.values())
    })

    @pytest.mark.parametrize("inp,ref", PARAM_PASS.values(), ids=PARAM_PASS.keys())
    def test_pass(self, inp: Any, ref: AnchorTup) -> None:
        out_tup = parse_anchors(inp)
        assertion.len_eq(out_tup, 1)
        out = out_tup[0]

        assertion.isinstance(out.mol, Chem.Mol)
        if out.angle_offset is not None:
            assertion.isclose(out.angle_offset, ref.angle_offset)
        assertion.eq(out._replace(mol=None, angle_offset=None), ref._replace(angle_offset=None))

    @pytest.mark.parametrize("split", [True, False], ids=["split", "no_split"])
    def test_rdkit_mol(self, split: bool) -> None:
        remove = (1,) if split else None
        mol = _smiles_to_rdmol("[O-]C")
        out = parse_anchors(mol, split=split)

        ref = AnchorTup(mol, group=None, group_idx=(0,), remove=remove)
        assertion.len_eq(out, 1)
        assertion.eq(out[0], ref)

    def test_anchor_tup(self) -> None:
        ref = AnchorTup(_smiles_to_rdmol("[O-]C"), group="[O-]C", group_idx=(0,))
        out = parse_anchors(ref)
        assertion.len_eq(out, 1)
        assertion.eq(out[0], ref)

    SPLIT_REF = [
        "C[N+].[F-]",
        "C[N+].[Cl-]",
        "C[N+].[Br-]",
        "C[N+].[I-]",
        "[H]NC",
        "[H]PC",
        "[H]OP",
        "[H]OC",
        "[H]SC",
        "[H]OS",
    ]
    NO_SPLIT_REF = [
        "C[N+]",
        "CN",
        "C[N-]",
        "CP",
        "C[P-]",
        "OP",
        "[O-]P",
        "CO",
        "C[O-]",
        "CS",
        "C[S-]",
        "OS",
        "[O-]S",
    ]

    @pytest.mark.parametrize("split,ref", [
        (True, SPLIT_REF),
        (False, NO_SPLIT_REF),
    ], ids=["split", "no_split"])
    def test_none(self, split: bool, ref: "list[str]") -> None:
        out = parse_anchors(split=split)
        smiles = [Chem.MolToSmiles(tup.mol) for tup in out]
        assertion.eq(smiles, ref)

    def test_atomic_symbol(self) -> None:
        """Test that groups consisting of a single atomic symbol are properly processed."""
        out_tup = parse_anchors({"group": "Cd", "group_idx": 0})
        assertion.len_eq(out_tup, 1)
        out = out_tup[0]

        mol = Molecule(PATH / "core" / "Cd68Se55.xyz")
        rdmol = to_rdmol(mol)
        match = np.ravel(rdmol.GetSubstructMatches(out.mol))
        match.sort()

        ref = np.fromiter((i for i, at in enumerate(mol) if at.symbol == "Cd"), dtype=np.int64)
        np.testing.assert_array_equal(match, ref)

    def test_invalid_smiles_atoms(self) -> None:
        """Check that only atom types that lack RDKit SMILES support are \
        included in ``INVALID_SMILES_ATOMS``."""
        smiles_parser = FormatEnum.SMILES.value
        symbol_set = set()
        for at in PT.symtonum:
            try:
                mol = smiles_parser(at)
            except ValueError:
                symbol_set.add(at)

        diff = symbol_set ^ INVALID_SMILES_ATOMS
        diff -= {"Xx"}
        assertion.not_(diff)
