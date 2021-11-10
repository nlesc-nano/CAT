"""Tests for :mod:`CAT.attachment.ligand_anchoring`."""

import os
import sys
import h5py
from shutil import rmtree
from os.path import join
from typing import Tuple, Generator, Dict, Any
from pathlib import Path

import pytest
import numpy as np
from unittest import mock
from rdkit import Chem
from scm.plams import from_smiles, Molecule
from assertionlib import assertion

from CAT.utils import get_template
from CAT.base import prep_input
from CAT.attachment.ligand_anchoring import (
    get_functional_groups, _smiles_to_rdmol, find_substructure, init_ligand_anchoring
)
from CAT.attachment.ligand_opt import optimize_ligand
from CAT.data_handling.anchor_parsing import parse_anchors

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
    def get_h5py_file(self) -> Generator[h5py.Group, None, None]:
        with h5py.File(PATH / "test_ligand_anchoring.hdf5", "r") as f:
            yield f["TestFindSubstructure"]

    OPTIONS_DICT = OrderedDict({
        "remove": dict(group="OC(=O)C", group_idx=1, remove=[0, 2, 3]),
        "kind_first": dict(group="OC(=O)C", group_idx=0, kind="FIRST"),
        "kind_mean": dict(group="OC(=O)C", group_idx=[0, 2], kind="MEAN"),
        "kind_mean_translate": dict(group="OC(=O)C", group_idx=[0, 2], kind="MEAN_TRANSLATE"),
    })

    @pytest.mark.parametrize("kwargs_id,kwargs", OPTIONS_DICT.items(), ids=OPTIONS_DICT.keys())
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
        try:
            np.testing.assert_allclose(coords, coords_ref, atol=10e-3)
        except AssertionError:
            if sys.version_info >= (3, 9):
                pytest.xfail("Geometries must be updated for RDKit >2019.09.2")
            raise

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
