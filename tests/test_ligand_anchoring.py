"""Tests for :mod:`CAT.attachment.ligand_anchoring`."""

import os
from shutil import rmtree
from os.path import join

import numpy as np
from unittest import mock

from rdkit import Chem
from scm.plams.interfaces.molecule.rdkit import from_smiles
from assertionlib import assertion

from CAT.utils import get_template
from CAT.base import prep_input
from CAT.attachment.ligand_anchoring import (
    get_functional_groups, _smiles_to_rdmol, find_substructure, init_ligand_anchoring
)

PATH = join('tests', 'test_files')


def test_get_functional_groups() -> None:
    """Tests for :meth:`CAT.attachment.ligand_anchoring.get_functional_groups`."""
    _func_groups = ['[O-]C', '[O-]CC', '[O-]CCC']
    at_count = (2, 3, 4)
    func_groups1 = get_functional_groups(_func_groups)
    for mol, v in zip(func_groups1, at_count):
        atoms = list(mol.GetAtoms())
        assertion.len_eq(atoms, v)
        assertion.len_eq(atoms[0].GetBonds(), 1)
        assertion.eq(atoms[0].GetSymbol(), 'O')
        assertion.eq(atoms[0].GetFormalCharge(), -1)

        for at in atoms[1:]:
            assertion.eq(at.GetSymbol(), 'C')
            assertion.eq(at.GetFormalCharge(), 0)

            if at is atoms[-1]:  # Terminal atom
                assertion.len_eq(at.GetBonds(), 1)
            else:
                assertion.len_eq(at.GetBonds(), 2)

    func_groups2 = get_functional_groups(split=True)
    assertion.len_eq(func_groups2, 10)
    for m in func_groups2:
        assertion.isinstance(m, Chem.Mol)

    func_groups3 = get_functional_groups(split=False)
    assertion.len_eq(func_groups3, 13)
    for m in func_groups3:
        assertion.isinstance(m, Chem.Mol)


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


def test_find_substructure() -> None:
    """Tests for :meth:`CAT.attachment.ligand_anchoring.find_substructure`."""
    func_groups = get_functional_groups(split=True)
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

    func_groups_nosplit = get_functional_groups(split=False)
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
