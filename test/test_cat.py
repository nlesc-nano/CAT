"""A module for CAT-related tests."""

import contextlib
import io
import os
import shutil

import yaml
import pandas as pd

from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule

from CAT import Database
from CAT.base import (prep_input, prep_core, prep_ligand)


# prepare input
ARG = Settings(yaml.load(
    """
    path: /Users/basvanbeek/Documents/GitHub/CAT/test

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC
        - OCC
        - OCCC
        - OCCCC

    optional:
        database:
            dirname: database
            read: True
            write: True
            overwrite: False
            mol_format: [xyz, pdb]
            mongodb: False

        core:
            dirname: core
            dummy: Cl

        ligand:
            dirname: ligand
            functional_groups: None
            optimize: True
            split: True
            cosmo-rs: False

        qd:
            dirname: QD
            optimize: False
            activation_strain: False
            dissociate: False
    """,
    Loader=yaml.FullLoader))
LIGAND_DF, CORE_DF = prep_input(ARG)
shutil.rmtree(ARG.optional.database.dirname)


def test_prep_core():
    """ Test the :func:`CAT.base.prep_core` function. """
    arg = ARG.copy()
    core_df = CORE_DF.copy()

    # Check the dataframe
    ret = prep_core(core_df, arg)
    assert isinstance(ret, pd.DataFrame)
    assert ret.shape == core_df.shape
    assert 'mol' in ret.columns

    # Check the molecule in the dataframe
    core = ret['mol'][0]
    assert isinstance(core, Molecule)
    assert len(core) == 123
    assert 'Cl' not in [at.symbol for at in core]
    assert len(core.properties.dummies) == 26
    assert len(set([at.symbol for at in core.properties.dummies])) == 1


def test_prep_ligand():
    """ Test the :func:`CAT.base.prep_ligand` function with **split** = *True*. """
    arg = ARG.copy()
    lig_df = LIGAND_DF.copy()
    if os.path.isdir(arg.optional.database.dirname):
        shutil.rmtree(arg.optional.database.dirname)
    os.mkdir(arg.optional.database.dirname)
    data = Database(path=arg.optional.database.dirname)

    # Check while splite=False
    arg.optional.ligand.split = False
    ret = prep_ligand(lig_df, arg)
    assert isinstance(ret, pd.DataFrame)
    assert 'mol' in ret.columns
    assert ret['mol'].shape[0] == lig_df.shape[0]

    # Check the molecules in the dataframe
    lig_list = ret['mol'].values.tolist()
    with data.open_csv_lig(data.csv_lig) as db:
        assert [6, 9, 12, 15] == sorted([len(lig) for lig in lig_list])
        for lig in lig_list:
            assert isinstance(lig, Molecule)
            assert 'O' in lig.properties.anchor
            assert lig.properties.charge == 0
            assert lig.properties.dummies.properties.charge == 0
            assert lig.properties.name + '.pdb' in os.listdir(arg.optional.ligand.dirname)
            assert lig.properties.name + '.xyz' in os.listdir(arg.optional.ligand.dirname)
            assert (lig.properties.smiles, lig.properties.anchor) in db.index

    # Check if previous structures can be pulled from the database
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        prep_ligand(lig_df, arg)
    print_list = f.getvalue().splitlines()
    for item in print_list:
        if item:
            assert 'has been pulled from the database' in item

    # Reset the directories
    shutil.rmtree(arg.optional.database.dirname)
    shutil.rmtree(arg.optional.ligand.dirname)
    os.mkdir(arg.optional.ligand.dirname)


def test_prep_ligand_split():
    """ Test the :func:`CAT.base.prep_ligand` function with **split** = *False*. """
    arg = ARG.copy()
    lig_df = LIGAND_DF.copy()
    if os.path.isdir(arg.optional.database.dirname):
        shutil.rmtree(arg.optional.database.dirname)
    os.mkdir(arg.optional.database.dirname)
    data = Database(path=arg.optional.database.dirname)

    # Check the dataframe
    ret = prep_ligand(lig_df, arg)
    assert isinstance(ret, pd.DataFrame)
    assert 'mol' in ret.columns
    assert ret['mol'].shape[0] == lig_df.shape[0]

    # Check the molecules in the dataframe
    lig_list = ret['mol'].values.tolist()
    with data.open_csv_lig(data.csv_lig) as db:
        assert [5, 8, 11, 14] == sorted([len(lig) for lig in lig_list])
        for lig in lig_list:
            assert isinstance(lig, Molecule)
            assert 'O' in lig.properties.anchor
            assert lig.properties.charge == -1
            assert lig.properties.dummies.properties.charge == -1
            assert lig.properties.name + '.pdb' in os.listdir(arg.optional.ligand.dirname)
            assert lig.properties.name + '.xyz' in os.listdir(arg.optional.ligand.dirname)
            assert (lig.properties.smiles, lig.properties.anchor) in db.index

    # Check if previous structures can be pulled from the database
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        prep_ligand(lig_df, arg)
    print_list = f.getvalue().splitlines()
    for item in print_list:
        if item:
            assert 'has been pulled from the database' in item

    # Reset the directories
    shutil.rmtree(arg.optional.database.dirname)
    shutil.rmtree(arg.optional.ligand.dirname)
    os.mkdir(arg.optional.ligand.dirname)
