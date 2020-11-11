"""Tests for :mod:`CAT.data_handling.validate_input`."""

import os
from os.path import join
from pathlib import Path
from itertools import chain

import yaml
from unittest import mock

from rdkit import Chem
from scm.plams import (Settings, AMSJob)
from assertionlib import assertion
from nanoutils import delete_finally

from CAT.data_handling.validate_input import validate_input
from dataCAT import Database

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_validate_input() -> None:
    """Test :func:`CAT.data_handling.validate_input.validate_input`."""
    with open(join(PATH, 'input1.yaml'), 'r') as f:
        s = Settings(yaml.load(f, Loader=yaml.FullLoader))
    s.path = PATH
    validate_input(s, validate_only=False)

    s2 = s.copy()
    validate_input(s2, validate_only=False)

    ref = Settings()
    ref.core.dirname = join(PATH, 'core')
    ref.core.anchor = 35
    ref.core.allignment = 'surface'
    ref.core.subset = None

    ref.database.dirname = join(PATH, 'database')
    ref.database.mol_format = ('pdb',)
    ref.database.mongodb = {}
    ref.database.overwrite = ()
    ref.database.read = ('core', 'ligand', 'qd')
    ref.database.write = ('core', 'ligand', 'qd')
    ref.database.thread_safe = False
    ref.database.db = Database(ref.database.dirname, **ref.database.mongodb)

    ref.ligand['cosmo-rs'] = False
    ref.ligand.dirname = join(PATH, 'ligand')
    ref.ligand.optimize = {'job1': None, 'job2': None, 's1': None, 's2': Settings(),
                           'use_ff': False, 'keep_files': True}
    ref.ligand.split = True
    ref.ligand.cdft = False

    ref.qd.bulkiness = False
    ref.qd.construct_qd = True
    ref.qd.activation_strain = False
    ref.qd.dirname = join(PATH, 'qd')
    ref.qd.dissociate = False
    ref.qd.multi_ligand = None
    ref.qd.optimize = {'job1': AMSJob, 'keep_files': True, 'use_ff': False, 's2': {'description': 'UFF with the default forcefield', 'input': {'uff': {'library': 'uff'}, 'ams': {'system': {'bondorders': {'_1': None}}}}}, 's1': {'description': 'UFF with the default forcefield', 'input': {'uff': {'library': 'uff'}, 'ams': {'system': {'bondorders': {'_1': None}}}}}, 'job2': AMSJob}  # noqa

    ref.forcefield = Settings()

    func_groups1 = s.optional.ligand.pop('anchor')
    func_groups2 = s2.optional.ligand.pop('anchor')

    for mol in chain(func_groups1, func_groups2):
        assertion.isinstance(mol, Chem.Mol)
    assertion.eq(s.optional, ref)
    assertion.eq(s2.optional, ref)
