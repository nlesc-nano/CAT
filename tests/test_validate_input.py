"""Tests for :mod:`CAT.data_handling.validate_input`."""

from os.path import join
from shutil import rmtree

import yaml

from scm.plams import (Settings, AMSJob)

from CAT.assertion_functions import assert_eq
from CAT.data_handling.validate_input import validate_input

PATH = 'tests/test_files'


def test_validate_input() -> None:
    """Test :func:`CAT.data_handling.validate_input.validate_input`."""
    with open(join(PATH, 'input1.yaml'), 'r') as f:
        s = Settings(yaml.load(f, Loader=yaml.FullLoader))
    s.path = PATH
    validate_input(s)

    ref = Settings()
    ref.core.dirname = join(PATH, 'core')
    ref.core.dummy = 17

    ref.database.dirname = join(PATH, 'database')
    ref.database.mol_format = ('pdb',)
    ref.database.mongodb = {}
    ref.database.overwrite = ()
    ref.database.read = ('core', 'ligand', 'qd')
    ref.database.write = ('core', 'ligand', 'qd')

    ref.ligand['cosmo-rs'] = False
    ref.ligand.dirname = join(PATH, 'ligand')
    ref.ligand.functional_groups = None
    ref.ligand.optimize = True
    ref.ligand.split = True

    ref.qd.activation_strain = False
    ref.qd.dirname = join(PATH, 'qd')
    ref.qd.dissociate = False
    ref.qd.optimize = {'job1': AMSJob, 's2': {'description': 'UFF with the default forcefield', 'input': {'uff': {'library': 'uff'}, 'ams': {'system': {'bondorders': {'_1': None}, 'charge': 0}}}}, 's1': {'description': 'UFF with the default forcefield', 'input': {'uff': {'library': 'uff'}, 'ams': {'system': {'bondorders': {'_1': None}, 'charge': 0}}}}, 'job2': AMSJob}  # noqa

    try:
        assert_eq(s.optional, ref)
    finally:
        rmtree(join(PATH, 'ligand'))
        rmtree(join(PATH, 'qd'))
        rmtree(join(PATH, 'database'))