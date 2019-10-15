"""Tests for :mod:`CAT.data_handling.validation_schemas`."""

import os
from os.path import join

from schema import SchemaError
from unittest import mock

from scm.plams import AMSJob, ADFJob, Settings

from CAT.utils import get_template
from CAT.assertion.assertion_manager import assertion
from CAT.data_handling.validation_schemas import (
    mol_schema, core_schema, ligand_schema, qd_schema, database_schema,
    mongodb_schema, bde_schema, qd_opt_schema, crs_schema
)

from nanoCAT.crs import CRSJob

PATH = join('tests', 'test_files')


def test_mol_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.mol_schema`."""
    mol_dict = {}
    args = SchemaError, mol_schema.validate, mol_dict

    assertion.eq(mol_schema.validate(mol_dict), {'guess_bonds': False})

    mol_dict['guess_bonds'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['guess_bonds'] = True  # Correct
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['is_core'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['is_core'] = True
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['column'] = -1  # Exception: value < 0
    assertion.exception(*args)
    mol_dict['column'] = 1.0  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['column'] = 1
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['row'] = -1  # Exception: value < 0
    assertion.exception(*args)
    mol_dict['row'] = 1.0  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['row'] = 1
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['indices'] = 1.0  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['indices'] = [1, 5, 6, 7.0]  # Exception: an element has an incorrect type
    assertion.exception(*args)
    mol_dict['indices'] = (i for i in range(10))  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['indices'] = -1  # Exception: value < 0
    assertion.exception(*args)
    mol_dict['indices'] = [-1, -2, -3, -4, -5]  # Exception: an element is < 0
    assertion.exception(*args)
    mol_dict['indices'] = [1, 1, 2]  # Exception: duplicate elements
    assertion.exception(*args)

    mol_dict['indices'] = 1
    assertion.eq(mol_schema.validate(mol_dict)['indices'], (1,))
    mol_dict['indices'] = [1, 2, 3, 4, 5]
    assertion.eq(mol_schema.validate(mol_dict)['indices'], (1, 2, 3, 4, 5))
    mol_dict['indices'] = {1, 2, 3, 4, 5}
    assertion.eq(mol_schema.validate(mol_dict)['indices'], (1, 2, 3, 4, 5))
    mol_dict['indices'] = (1, 2, 3, 4, 5)

    mol_dict['type'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['type'] = 'bob'
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['name'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    mol_dict['name'] = 'bob'
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)


def test_database_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.database_schema`."""
    db_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'read': ('core', 'ligand', 'qd'),
        'write': ('core', 'ligand', 'qd'),
        'overwrite': (),
        'mongodb': {},
        'mol_format': ('pdb', 'xyz')

    }

    assertion.eq(database_schema.validate(db_dict), ref)

    for key in ('read', 'write', 'overwrite', 'mol_format', 'mongodb'):
        _db_dict = db_dict.copy()
        args = SchemaError, database_schema.validate, _db_dict

        _db_dict[key] = 1  # Exception: incorrect type
        assertion.exception(*args)
        _db_dict[key] = 'bob'  # Exception: incorrect value
        assertion.exception(*args)
        _db_dict[key] = [1]  # Exception: element has incorrect type
        assertion.exception(*args)
        _db_dict[key] = ['bob']  # Exception: element has incorrect value
        assertion.exception(*args)

    args = SchemaError, database_schema.validate, db_dict
    db_dict['mongodb'] = True  # Exception: incorrect value
    assertion.exception(*args)
    db_dict['mongodb'] = False
    assertion.eq(database_schema.validate(db_dict), ref)
    db_dict['mongodb'] = {}
    assertion.eq(database_schema.validate(db_dict), ref)


def test_ligand_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.ligand_schema`."""
    lig_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'functional_groups': None,
        'optimize': True,
        'split': True,
        'cosmo-rs': False
    }
    args = SchemaError, ligand_schema.validate, lig_dict

    assertion.eq(ligand_schema.validate(lig_dict), ref)

    lig_dict['optimize'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    lig_dict['optimize'] = True

    lig_dict['split'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    lig_dict['split'] = True

    lig_dict['cosmo-rs'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    lig_dict['cosmo-rs'] = {}
    assertion.eq(ligand_schema.validate(lig_dict)['cosmo-rs'], {})
    lig_dict['cosmo-rs'] = False
    assertion.is_(ligand_schema.validate(lig_dict)['cosmo-rs'], False)
    lig_dict['cosmo-rs'] = True
    assertion.eq(ligand_schema.validate(lig_dict)['cosmo-rs'], {'job1': 'AMSJob'})

    lig_dict['functional_groups'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    lig_dict['functional_groups'] = 'CO'
    assertion.eq(ligand_schema.validate(lig_dict)['functional_groups'], ('CO',))
    lig_dict['functional_groups'] = ['CO']
    assertion.eq(ligand_schema.validate(lig_dict)['functional_groups'], ('CO',))
    lig_dict['functional_groups'] = ['CO', 'CO']  # Exception: duplicate elements
    assertion.exception(*args)


def test_core_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.core_schema`."""
    core_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'dummy': 17
    }
    args = SchemaError, core_schema.validate, core_dict

    assertion.eq(core_schema.validate(core_dict), ref)

    core_dict['dummy'] = 1.0  # Exception: incorrect type
    assertion.exception(*args)
    core_dict['dummy'] = 'H'
    assertion.eq(core_schema.validate(core_dict)['dummy'], 1)
    core_dict['dummy'] = 1
    assertion.eq(core_schema.validate(core_dict)['dummy'], 1)


def test_qd_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.qd_schema`."""
    qd_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'activation_strain': False,
        'optimize': False,
        'dissociate': False,
        'bulkiness': False
    }
    args = SchemaError, qd_schema.validate, qd_dict

    assertion.eq(qd_schema.validate(qd_dict), ref)

    qd_dict['activation_strain'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    qd_dict['activation_strain'] = True

    qd_dict['bulkiness'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    qd_dict['bulkiness'] = False

    qd_dict['optimize'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    qd_dict['optimize'] = True
    assertion.eq(qd_schema.validate(qd_dict)['optimize'], {'job1': 'AMSJob'})
    qd_dict['optimize'] = False
    assertion.is_(qd_schema.validate(qd_dict)['optimize'], False)

    qd_dict['dissociate'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    qd_dict['dissociate'] = True  # Exception: incorrect value
    assertion.exception(*args)
    qd_dict['dissociate'] = False
    assertion.is_(qd_schema.validate(qd_dict)['dissociate'], False)


def test_mongodb_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.mongodb_schema`."""
    mongodb_dict = {}
    ref = {
        'host': 'localhost',
        'port': 27017
    }
    args = SchemaError, mongodb_schema.validate, mongodb_dict

    assertion.eq(mongodb_schema.validate(mongodb_dict), ref)

    mongodb_dict['port'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    mongodb_dict['port'] = 27017

    mongodb_dict['host'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    mongodb_dict['host'] = 'localhost'
    assertion.eq(mongodb_schema.validate(mongodb_dict)['host'], 'localhost')
    mongodb_dict['host'] = 51
    assertion.eq(mongodb_schema.validate(mongodb_dict)['host'], 51)

    mongodb_dict['username'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    mongodb_dict['username'] = 'bob'
    assertion.eq(mongodb_schema.validate(mongodb_dict)['username'], 'bob')
    mongodb_dict['username'] = 52
    assertion.eq(mongodb_schema.validate(mongodb_dict)['username'], 52)

    mongodb_dict['password'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    mongodb_dict['password'] = 'secret'
    assertion.eq(mongodb_schema.validate(mongodb_dict)['password'], 'secret')
    mongodb_dict['password'] = 53
    assertion.eq(mongodb_schema.validate(mongodb_dict)['password'], 53)


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
def test_qd_opt_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.qd_opt_schema`."""
    _qd_opt_s1_default = get_template('qd.yaml')['UFF']
    _qd_opt_s2_default = _qd_opt_s1_default

    qd_opt_dict = Settings()
    ref = Settings({
        'job1': AMSJob,
        's1': _qd_opt_s1_default,
        'job2': AMSJob,
        's2': _qd_opt_s2_default,
        'use_ff': False
    })
    args = SchemaError, qd_opt_schema.validate, qd_opt_dict

    assertion.eq(qd_opt_schema.validate(qd_opt_dict), ref)

    for job in ('job1', 'job2'):
        qd_opt_dict[job] = 1  # Exception: incorrect type
        assertion.exception(*args)
        qd_opt_dict[job] = int  # Exception: incorrect value
        assertion.exception(*args)
        qd_opt_dict[job] = 'bob'  # Exception: incorrect value
        assertion.exception(*args)
        qd_opt_dict[job] = 'ADFJob'
        assertion.is_(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)
        qd_opt_dict[job] = 'ADFJOB'
        assertion.is_(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)
        qd_opt_dict[job] = ADFJob
        assertion.is_(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        qd_opt_dict[s] = 1  # Exception: incorrect type
        assertion.exception(*args)
        qd_opt_dict[s] = {'key1': {'key2': {'key3': True}}}
        assertion.eq(qd_opt_schema.validate(qd_opt_dict)[s], ref)
        qd_opt_dict[s] = join(PATH, 'settings.yaml')
        assertion.eq(qd_opt_schema.validate(qd_opt_dict)[s], ref)


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
def test_crs_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.crs_schema`."""
    _crs_s1_default = get_template('qd.yaml')['COSMO-MOPAC']
    _crs_s2_default = get_template('qd.yaml')['COSMO-RS activity coefficient']
    _crs_s2_default.update(get_template('crs.yaml')['MOPAC PM6'])

    crs_dict = Settings()
    ref = Settings({
        'keep_files': True,
        'job1': AMSJob,
        's1': _crs_s1_default,
        'job2': CRSJob,
        's2': _crs_s2_default
    })
    args = SchemaError, crs_schema.validate, crs_dict

    assertion.eq(crs_schema.validate(crs_dict), ref)

    crs_dict['keep_files'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    crs_dict['keep_files'] = False
    assertion.is_(crs_schema.validate(crs_dict)['keep_files'], False)
    crs_dict['keep_files'] = True

    for job in ('job1', 'job2'):
        crs_dict[job] = 1  # Exception: incorrect type
        assertion.exception(*args)
        crs_dict[job] = int  # Exception: incorrect value
        assertion.exception(*args)
        crs_dict[job] = 'bob'  # Exception: incorrect value
        assertion.exception(*args)
        crs_dict[job] = 'ADFJob'
        assertion.is_(crs_schema.validate(crs_dict)[job], ADFJob)
        crs_dict[job] = 'ADFJOB'
        assertion.is_(crs_schema.validate(crs_dict)[job], ADFJob)
        crs_dict[job] = ADFJob
        assertion.is_(crs_schema.validate(crs_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        crs_dict[s] = 1  # Exception: incorrect type
        assertion.exception(*args)
        crs_dict[s] = {'key1': {'key2': {'key3': True}}}
        assertion.eq(crs_schema.validate(crs_dict)[s], ref)
        crs_dict[s] = join(PATH, 'settings.yaml')
        assertion.eq(crs_schema.validate(crs_dict)[s], ref)


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
def test_bde_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.bde_schema`."""
    _bde_s1_default = get_template('qd.yaml')['MOPAC']
    _bde_s2_default = get_template('qd.yaml')['UFF']

    bde_dict = Settings({'core_atom': 'Cd', 'lig_count': 2})
    ref = Settings({
        'keep_files': True,
        'use_ff': False,
        'core_atom': 48,
        'lig_count': 2,
        'core_core_dist': 0.0,
        'lig_core_dist': 5.0,
        'topology': {},
        'job1': AMSJob,
        's1': _bde_s1_default
    })
    args = SchemaError, bde_schema.validate, bde_dict

    assertion.eq(bde_schema.validate(bde_dict), ref)

    bde_dict['keep_files'] = 1  # Exception: incorrect type
    assertion.exception(*args)
    bde_dict['keep_files'] = False
    assertion.is_(bde_schema.validate(bde_dict)['keep_files'], False)
    bde_dict['keep_files'] = True

    bde_dict['core_atom'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    bde_dict['core_atom'] = 'H'
    assertion.eq(bde_schema.validate(bde_dict)['core_atom'], 1)
    bde_dict['core_atom'] = 1
    assertion.eq(bde_schema.validate(bde_dict)['core_atom'], 1)

    bde_dict['lig_count'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    bde_dict['lig_count'] = -1  # Exception: incorrect value
    assertion.exception(*args)
    bde_dict['lig_count'] = 3
    assertion.eq(bde_schema.validate(bde_dict)['lig_count'], 3)

    bde_dict['core_index'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    bde_dict['core_index'] = [1, 2, 3, 4, 5.0]  # Exception: incorrect element type
    assertion.exception(*args)
    bde_dict['core_index'] = [1, 2, 3, 4, 4]  # Exception: duplicate elements
    assertion.exception(*args)
    bde_dict['core_index'] = 1
    assertion.eq(bde_schema.validate(bde_dict)['core_index'], (1,))
    bde_dict['core_index'] = [1, 2, 3]
    assertion.eq(bde_schema.validate(bde_dict)['core_index'], (1, 2, 3))
    bde_dict['core_index'] = {1, 2, 3}
    assertion.eq(bde_schema.validate(bde_dict)['core_index'], (1, 2, 3))

    bde_dict['topology'] = 5.0  # Exception: incorrect type
    assertion.exception(*args)
    bde_dict['topology'] = {'key': 'value'}  # Exception: incorrect value
    assertion.exception(*args)
    bde_dict['topology'] = {1: 'value'}
    assertion.eq(bde_schema.validate(bde_dict)['topology'], {1: 'value'})

    for dist in ('core_core_dist', 'lig_core_dist'):
        bde_dict[dist] = 'bob'  # Exception: incorrect type
        assertion.exception(*args)
        bde_dict[dist] = -1  # Exception: incorrect value
        assertion.exception(*args)
        bde_dict[dist] = 4
        assertion.eq(bde_schema.validate(bde_dict)[dist], 4.0)
        bde_dict[dist] = 4.0
        assertion.eq(bde_schema.validate(bde_dict)[dist], 4.0)

    for job in ('job1', 'job2'):
        bde_dict[job] = 1  # Exception: incorrect type
        assertion.exception(*args)
        bde_dict[job] = int  # Exception: incorrect value
        assertion.exception(*args)
        bde_dict[job] = 'bob'  # Exception: incorrect value
        assertion.exception(*args)
        bde_dict[job] = 'ADFJob'
        assertion.is_(bde_schema.validate(bde_dict)[job], ADFJob)
        bde_dict[job] = 'ADFJOB'
        assertion.is_(bde_schema.validate(bde_dict)[job], ADFJob)
        bde_dict[job] = ADFJob
        assertion.is_(bde_schema.validate(bde_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        bde_dict[s] = 1  # Exception: incorrect type
        assertion.exception(*args)
        bde_dict[s] = {'key1': {'key2': {'key3': True}}}
        assertion.eq(bde_schema.validate(bde_dict)[s], ref)
        bde_dict[s] = join(PATH, 'settings.yaml')
        assertion.eq(bde_schema.validate(bde_dict)[s], ref)
