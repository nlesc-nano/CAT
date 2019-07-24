"""Tests for :mod:`CAT.data_handling.validation_schemas`."""

from os.path import join

from schema import SchemaError

from scm.plams import AMSJob, ADFJob, Settings

from CAT.utils import get_template
from CAT.assertion_functions import (assert_eq, assert_id, assert_exception)
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

    assert_eq(mol_schema.validate(mol_dict), {'guess_bonds': False})

    mol_dict['guess_bonds'] = 1  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['guess_bonds'] = True  # Correct
    assert_eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['is_core'] = 1  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['is_core'] = True
    assert_eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['column'] = -1  # Exception: value < 0
    assert_exception(*args)
    mol_dict['column'] = 1.0  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['column'] = 1
    assert_eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['row'] = -1  # Exception: value < 0
    assert_exception(*args)
    mol_dict['row'] = 1.0  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['row'] = 1
    assert_eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['indices'] = 1.0  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['indices'] = [1, 5, 6, 7.0]  # Exception: an element has an incorrect type
    assert_exception(*args)
    mol_dict['indices'] = (i for i in range(10))  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['indices'] = -1  # Exception: value < 0
    assert_exception(*args)
    mol_dict['indices'] = [-1, -2, -3, -4, -5]  # Exception: an element is < 0
    assert_exception(*args)
    mol_dict['indices'] = [1, 1, 2]  # Exception: duplicate elements
    assert_exception(*args)

    mol_dict['indices'] = 1
    assert_eq(mol_schema.validate(mol_dict)['indices'], (1,))
    mol_dict['indices'] = [1, 2, 3, 4, 5]
    assert_eq(mol_schema.validate(mol_dict)['indices'], (1, 2, 3, 4, 5))
    mol_dict['indices'] = {1, 2, 3, 4, 5}
    assert_eq(mol_schema.validate(mol_dict)['indices'], (1, 2, 3, 4, 5))
    mol_dict['indices'] = (1, 2, 3, 4, 5)

    mol_dict['type'] = 1  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['type'] = 'bob'
    assert_eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['name'] = 1  # Exception: incorrect type
    assert_exception(*args)
    mol_dict['name'] = 'bob'
    assert_eq(mol_schema.validate(mol_dict), mol_dict)


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

    assert_eq(database_schema.validate(db_dict), ref)

    for key in ('read', 'write', 'overwrite', 'mol_format', 'mongodb'):
        _db_dict = db_dict.copy()
        args = SchemaError, database_schema.validate, _db_dict

        _db_dict[key] = 1  # Exception: incorrect type
        assert_exception(*args)
        _db_dict[key] = 'bob'  # Exception: incorrect value
        assert_exception(*args)
        _db_dict[key] = [1]  # Exception: element has incorrect type
        assert_exception(*args)
        _db_dict[key] = ['bob']  # Exception: element has incorrect value
        assert_exception(*args)

    args = SchemaError, database_schema.validate, db_dict
    db_dict['mongodb'] = True  # Exception: incorrect value
    assert_exception(*args)
    db_dict['mongodb'] = False
    assert_eq(database_schema.validate(db_dict), ref)
    db_dict['mongodb'] = {}
    assert_eq(database_schema.validate(db_dict), ref)


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

    assert_eq(ligand_schema.validate(lig_dict), ref)

    lig_dict['optimize'] = 1  # Exception: incorrect type
    assert_exception(*args)
    lig_dict['optimize'] = True

    lig_dict['split'] = 1  # Exception: incorrect type
    assert_exception(*args)
    lig_dict['split'] = True

    lig_dict['cosmo-rs'] = 1  # Exception: incorrect type
    assert_exception(*args)
    lig_dict['cosmo-rs'] = {}
    assert_eq(ligand_schema.validate(lig_dict)['cosmo-rs'], {})
    lig_dict['cosmo-rs'] = False
    assert_id(ligand_schema.validate(lig_dict)['cosmo-rs'], False)
    lig_dict['cosmo-rs'] = True
    assert_eq(ligand_schema.validate(lig_dict)['cosmo-rs'], {'job1': AMSJob})

    lig_dict['functional_groups'] = 1  # Exception: incorrect type
    assert_exception(*args)
    lig_dict['functional_groups'] = 'CO'
    assert_eq(ligand_schema.validate(lig_dict)['functional_groups'], ('CO',))
    lig_dict['functional_groups'] = ['CO']
    assert_eq(ligand_schema.validate(lig_dict)['functional_groups'], ('CO',))
    lig_dict['functional_groups'] = ['CO', 'CO']  # Exception: duplicate elements
    assert_exception(*args)


def test_core_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.core_schema`."""
    core_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'dummy': 17
    }
    args = SchemaError, core_schema.validate, core_dict

    assert_eq(core_schema.validate(core_dict), ref)

    core_dict['dummy'] = 1.0  # Exception: incorrect type
    assert_exception(*args)
    core_dict['dummy'] = 'H'
    assert_eq(core_schema.validate(core_dict)['dummy'], 1)
    core_dict['dummy'] = 1
    assert_eq(core_schema.validate(core_dict)['dummy'], 1)


def test_qd_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.qd_schema`."""
    qd_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'activation_strain': False,
        'optimize': False,
        'dissociate': False
    }
    args = SchemaError, qd_schema.validate, qd_dict

    assert_eq(qd_schema.validate(qd_dict), ref)

    qd_dict['activation_strain'] = 1  # Exception: incorrect type
    assert_exception(*args)
    qd_dict['activation_strain'] = True

    qd_dict['optimize'] = 1  # Exception: incorrect type
    assert_exception(*args)
    qd_dict['optimize'] = True
    assert_eq(qd_schema.validate(qd_dict)['optimize'], {'job1': AMSJob})
    qd_dict['optimize'] = False
    assert_id(qd_schema.validate(qd_dict)['optimize'], False)

    qd_dict['dissociate'] = 1  # Exception: incorrect type
    assert_exception(*args)
    qd_dict['dissociate'] = True  # Exception: incorrect value
    assert_exception(*args)
    qd_dict['dissociate'] = False
    assert_id(qd_schema.validate(qd_dict)['dissociate'], False)


def test_mongodb_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.mongodb_schema`."""
    mongodb_dict = {}
    ref = {
        'host': 'localhost',
        'port': 27017
    }
    args = SchemaError, mongodb_schema.validate, mongodb_dict

    assert_eq(mongodb_schema.validate(mongodb_dict), ref)

    mongodb_dict['port'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    mongodb_dict['port'] = 27017

    mongodb_dict['host'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    mongodb_dict['host'] = 'localhost'
    assert_eq(mongodb_schema.validate(mongodb_dict)['host'], 'localhost')
    mongodb_dict['host'] = 51
    assert_eq(mongodb_schema.validate(mongodb_dict)['host'], 51)

    mongodb_dict['username'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    mongodb_dict['username'] = 'bob'
    assert_eq(mongodb_schema.validate(mongodb_dict)['username'], 'bob')
    mongodb_dict['username'] = 52
    assert_eq(mongodb_schema.validate(mongodb_dict)['username'], 52)

    mongodb_dict['password'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    mongodb_dict['password'] = 'secret'
    assert_eq(mongodb_schema.validate(mongodb_dict)['password'], 'secret')
    mongodb_dict['password'] = 53
    assert_eq(mongodb_schema.validate(mongodb_dict)['password'], 53)


def test_qd_opt_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.qd_opt_schema`."""
    _qd_opt_s1_default = get_template('qd.yaml')['UFF']
    _qd_opt_s2_default = _qd_opt_s1_default

    qd_opt_dict = Settings()
    ref = Settings({
        'job1': AMSJob,
        's1': _qd_opt_s1_default,
        'job2': AMSJob,
        's2': _qd_opt_s2_default
    })
    args = SchemaError, qd_opt_schema.validate, qd_opt_dict

    assert_eq(qd_opt_schema.validate(qd_opt_dict), ref)

    for job in ('job1', 'job2'):
        qd_opt_dict[job] = 1  # Exception: incorrect type
        assert_exception(*args)
        qd_opt_dict[job] = int  # Exception: incorrect value
        assert_exception(*args)
        qd_opt_dict[job] = 'bob'  # Exception: incorrect value
        assert_exception(*args)
        qd_opt_dict[job] = 'ADFJob'
        assert_id(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)
        qd_opt_dict[job] = 'ADFJOB'
        assert_id(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)
        qd_opt_dict[job] = ADFJob
        assert_id(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        qd_opt_dict[s] = 1  # Exception: incorrect type
        assert_exception(*args)
        qd_opt_dict[s] = {'key1': {'key2': {'key3': True}}}
        assert_eq(qd_opt_schema.validate(qd_opt_dict)[s], ref)
        qd_opt_dict[s] = join(PATH, 'settings.yaml')
        assert_eq(qd_opt_schema.validate(qd_opt_dict)[s], ref)


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

    assert_eq(crs_schema.validate(crs_dict), ref)

    crs_dict['keep_files'] = 1  # Exception: incorrect type
    assert_exception(*args)
    crs_dict['keep_files'] = False
    assert_id(crs_schema.validate(crs_dict)['keep_files'], False)
    crs_dict['keep_files'] = True

    for job in ('job1', 'job2'):
        crs_dict[job] = 1  # Exception: incorrect type
        assert_exception(*args)
        crs_dict[job] = int  # Exception: incorrect value
        assert_exception(*args)
        crs_dict[job] = 'bob'  # Exception: incorrect value
        assert_exception(*args)
        crs_dict[job] = 'ADFJob'
        assert_id(crs_schema.validate(crs_dict)[job], ADFJob)
        crs_dict[job] = 'ADFJOB'
        assert_id(crs_schema.validate(crs_dict)[job], ADFJob)
        crs_dict[job] = ADFJob
        assert_id(crs_schema.validate(crs_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        crs_dict[s] = 1  # Exception: incorrect type
        assert_exception(*args)
        crs_dict[s] = {'key1': {'key2': {'key3': True}}}
        assert_eq(crs_schema.validate(crs_dict)[s], ref)
        crs_dict[s] = join(PATH, 'settings.yaml')
        assert_eq(crs_schema.validate(crs_dict)[s], ref)


def test_bde_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.bde_schema`."""
    _bde_s1_default = get_template('qd.yaml')['MOPAC']
    _bde_s2_default = get_template('qd.yaml')['UFF']

    bde_dict = Settings({'core_atom': 'Cd', 'lig_count': 2})
    ref = Settings({
        'keep_files': True,
        'core_atom': 48,
        'lig_count': 2,
        'core_core_dist': 0.0,
        'lig_core_dist': 5.0,
        'topology': {},
        'job1': AMSJob,
        's1': _bde_s1_default
    })
    args = SchemaError, bde_schema.validate, bde_dict

    assert_eq(bde_schema.validate(bde_dict), ref)

    bde_dict['keep_files'] = 1  # Exception: incorrect type
    assert_exception(*args)
    bde_dict['keep_files'] = False
    assert_id(bde_schema.validate(bde_dict)['keep_files'], False)
    bde_dict['keep_files'] = True

    bde_dict['core_atom'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    bde_dict['core_atom'] = 'H'
    assert_eq(bde_schema.validate(bde_dict)['core_atom'], 1)
    bde_dict['core_atom'] = 1
    assert_eq(bde_schema.validate(bde_dict)['core_atom'], 1)

    bde_dict['lig_count'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    bde_dict['lig_count'] = -1  # Exception: incorrect value
    assert_exception(*args)
    bde_dict['lig_count'] = 3
    assert_eq(bde_schema.validate(bde_dict)['lig_count'], 3)

    bde_dict['core_index'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    bde_dict['core_index'] = [1, 2, 3, 4, 5.0]  # Exception: incorrect element type
    assert_exception(*args)
    bde_dict['core_index'] = [1, 2, 3, 4, 4]  # Exception: duplicate elements
    assert_exception(*args)
    bde_dict['core_index'] = 1
    assert_eq(bde_schema.validate(bde_dict)['core_index'], (1,))
    bde_dict['core_index'] = [1, 2, 3]
    assert_eq(bde_schema.validate(bde_dict)['core_index'], (1, 2, 3))
    bde_dict['core_index'] = {1, 2, 3}
    assert_eq(bde_schema.validate(bde_dict)['core_index'], (1, 2, 3))

    bde_dict['topology'] = 5.0  # Exception: incorrect type
    assert_exception(*args)
    bde_dict['topology'] = {'key': 'value'}  # Exception: incorrect value
    assert_exception(*args)
    bde_dict['topology'] = {1: 'value'}
    assert_eq(bde_schema.validate(bde_dict)['topology'], {1: 'value'})

    for dist in ('core_core_dist', 'lig_core_dist'):
        bde_dict[dist] = 'bob'  # Exception: incorrect type
        assert_exception(*args)
        bde_dict[dist] = -1  # Exception: incorrect value
        assert_exception(*args)
        bde_dict[dist] = 4
        assert_eq(bde_schema.validate(bde_dict)[dist], 4.0)
        bde_dict[dist] = 4.0
        assert_eq(bde_schema.validate(bde_dict)[dist], 4.0)

    for job in ('job1', 'job2'):
        bde_dict[job] = 1  # Exception: incorrect type
        assert_exception(*args)
        bde_dict[job] = int  # Exception: incorrect value
        assert_exception(*args)
        bde_dict[job] = 'bob'  # Exception: incorrect value
        assert_exception(*args)
        bde_dict[job] = 'ADFJob'
        assert_id(bde_schema.validate(bde_dict)[job], ADFJob)
        bde_dict[job] = 'ADFJOB'
        assert_id(bde_schema.validate(bde_dict)[job], ADFJob)
        bde_dict[job] = ADFJob
        assert_id(bde_schema.validate(bde_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        bde_dict[s] = 1  # Exception: incorrect type
        assert_exception(*args)
        bde_dict[s] = {'key1': {'key2': {'key3': True}}}
        assert_eq(bde_schema.validate(bde_dict)[s], ref)
        bde_dict[s] = join(PATH, 'settings.yaml')
        assert_eq(bde_schema.validate(bde_dict)[s], ref)
