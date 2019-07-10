"""Tests for :mod:`CAT.data_handling.validation_schemas`."""

from schema import SchemaError

from scm.plams import Settings

from CAT.assertion_functions import (assert_eq, assert_exception)
from CAT.data_handling.validation_schemas import (
    mol_schema, core_schema, ligand_schema, qd_schema, database_schema,
    mongodb_schema, bde_schema, qd_opt_schema, crs_schema
)


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
        'overwrite': ('core', 'ligand', 'qd'),
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

"""
ligand_schema = Schema({
    # path+directory name of the ligand directory
    'dirname':
        str,

    Optional_('functional_groups', default=None):
        Or(
            And(str, Use(lambda n: (n,))),
            And(abc.Collection, lambda n: all(isinstance(i, str) for i in n), Use(tuple))
        ),

    Optional_('optimize', default=True):  # Optimize the ligands
        bool,

    Optional_('split', default=True):  # Remove a counterion from the function group
        bool,

    Optional_('cosmo-rs', default=False):  # Settings specific to ligand COSMO-RS calculations
        Or(
            dict,
            And(bool, Use(lambda n: {'job1': AMSJob} if n else False))
        ),
})
"""