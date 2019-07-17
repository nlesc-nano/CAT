"""
CAT.data_handling.validation_schemas
====================================

A module designed for sanitizing and interpreting the input file.

Index
-----
.. currentmodule:: CAT.data_handling.validation_schemas
.. autosummary::
    mol_schema
    core_schema
    ligand_schema
    qd_schema
    database_schema
    mongodb_schema
    bde_schema
    qd_opt_schema
    crs_schema

API
---
.. autodata:: mol_schema
    :annotation: = schema.Schema
.. autodata:: core_schema
    :annotation: = schema.Schema
.. autodata:: ligand_schema
    :annotation: = schema.Schema
.. autodata:: qd_schema
    :annotation: = schema.Schema
.. autodata:: database_schema
    :annotation: = schema.Schema
.. autodata:: mongodb_schema
    :annotation: = schema.Schema
.. autodata:: bde_schema
    :annotation: = schema.Schema
.. autodata:: qd_opt_schema
    :annotation: = schema.Schema
.. autodata:: crs_schema
    :annotation: = schema.Schema

"""

from typing import (Dict, Collection)
from collections import abc

from schema import (Or, And, Use, Schema)
from schema import Optional as Optional_

from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.uff import UFFJob
from scm.plams.interfaces.adfsuite.band import BANDJob
from scm.plams.interfaces.adfsuite.dftb import DFTBJob
from scm.plams.interfaces.adfsuite.mopac import MOPACJob
from scm.plams.interfaces.adfsuite.reaxff import ReaxFFJob

from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob
from scm.plams.interfaces.thirdparty.dftbplus import DFTBPlusJob

from scm.plams.core.basejob import Job

from ..utils import get_template, validate_path
from ..mol_utils import to_atnum

try:
    from nanoCAT.crs import CRSJob
    NANO_CAT = True
except ImportError:
    CRSJob = Job
    NANO_CAT = False

__all__ = ['mol_schema', 'core_schema', 'ligand_schema', 'qd_schema', 'database_schema',
           'mongodb_schema', 'bde_schema', 'qd_opt_schema', 'crs_schema']


def to_tuple(collection: Collection) -> tuple:
    """Convert a collection into a sorted tuple."""
    try:
        ret = sorted(collection)
    except TypeError:  # The collection contains a mix of sorting-incompatibl objects
        ret = sorted(collection, key=str)
    finally:
        return tuple(ret)


# The **default** parameter of schema.Optional() will automatically call any callable
# Solution: provide a callable that returns another callable
def _get_amsjob() -> type:
    """Return a type object of :class:`.AMSJob`."""
    return AMSJob


def _get_crsjob() -> type:
    """Return a type object of :class:`.CRSJob`."""
    return CRSJob


# Default settings templates
_bde_s1_default = get_template('qd.yaml')['MOPAC']
_bde_s2_default = get_template('qd.yaml')['UFF']

_qd_opt_s1_default = get_template('qd.yaml')['UFF']
_qd_opt_s2_default = _qd_opt_s1_default

_crs_s1_default = get_template('qd.yaml')['COSMO-MOPAC']
_crs_s2_default = get_template('qd.yaml')['COSMO-RS activity coefficient']
_crs_s2_default.update(get_template('crs.yaml')['MOPAC PM6'])


# A dictionary for translating strings into :class:`plams.Job` types
_class_dict: Dict[str, type] = {
    'adf': ADFJob, 'adfjob': ADFJob,
    'ams': AMSJob, 'amsjob': AMSJob,
    'uff': UFFJob, 'uffjob': UFFJob,
    'band': BANDJob, 'bandjob': BANDJob,
    'dftb': DFTBJob, 'dftbjob': DFTBJob,
    'mopac': MOPACJob, 'mopacjob': MOPACJob,
    'reaxff': ReaxFFJob, 'reaxffjob': ReaxFFJob,
    'cp2k': Cp2kJob, 'cp2kjob': Cp2kJob,
    'orca': ORCAJob, 'orcajob': ORCAJob,
    'dirac': DiracJob, 'diracjob': DiracJob,
    'gamess': GamessJob, 'gamessjob': GamessJob,
    'dftbplus': DFTBPlusJob, 'dftbplusjob': DFTBPlusJob,
    'crs': CRSJob, 'cosmo-rs': CRSJob, 'crsjob': CRSJob
}


#: Schema for validating the ``['input_ligands']`` and ``['input_cores']`` blocks.
mol_schema: Schema = Schema({
    Optional_('guess_bonds', default=False):
        bool,

    Optional_('is_core'):
        bool,

    Optional_('column'):
        And(int, lambda n: n >= 0),

    Optional_('row'):
        And(int, lambda n: n >= 0),

    Optional_('indices'):
        Or(
            And(int, lambda n: n >= 0, Use(lambda n: (n,))),
            And(
                abc.Collection,
                lambda n: all(isinstance(i, int) and i >= 0 for i in n),
                lambda n: len(n) == len(set(n)),
                Use(tuple)
            ),
        ),

    Optional_('type'):
        str,

    Optional_('name'):
        str,

    Optional_('path'):
        Use(validate_path)
})

#: Schema for validating the ``['optional']['core']`` block.
core_schema: Schema = Schema({
    'dirname':
        str,

    Optional_('dummy', default=17):  # Return a tuple of atomic numbers
        Or(
            And(int, Use(to_atnum)),
            And(str, Use(to_atnum))
        )
})

_db_names = ('core', 'ligand', 'qd')
_format_names = ('pdb', 'xyz')

#: Schema for validating the ``['optional']['database']`` block.
database_schema: Schema = Schema({
    # path+directory name of the database
    'dirname':
        str,

    Optional_('read', default=_db_names):  # Attempt to pull structures from the database
        Or(
            And(bool, Use(lambda n: _db_names if n is True else ())),
            And(str, lambda n: n in _db_names, Use(lambda n: (n,))),
            And(abc.Collection,
                lambda n: all(i in _db_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple))
        ),

    Optional_('write', default=_db_names):  # Attempt to write structures to the database
        Or(
            And(bool, Use(lambda n: _db_names if n is True else ())),
            And(str, lambda n: n in _db_names, Use(lambda n: (n,))),
            And(abc.Collection,
                lambda n: all(i in _db_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple))
        ),

    Optional_('overwrite', default=tuple):  # Allow previous entries to be overwritten
        Or(
            And(bool, Use(lambda n: _db_names if n is True else ())),
            And(str, lambda n: n in _db_names, Use(lambda n: (n,))),
            And(abc.Collection,
                lambda n: all(i in _db_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple))
        ),

    Optional_('mongodb', default=dict):  # Settings specific to MongoDB
        Or(
            dict,
            And(bool, lambda n: n is False, Use(lambda n: {}))
        ),

    Optional_('mol_format', default=_format_names):  # Return a tuple of file formats
        Or(
            And(bool, Use(lambda n: _format_names if n is True else ())),
            And(str, lambda n: n in _format_names),
            And(abc.Collection,
                lambda n: all(i in _format_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple))
        )
})


#: Schema for validating the ``['optional']['ligand']`` block.
ligand_schema: Schema = Schema({
    # path+directory name of the ligand directory
    'dirname':
        str,

    Optional_('functional_groups', default=None):
        Or(
            And(str, Use(lambda n: (n,))),
            And(abc.Collection,
                lambda n: all(isinstance(i, str) for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple))
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


#: Schema for validating the ``['optional']['qd']`` block.
qd_schema: Schema = Schema({
    # path+directory name of the quantum dot directory
    'dirname':
        str,

    # Settings specific to a quantum dot activation strain analyses
    Optional_('activation_strain', default=False):
        bool,

    Optional_('optimize', default=False):  # Settings for quantum dot geometry optimizations
        Or(
            dict,
            And(bool, Use(lambda n: ({'job1': AMSJob} if n else False)))
        ),

    # Settings for quantum dot ligand dissociation calculations
    Optional_('dissociate', default=False):
        Or(
            dict,
            And(bool, lambda n: n is False)
        )
})


#: Schema for validating the ``['optional']['database']['mongodb']`` block.
mongodb_schema: Schema = Schema({
    # Optional username for the MongoDB host
    Optional_('username'):
        Or(str, int),

    Optional_('password'):  # Optional password for the MongoDB host
        Or(str, int),

    Optional_('host', default='localhost'):  # Name of the MongoDB host
        Or(str, int),

    Optional_('port', default=27017):  # Port of the MongoDB host
        int,

    Optional_(str):  # Other keyword arguments for :class:`pymongo.MongoClient`
        object
})


#: Schema for validating the ``['optional']['qd']['dissociate']`` block.
bde_schema: Schema = Schema({
    # Atom type of the to-be dissociated core atom
    'core_atom':
        And(Or(int, str), Use(to_atnum)),

    'lig_count':  # The number of ligands per core_atom
        And(int, lambda n: n >= 0),

    Optional_('keep_files', default=True):  # Delete files after the calculations are finished
        bool,

    Optional_('core_core_dist', default=0.0):
        And(Or(int, float), lambda n: n >= 0.0, Use(float)),

    Optional_('lig_core_dist', default=5.0):
        And(Or(int, float), lambda n: n >= 0.0, Use(float)),

    Optional_('core_index'):
        Or(
            And(int, lambda n: n >= 0, Use(lambda n: (n,))),
            And(
                abc.Collection,
                lambda n: all(isinstance(i, int) and i >= 0 for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple)
            )
        ),

    Optional_('topology', default=dict):
        And(dict, lambda n: all(isinstance(k, int) for k in n)),

    Optional_('job1', default=_get_amsjob):
        Or(
            And(type, lambda n: issubclass(n, Job)),
            And(str, lambda n: n.lower() in _class_dict, Use(lambda n: _class_dict[n.lower()]))
        ),

    Optional_('s1', default=_bde_s1_default):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False)))
        ),

    Optional_('job2'):
        Or(
            And(type, lambda n: issubclass(n, Job)),
            And(str, lambda n: n.lower() in _class_dict, Use(lambda n: _class_dict[n.lower()]))
        ),

    Optional_('s2'):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False)))
        )
})

#: Schema for validating the ``['optional']['qd']['optimize']`` block.
qd_opt_schema: Schema = Schema({
    # The job type for the first half of the optimization
    Optional_('job1', default=_get_amsjob):
        Or(
            And(type, lambda n: issubclass(n, Job)),
            And(str, lambda n: n.lower() in _class_dict, Use(lambda n: _class_dict[n.lower()]))
        ),

    # The job settings for the first half of the optimization
    Optional_('s1', default=_qd_opt_s1_default):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False)))
        ),

    # The job type for the second half of the optimization
    Optional_('job2', default=_get_amsjob):
        Or(
            And(type, lambda n: issubclass(n, Job)),
            And(str, lambda n: n.lower() in _class_dict, Use(lambda n: _class_dict[n.lower()]))
        ),

    # The job settings for the second half of the optimization
    Optional_('s2', default=_qd_opt_s2_default):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False)))
        )
})

#: Schema for validating the ``['optional']['ligand']['cosmo-rs']`` block.
crs_schema: Schema = Schema({
    # Delete files after the calculations are finished
    Optional_('keep_files', default=True):
        bool,

    # The job type for constructing the COSMO surface
    Optional_('job1', default=_get_amsjob):
        Or(
            And(type, lambda n: issubclass(n, Job)),
            And(str, lambda n: n.lower() in _class_dict, Use(lambda n: _class_dict[n.lower()]))
        ),

    # The settings for constructing the COSMO surface
    Optional_('s1', default=_crs_s1_default):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False)))
        ),

    Optional_('job2', default=_get_crsjob):  # The job type for the actual COSMO-RS calculation
        Or(
            And(type, lambda n: issubclass(n, Job)),
            And(str, lambda n: n.lower() in _class_dict, Use(lambda n: _class_dict[n.lower()]))
        ),

    Optional_('s2', default=_crs_s2_default):  # The settings for the actual COSMO-RS calculation
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False)))
        )
})
