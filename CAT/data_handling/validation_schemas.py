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
    _class_dict
    _class_dict_scm

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
.. autodata:: _class_dict
    :annotation: = dict[str, type]
.. autodata:: _class_dict_scm
    :annotation: = dict[str, type]

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

from ..utils import (get_template, validate_path, validate_core_atom, check_sys_var)
from ..mol_utils import to_atnum

try:
    from nanoCAT.crs import CRSJob
    NANO_CAT = True
except ImportError:
    CRSJob = Job
    NANO_CAT = False

__all__ = ['mol_schema', 'core_schema', 'ligand_schema', 'qd_schema', 'database_schema',
           'mongodb_schema', 'bde_schema', 'qd_opt_schema', 'crs_schema']


def val_job_type(value: type) -> type:
    """Call :func:`.check_sys_var` if value is in :data:`._class_dict_scm`"""
    if value in _class_dict_scm.values():
        check_sys_var()
    return value


def str_to_job_type(key: str) -> type:
    """Convert a string into a type object.

    Parameters
    ----------
    key : str
        An alias for a :class:`type` object of a PLAMS :class:`Job`.
        **key** is converted by passing it into :data:`._class_dict` or :data:`._class_dict_scm`.
        The :func:`check_sys_var` function is called if **key** is present
        in :data:`._class_dict_scm`.

    Returns
    -------
    |type|_
        A type object constructed from **key**.

    Raises
    ------
    KeyError
        Raised if the decapitalized **key** is available in neither :data:`._class_dict` nor
        :data:`._class_dict_scm`.

    """
    _key = key.lower()
    if _key in _class_dict:
        return _class_dict[_key]
    elif _key in _class_dict_scm:
        check_sys_var()
        return _class_dict_scm[_key]
    raise KeyError(f'No Job type alias available for {repr(_key)}')


def to_tuple(collection: Collection) -> tuple:
    """Convert a collection to a sorted tuple."""
    try:
        ret = sorted(collection)
    except TypeError:  # The collection contains a mix of sorting-incompatible objects
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


#: A dictionary for translating strings into :class:`plams.Job` types for third party software
_class_dict: Dict[str, type] = {
    'cp2k': Cp2kJob, 'cp2kjob': Cp2kJob,
    'orca': ORCAJob, 'orcajob': ORCAJob,
    'dirac': DiracJob, 'diracjob': DiracJob,
    'gamess': GamessJob, 'gamessjob': GamessJob,
    'dftbplus': DFTBPlusJob, 'dftbplusjob': DFTBPlusJob,
}

#: A dictionary for translating strings into :class:`plams.Job` types for ADF
_class_dict_scm: Dict[str, type] = {
    'adf': ADFJob, 'adfjob': ADFJob,
    'ams': AMSJob, 'amsjob': AMSJob,
    'uff': UFFJob, 'uffjob': UFFJob,
    'band': BANDJob, 'bandjob': BANDJob,
    'dftb': DFTBJob, 'dftbjob': DFTBJob,
    'mopac': MOPACJob, 'mopacjob': MOPACJob,
    'reaxff': ReaxFFJob, 'reaxffjob': ReaxFFJob,
    'crs': CRSJob, 'cosmo-rs': CRSJob, 'crsjob': CRSJob
}


#: Schema for validating the ``['input_ligands']`` and ``['input_cores']`` blocks.
mol_schema: Schema = Schema({
    Optional_('guess_bonds', default=False):
        And(bool, error=".guess_bonds expects a boolean"),

    Optional_('is_core'):
        And(bool, error=".is_core expects a boolean"),

    Optional_('column'):
        And(int, lambda n: n >= 0, error=".column expects an integer larger than or equal to 0"),

    Optional_('row'):
        And(int, lambda n: n >= 0, error=".row expects an integer larger than or equal to 0"),

    Optional_('indices'):
        Or(
            And(
                int, lambda n: n >= 0, Use(lambda n: (n,)),
                error=".indices expects an integer larger than or equal to 0"
            ),
            And(
                abc.Collection,
                lambda n: all(isinstance(i, int) and i >= 0 for i in n),
                lambda n: len(n) == len(set(n)),
                Use(tuple),
                error=".indices expects one or more unique integers larger than or equal to 0"
            ),
            error=".indices expects an atomic index (int) or a list unique atomic indices"
        ),

    Optional_('type'):
        And(str, error='.type expects a string'),

    Optional_('name'):
        And(str, error='.name expects a string'),

    Optional_('path'):
        And(
            Use(validate_path),
            error=".path expects a None or a string pointing to an existing directory"
        ),
})

#: Schema for validating the ``['optional']['core']`` block.
core_schema: Schema = Schema({
    'dirname':
        And(str, error='optional.core.dirname expects a string'),

    Optional_('dummy', default=17):  # Return a tuple of atomic numbers
        Or(
            And(int, Use(to_atnum)),
            And(str, Use(to_atnum)),
            error='optional.core.dummy expects a valid atomic number (int) or symbol (string)'
        )
})

_db_names = ('core', 'ligand', 'qd')
_format_names = ('pdb', 'xyz')
_format_names2 = ('pdb', 'xyz', 'mol', 'mol2')

#: Schema for validating the ``['optional']['database']`` block.
database_schema: Schema = Schema({
    # path+directory name of the database
    'dirname':
        And(str, error='optional.database.dirname expects a string'),

    Optional_('read', default=_db_names):  # Attempt to pull structures from the database
        Or(
            And(bool, Use(lambda n: _db_names if n is True else ())),
            And(
                str,
                lambda n: n in _db_names,
                Use(lambda n: (n,)),
                error=f'allowed values for optional.database.read are: {repr(_db_names)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in _db_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.read are: {repr(_db_names)}'
            ),
            error='optional.database.read expects a boolean, string or list of unique strings'
        ),

    Optional_('write', default=_db_names):  # Attempt to write structures to the database
        Or(
            And(bool, Use(lambda n: _db_names if n is True else ())),
            And(
                str,
                lambda n: n in _db_names,
                Use(lambda n: (n,)),
                error=f'allowed values for optional.database.write are: {repr(_db_names)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in _db_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.write are: {repr(_db_names)}'
            ),
            error='optional.database.write expects a boolean, string or list of unique strings'
        ),

    Optional_('overwrite', default=tuple):  # Allow previous entries to be overwritten
        Or(
            And(bool, Use(lambda n: _db_names if n is True else ())),
            And(
                str,
                lambda n: n in _db_names,
                Use(lambda n: (n,)),
                error=f'allowed values for optional.database.overwrite are: {repr(_db_names)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in _db_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.overwrite are: {repr(_db_names)}'
            ),
            error='optional.database.overwrite expects a boolean, string or list of unique strings'
        ),

    Optional_('mongodb', default=dict):  # Settings specific to MongoDB
        Or(
            dict,
            And(bool, lambda n: n is False, Use(lambda n: {})),
            error='optional.database.mongodb expects True (boolean) or a dictionary'
        ),

    Optional_('mol_format', default=_format_names):  # Return a tuple of file formats
        Or(
            And(bool, Use(lambda n: _format_names if n is True else ())),
            And(
                str,
                lambda n: n in _format_names,
                error=f'allowed values for optional.database.mol_format are: {repr(_format_names2)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in _format_names for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.mol_format are: {repr(_format_names2)}'
            ),
            error='optional.database.mol_format expects a boolean, string or list of unique strings'
        )
})


#: Schema for validating the ``['optional']['ligand']`` block.
ligand_schema: Schema = Schema({
    # path+directory name of the ligand directory
    'dirname':
        And(str, error='optional.ligand.dirname expects a string'),

    Optional_('functional_groups', default=None):
        Or(
            None,
            And(str, Use(lambda n: (n,))),
            And(
                abc.Collection,
                lambda n: all(isinstance(i, str) for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error='optional.ligand.functional_groups expects a list of unique SMILES strings'
            ),
            error=('optional.ligand.functional_groups expects None (NoneType), a SMILES string, '
                   'or a list of unique SMILES string')
        ),

    Optional_('optimize', default=True):  # Optimize the ligands
        And(bool, error='optional.ligand.optimize expects a boolean'),

    Optional_('split', default=True):  # Remove a counterion from the function group
        And(bool, error='optional.ligand.split expects a boolean'),

    Optional_('cosmo-rs', default=False):  # Settings specific to ligand COSMO-RS calculations
        Or(
            dict,
            And(bool, Use(lambda n: {'job1': AMSJob} if n else False)),
            error='optional.ligand.cosmo-rs expects a boolean or dictionary'
        ),
})


#: Schema for validating the ``['optional']['qd']`` block.
qd_schema: Schema = Schema({
    # path+directory name of the quantum dot directory
    'dirname':
        And(str, error='optional.qd.dirname expects a string'),

    # Settings specific to a quantum dot activation strain analyses
    Optional_('activation_strain', default=False):
        And(bool, error='optional.qd.activation_strain expects a boolean'),

    Optional_('optimize', default=False):  # Settings for quantum dot geometry optimizations
        Or(
            dict,
            And(bool, Use(lambda n: ({'job1': AMSJob} if n else False))),
            error='optional.ligand.optimize expects a boolean or dictionary'
        ),

    # Settings for quantum dot ligand dissociation calculations
    Optional_('dissociate', default=False):
        Or(
            dict,
            And(bool, lambda n: n is False),
            error='optional.ligand.dissociate expects False (boolean) or a dictionary'
        )
})


#: Schema for validating the ``['optional']['database']['mongodb']`` block.
mongodb_schema: Schema = Schema({
    # Optional username for the MongoDB host
    Optional_('username'):
        Or(str, int, error='optional.database.mongodb.username expects a string or integer'),

    Optional_('password'):  # Optional password for the MongoDB host
        Or(str, int, error='optional.database.mongodb.password expects a string or integer'),

    Optional_('host', default='localhost'):  # Name of the MongoDB host
        Or(str, int, error='optional.database.mongodb.host expects a string or integer'),

    Optional_('port', default=27017):  # Port of the MongoDB host
        And(int, error='optional.database.mongodb.port expects an integer'),

    Optional_(str):  # Other keyword arguments for :class:`pymongo.MongoClient`
        object
})


#: Schema for validating the ``['optional']['qd']['dissociate']`` block.
bde_schema: Schema = Schema({
    # Atom type of the to-be dissociated core atom
    'core_atom':
        And(
            Or(int, str), Use(validate_core_atom),
            error=('optional.qd.dissociate.core_atom expects a SMILES string, '
                   'atomic number (int) or atomic symbol (str)')
        ),

    'lig_count':  # The number of ligands per core_atom
        And(
            int, lambda n: n >= 0,
            error='optional.qd.dissociate.lig_count expects an integer larger than or equal to 0'
        ),

    Optional_('keep_files', default=True):  # Delete files after the calculations are finished
        And(bool, error='optional.qd.dissociate.keep_files expects a boolean'),

    Optional_('core_core_dist', default=0.0):
        And(
            Or(int, float), lambda n: n >= 0.0, Use(float),
            error=('optional.qd.dissociate.core_core_dist expects an integer or float '
                   'larger than or equal to 0.0')
        ),

    Optional_('lig_core_dist', default=5.0):
        And(
            Or(int, float), lambda n: n >= 0.0, Use(float),
            error=('optional.qd.dissociate.lig_core_dist expects an integer or float '
                   'larger than or equal to 0.0')
        ),

    Optional_('core_index'):
        Or(
            And(
                int, lambda n: n >= 0, Use(lambda n: (n,)),
                error=('optional.qd.dissociate.core_index expects an integer '
                       'larger than or equal to 0')
            ),
            And(
                abc.Collection,
                lambda n: all(isinstance(i, int) and i >= 0 for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=('optional.qd.dissociate.core_index expects a list of unique integers '
                       'larger than or equal to 0')
            ),
            error=('optional.qd.dissociate.core_index expects an integer or list of unique integers'
                   'larger than or equal to 0')
        ),

    Optional_('topology', default=dict):
        And(
            dict, lambda n: all(isinstance(k, int) for k in n),
            error='optional.qd.dissociate.topology expects a dictionary with integers as keys'
        ),

    Optional_('job1', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.dissociate.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, lambda n: n.lower() in _class_dict, Use(str_to_job_type),
                error='optional.qd.dissociate.job1 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.dissociate.job1 expects a string or a type object'

        ),

    Optional_('s1', default=_bde_s1_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.dissociate.s1 expects a string or a dictionary'
        ),

    Optional_('job2'):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.dissociate.job2 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, lambda n: n.lower() in _class_dict, Use(str_to_job_type),
                error='optional.qd.dissociate.job2 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.dissociate.job2 expects a string or a type object'
        ),

    Optional_('s2'):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.dissociate.s2 expects a string or a dictionary'
        )
})

#: Schema for validating the ``['optional']['qd']['optimize']`` block.
qd_opt_schema: Schema = Schema({
    # The job type for the first half of the optimization
    Optional_('job1', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.opt.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, lambda n: n.lower() in _class_dict, Use(str_to_job_type),
                error='optional.qd.opt.job1 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.opt.job1 expects a string or a type object'
        ),

    # The job settings for the first half of the optimization
    Optional_('s1', default=_qd_opt_s1_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.opt.s1 expects a string or a dictionary'
        ),

    # The job type for the second half of the optimization
    Optional_('job2', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.opt.job2 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, lambda n: n.lower() in _class_dict, Use(str_to_job_type),
                error='optional.qd.opt.job2 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.opt.job2 expects a string or a type object'
        ),

    # The job settings for the second half of the optimization
    Optional_('s2', default=_qd_opt_s2_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.opt.s2 expects a string or a dictionary'
        )
})

#: Schema for validating the ``['optional']['ligand']['cosmo-rs']`` block.
crs_schema: Schema = Schema({
    # Delete files after the calculations are finished
    Optional_('keep_files', default=True):
        And(bool, error='optional.ligand.cosmo-rs.keep_files expects a boolean'),

    # The job type for constructing the COSMO surface
    Optional_('job1', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.ligand.cosmo-rs.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, lambda n: n.lower() in _class_dict, Use(str_to_job_type),
                error=('optional.ligand.cosmo-rs.job1 expects a string '
                       'that is a valid plams.Job alias')
            ),
            error='optional.ligand.cosmo-rs.job1 expects a string or a type object'
        ),

    # The settings for constructing the COSMO surface
    Optional_('s1', default=_crs_s1_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.ligand.cosmo-rs.s1 expects a string or a dictionary'
        ),

    Optional_('job2', default=_get_crsjob):  # The job type for the actual COSMO-RS calculation
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.ligand.cosmo-rs.job2 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, lambda n: n.lower() in _class_dict, Use(str_to_job_type),
                error=('optional.ligand.cosmo-rs.job2 expects a string '
                       'that is a valid plams.Job alias')
            ),
            error='optional.ligand.cosmo-rs.job2 expects a string or a type object'
        ),

    # The settings for the actual COSMO-RS calculation
    Optional_('s2', default=_crs_s2_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.ligand.cosmo-rs.s2 expects a string or a dictionary'
        )
})
