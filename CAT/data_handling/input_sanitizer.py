""" A module designed for sanitizing and interpreting the input file. """

__all__ = ['sanitize_optional', 'sanitize_input_mol', 'sanitize_path']

import os
from os.path import (join, isdir, isfile, exists)
from itertools import chain

import yaml
from schema import (Schema, Or, And, Use)

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

from scm.plams.mol.molecule import Molecule
from scm.plams.core.basejob import Job
from scm.plams.core.settings import Settings
from scm.plams.tools.periodic_table import PeriodicTable
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from .. import utils as CAT

from ..utils import get_time
from ..mol_utils import to_atnum
from ..analysis.crs import CRSJob


""" ###################################  Sanitize path  ####################################### """


def sanitize_path(arg):
    """ Sanitize and return the settings of arg.path. """
    if arg.path is None:
        arg.path = os.getcwd()
        return arg
    elif isinstance(arg.path, str):
        if arg.path.lower() in ('none', '.', 'pwd', '$pwd', 'cwd'):
            arg.path = os.getcwd()
        else:
            if not os.path.exists(arg.path):
                raise FileNotFoundError(get_time() + 'path ' + arg.path + ' not found')
            elif os.path.isfile(arg.path):
                raise TypeError(get_time() + 'path ' + arg.path + ' is a file, not a directory')
        return arg

    else:
        error = 'arg.path should be None or a string, ' + str(type(arg.path))
        error += ' is not a valid type'
        raise TypeError(error)


""" ##########################  Sanitize input_ligands & input_cores  ######################## """


def sanitize_input_mol(arg):
    """ Sanitize and return the settings of arg.input_cores & arg.input_ligands. """
    core_path = arg.optional.core.dirname
    arg.input_cores = get_mol_defaults(arg.input_cores, path=core_path, core=True)
    arg.input_cores = sanitize_mol_type(arg.input_cores)

    ligand_path = arg.optional.ligand.dirname
    arg.input_ligands = get_mol_defaults(arg.input_ligands, path=ligand_path, core=False)
    arg.input_ligands = sanitize_mol_type(arg.input_ligands)

    return arg


def get_mol_defaults(mol_list, path=None, core=False):
    """ Prepare the default input settings for a molecule. """
    key_dict = {
        'guess_bonds': val_bool,
        'is_core': val_bool,
        'column': val_int,
        'row': val_int,
        'indices': val_indices,
        'type': val_type,
        'name': val_string,
        'path': val_string,
    }

    ret = []
    for mol in mol_list:
        tmp = get_default_input_mol()
        tmp.mol = mol
        tmp.path = path
        tmp.is_core = core

        if isinstance(mol, dict):
            for key1 in mol:
                tmp.mol = key1
                for key2 in mol[key1]:
                    try:
                        tmp[key2] = key_dict[key2](mol[key1][key2])
                    except KeyError:
                        raise KeyError(str(key2) + ' is not a valid argument for ' + str(key1))
                    if key2 == 'guess_bonds':
                        tmp.tmp_guess = True

        ret.append(tmp)
    return ret


def sanitize_mol_type(input_mol):
    """ Sanitize and return the (file) type of the input molecule (SMILES, .xyz, dir, etc...). """
    for mol in input_mol:
        # Figure out the (file) type and mol name
        try:
            if isfile(join(mol.path, mol.mol)):
                mol.type = mol.mol.rsplit('.', 1)[-1]
                mol.name = mol.mol.rsplit('.', 1)[0]
                mol.mol = join(mol.path, mol.mol)
                if mol.type == '.xyz' and not mol.get('tmp_guess'):
                    mol.guess_bonds = False
            elif isdir(join(mol.path, mol.mol)):
                mol.type = 'folder'
                mol.name = mol.mol
                mol.mol = join(mol.path, mol.mol)
            elif isfile(mol.mol):
                mol.type = mol.mol.rsplit('.', 1)[-1]
                mol.name = mol.mol.rsplit('.', 1)[0].rsplit('/', 1)[-1].rsplit('\\', 1)[-1]
            elif isdir(mol.mol):
                mol.type = 'folder'
                mol.name = mol.mol
            else:
                mol.type = 'smiles'
                mol.name = santize_smiles(mol.mol)
        except TypeError:
            if isinstance(mol.mol, Molecule):
                mol.type = 'plams_mol'
                if not mol.properties.name:
                    mol.name = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol.mol)))
                    mol.name = Chem.CanonSmiles(mol.name)
                else:
                    mol.name = mol.properties.name
            elif isinstance(mol.mol, Chem.rdchem.Mol):
                mol.type = 'rdmol'
                mol.name = Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol.mol)))

    return input_mol


def get_default_input_mol():
    """ Return the default settings of arg.input_cores & arg.input_ligands. """
    ret = yaml.load("""
        mol: None
        name: None
        path: None
        guess_bonds: False
        is_core: False
        column: 0
        row: 0
        indices: None
        type: None
    """)

    for key in ret:
        if ret[key] == 'None':
            ret[key] = None

    return Settings(ret)


def santize_smiles(string):
    """ Sanitize a SMILES string: turn it into a valid filename. """
    name = string.replace('(', '[').replace(')', ']')
    cis_trans = [item for item in string if item == '/' or item == '\\']
    if cis_trans:
        cis_trans = [item + cis_trans[i*2+1] for i, item in enumerate(cis_trans[::2])]
        cis_trans_dict = {'//': 'trans-', '/\\': 'cis-'}
        for item in cis_trans[::-1]:
            name = cis_trans_dict[item] + name
        name = name.replace('/', '').replace('\\', '')

    return name


""" ####################################  Sanitize optional  ################################## """


def sanitize_optional(arg_dict):
    """ Sanitize and return the settings of arg.optional. """
    arg = get_default_optional()
    arg.update(arg_dict)

    mol_format = ('xyz', 'pdb')

    # Validate arguments consisting of booleans, integers, strings and/or iterables
    arg.optional.core.dirname = val_dir_names(arg.optional.core.dirname, arg.path)
    arg.optional.core.dummy = val_atnum(arg.optional.core.dummy)
    arg.optional.database.dirname = val_dir_names(arg.optional.database.dirname, arg.path)
    arg.optional.database.read = val_data(arg.optional.database.read)
    arg.optional.database.write = val_data(arg.optional.database.write)
    arg.optional.database.overwrite = val_data(arg.optional.database.overwrite)
    arg.optional.database.mol_format = val_format(arg.optional.database.mol_format, mol_format)
    arg.optional.database.mongodb = False
    arg.optional.ligand.dirname = val_dir_names(arg.optional.ligand.dirname, arg.path)
    arg.optional.ligand.optimize = val_bool(arg.optional.ligand.optimize)
    arg.optional.ligand.split = val_bool(arg.optional.ligand.split)
    arg.optional.qd.dirname = val_dir_names(arg.optional.qd.dirname, arg.path)
    arg.optional.qd.activation_strain = val_bool(arg.optional.qd.activation_strain)

    # Prepares COSMO-RS default settings
    crs = CAT.get_template('qd.yaml')['COSMO-RS activity coefficient']
    crs.update(CAT.get_template('crs.yaml')['MOPAC PM6'])

    # Validate arguments containing job recipes
    arg.optional.ligand.crs = val_job(arg.optional.ligand['cosmo-rs'],
                                      job1=AMSJob,
                                      job2=CRSJob,
                                      s1=CAT.get_template('qd.yaml')['COSMO-MOPAC'],
                                      s2=crs)
    del arg.optional.ligand['cosmo-rs']

    arg.optional.qd.optimize = val_job(arg.optional.qd.optimize,
                                       job1=AMSJob,
                                       job2=AMSJob,
                                       s1=CAT.get_template('qd.yaml')['UFF'],
                                       s2=CAT.get_template('qd.yaml')['UFF'])

    arg.optional.qd.dissociate = val_dissociate(arg.optional.qd.dissociate)

    del arg.path
    return arg


def get_default_optional():
    """ Return the default settings of arg.optional. """
    ret = yaml.load("""
        optional:
            database:
                dirname: database
                read: True
                write: True
                overwrite: False
                mol_format: [pdb, xyz]
                mongodb: False

            core:
                dirname: core
                dummy: Cl

            ligand:
                dirname: ligand
                optimize: True
                cosmo-rs: False
                split: True

            qd:
                dirname: QD
                optimize: False
                activation_strain: False
                dissociate: False
    """)

    return Settings(ret)


def get_default_dissociate():
    """ Return the default settings of arg.optional. """
    ret = yaml.load("""
        core_atom: Cd
        lig_count: 2
        core_core_dist: 5.0
        lig_core_dist: 5.0
        topology:
            7: vertice
            8: edge
            10: face

        job1: AMSJob
        s1: True
        job2: AMSJob
        s2: True
    """)

    return Settings(ret)


str_to_class = {
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


def val_format(arg, ref):
    """ Validate database.mol_format & database_format. """
    schema = Schema(Or(
            And(None, Use(bool)),
            And(bool, lambda n: n is False),
            And(str, lambda n: not n, Use(bool)),
            And(str, lambda n: n.lower().rsplit('.', 1)[-1] in ref),
            And([str], lambda n: [i.lower().rsplit('.', 1)[-1] in ref for i in n], Use(list))
    ))

    # Decapitalize and remove any periods.
    ret = schema.validate(arg)
    if isinstance(ret, list):
        for i, item in enumerate(ret):
            ret[i] = item.lower().rsplit('.', 1)[-1]
        ret = tuple(ret)
    elif isinstance(ret, str):
        ret = (ret.lower().rsplit('.', 1)[-1])
    elif not ret:
        ret = ()

    return ret


def val_data(arg):
    """ Validate the input arguments for database.read, write and overwrite.
    Returns *False* or tuple with *ligand*, *core* and/or *qd*.
    """
    ref = ('ligand', 'core', 'qd')

    def get_arg(n):
        if n:
            return ref
        else:
            return False

    def get_false(n):
        return False

    schema = Schema(Or(
            And(bool, Use(get_arg)),
            And(str, lambda n: not n, Use(bool)),
            And(str, lambda n: n.lower() in ref, Use(list)),
            And([str], lambda n: not any([bool(i) for i in n]), Use(get_false)),
            And([str], lambda n: [i.lower() in ref for i in n], Use(list))
    ))

    # Decapitalize
    ret = schema.validate(arg)
    if isinstance(ret, list):
        for i, item in enumerate(ret):
            ret[i] = item.lower()
        ret = tuple(ret)
    elif not ret:
        ret = ()

    return ret


def val_type(file_type):
    """ Validate a the fle type, returns a <str> or <None>. """
    return Schema(Or(str, None, Molecule, Chem.rdchem.Mol)).validate(file_type)


def val_int(integer):
    """ Validate a positive integer; returns an <int>. """
    schema = Schema(And([int], lambda n: n >= 0))
    return schema.validate(integer)


def val_string(string):
    """ Validate a string; returns a <str>. """
    return Schema(str).validate(string)


def val_indices(indices):
    """ Validate an iterable consisting if integers; returns a <tuple> consisting of 3 <str>. """
    if indices is None:
        return tuple()
    schema = Schema(And([int], Use(tuple), lambda n: [i >= 0 for i in n]))
    return schema.validate(list(indices))


def val_dir_names(dirname, path):
    """ Validate a str; returns a str.
    Creates a directory at path/dirname if it does not yet exist. """
    ret = join(path, Schema(str).validate(dirname))
    if not exists(ret):
        os.makedirs(ret)
    else:
        assert isdir(ret)
    return ret


def val_atnum(atnum):
    """ Validate an atomic number or symbol; returns an atomic number <int>. """
    at_gen = chain.from_iterable([[i, j[0]] for i, j in enumerate(PeriodicTable.data)])
    schema = Schema(And(Or(int, str), lambda n: n in at_gen, Use(to_atnum)))
    return schema.validate(atnum)


def val_bool(my_bool):
    """ Validate a boolean; returns a <bool>. """
    return Schema(bool).validate(my_bool)


def val_job(job, job1=None, job2=None, s1=None, s2=None):
    """ Validate a job recipe.
    Returns a dictionary: {'job1': <Job>, 'job2': <Job>, 's1': <Settings>, 's2': <Settings>}. """
    # Validate the object type
    Schema(Or(bool, dict)).is_valid(job)
    if isinstance(job, bool):
        if job is False:
            return job
        job = {'job1': None, 's1': None,
               'job2': None, 's2': None}

    # Validate the object types of the various elements
    schema = Schema({'job1': Or(None, Job, str),
                     's1': Or(None, dict),
                     'job2': Or(None, Job, str),
                     's2': Or(None, dict)})
    schema.is_valid(job)

    # Assign proper default settings
    str_to_def = {'job1': job1, 'job2': job2, 's1': s1, 's2': s2}
    for key in job:
        if job[key] is None or 'None':
            job[key] = str_to_def[key]
        if not job[key]:
            job[key] = False
        elif isinstance(job[key], str):
            try:
                job[key] = str_to_class[job[key].lower()]
            except KeyError:
                raise KeyError(get_time() + 'No Job-derived object exists for the string:', job[key]
                               + ', please provide the actual <Job> object instead of <str>')
        elif isinstance(job[key], type):
            pass
        elif isinstance(job[key], dict):
            job[key] = Settings(job[key])
        else:
            raise TypeError(get_time() + str(type(job[key])), 'is an unspported object type')
    return job


def val_dissociate(dissociate):
    """ Validate the optional.qd.dissociate block in the input file. """
    ret = get_default_dissociate()
    if dissociate is True:
        dissociate = Settings()
    elif dissociate is False:
        return False

    ret.update(dissociate)
    if ret.job1 is False or ret.s1 is False:
        return False

    ret.core_atom = to_atnum(ret.core_atom)
    ret.lig_count = int(ret.lig_count)
    ret.core_core_dist = float(ret.core_core_dist)
    ret.lig_core_dist = float(ret.lig_core_dist)
    assert isinstance(ret.topology, dict)
    for key in ret.topology:
        assert isinstance(key, int)
        assert isinstance(ret.topology[key], str)

    # Interpret job1
    assert isinstance(ret.job1, (bool, type, str))
    if ret.job1 is True:
        ret.job1 = AMSJob
    elif isinstance(ret.job1, str):
        ret.job1 = str_to_class[ret.job1.lower()]

    # Interpret job2
    assert isinstance(ret.job2, (bool, type, str))
    if ret.job2 is True:
        ret.job2 = AMSJob
    elif ret.job2 is False:
        ret.s2 = False
    elif isinstance(ret.job2, str):
        ret.job2 = str_to_class[ret.job2.lower()]

    # Interpret s1
    assert isinstance(ret.s1, (bool, dict, str))
    if ret.s1 is True:
        ret.s1 = CAT.get_template('qd.yaml')['MOPAC']
    elif isinstance(ret.s1, str):
        if isfile(ret.s1):
            ret.s1 = CAT.get_template(ret.s1, from_cat_data=False)
        else:
            raise FileNotFoundError(get_time() + str(ret.s1) + ' was not found')

    # Interpret s2
    assert isinstance(ret.s2, (bool, dict, str))
    if ret.s2 is True:
        ret.s2 = CAT.get_template('qd.yaml')['UFF']
    elif ret.s2 is False:
        ret.job2 = False
    elif isinstance(ret.s2, str):
        if isfile(ret.s2):
            ret.s2 = CAT.get_template(ret.s2, from_cat_dataj=False)
        else:
            raise FileNotFoundError(get_time() + str(ret.s1) + ' was not found')

    return ret
