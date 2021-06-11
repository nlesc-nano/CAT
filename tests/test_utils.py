"""Tests for :mod:`CAT.utils`."""

import os
from os.path import join

from unittest import mock

from scm.plams import config
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob
from assertionlib import assertion

from CAT.utils import (
    type_to_string, dict_concatenate, get_template, validate_path, check_sys_var, restart_init
)

PATH = join('tests', 'test_files')
FOLDER = 'test_plams_workdir'


def test_type_to_string() -> None:
    """Test :func:`CAT.utils.type_to_string`."""
    ref = {
        ADFJob: 'adf',
        AMSJob: 'ams',
        DiracJob: 'dirac',
        Cp2kJob: 'cp2k',
        GamessJob: 'gamess',
        ORCAJob: 'orca'
    }

    for k, j in ref.items():
        i = type_to_string(k)
        assertion.eq(i, j)
    assertion.eq(type_to_string('bob'), '')


def test_dict_concatenate() -> None:
    """Test :func:`CAT.utils.dict_concatenate`."""
    ref = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    dict_list = [{'a': 1}, {'b': 2}, {'c': 3}, {'d': 4}]
    out = dict_concatenate(dict_list)
    assertion.eq(out, ref)


def test_get_template() -> None:
    """Test :func:`CAT.utils.get_template`."""
    ref = {'CRSParameters': {'_1': 'HB_HNOF', '_2': 'HB_TEMP', '_3': 'FAST', '_4': 'COMBI2005', 'rav': 0.4, 'aprime': 1550.0, 'fcorr': 2.802, 'chb': 0.0, 'sigmahbond': 0.00978, 'aeff': 5.96, 'Lambda': 0.135, 'omega': -0.212, 'Eta': -9.65, 'chortf': 0.816}, 'Dispersion': {'H': -0.034, 'C': -0.0356, 'N': -0.0224, 'O': -0.0333, 'F': -0.026, 'Si': -0.04, 'P': -0.045, 'S': -0.052, 'Cl': -0.0485, 'Br': -0.055, 'I': -0.062}, 'Technical': {'rsconv': '1e-7', 'maxiter': 10000, 'bpconv': '1e-6', 'bpmaxiter': 40, 'solconv': '1e-5', 'solmaxiter': 40, 'solxilarge': 0.99, 'ehdeltaT': 1.0}}  # noqa
    out = get_template('crs.yaml')['MOPAC PM6']['input'].as_dict()
    assertion.eq(out, ref)


def test_validate_path() -> None:
    """Test :func:`CAT.utils.validate_path`."""
    assertion.eq(validate_path(None), os.getcwd())
    assertion.eq(validate_path(''), os.getcwd())
    assertion.eq(validate_path('.'), os.getcwd())
    assertion.eq(validate_path(PATH), PATH)
    assertion.assert_(validate_path, join(PATH, 'bob'), exception=FileNotFoundError)
    assertion.assert_(validate_path, join(PATH, 'Methanol.xyz'), exception=NotADirectoryError)


def test_check_sys_var() -> None:
    """Test :func:`CAT.utils.validate_path`."""
    @mock.patch.dict(os.environ,
                     {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
    def test1() -> None:
        assertion.is_(check_sys_var(), None)

    @mock.patch.dict(os.environ,
                     {'ADFBIN': '', 'ADFHOME': '2019', 'ADFRESOURCES': '', 'SCMLICENSE': ''})
    def test2() -> None:
        assertion.assert_(check_sys_var, exception=EnvironmentError)

    test1()
    test2()


def test_restart_init() -> None:
    """Test :func:`CAT.utils.restart_init`."""
    restart_init(PATH, FOLDER)
    manager = config.default_jobmanager

    _hash = '0da9b13507022986d26bbc57b4c366cf1ead1fe70ff750e071e79e393b14dfb5'
    assertion.contains(manager.hashes, _hash)
