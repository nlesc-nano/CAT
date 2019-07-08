"""Tests for :mod:`CAT.utils`."""

from os.path import join

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob

from CAT.assertion_functions import (assert_value, assert_lt)
from CAT.utils import (type_to_string, get_time, dict_concatenate, get_template)

PATH = 'test/test_files'


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
        assert_value(i, j)
    assert_value(type_to_string('bob'), '')


def test_get_time() -> None:
    """Test :func:`CAT.utils.get_time`."""
    time = get_time()
    assert_value(time[0], '[')
    assert_value(time[9:], '] ')
    assert_value(time[3], ':')
    assert_value(time[6], ':')
    assert_lt(int(time[1:3]), 25)
    assert_lt(int(time[4:6]), 61)
    assert_lt(int(time[7:9]), 61)


def test_dict_concatenate() -> None:
    """Test :func:`CAT.utils.dict_concatenate`."""
    ref = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    dict_list = [{'a': 1}, {'b': 2}, {'c': 3}, {'d': 4}]
    out = dict_concatenate(dict_list)
    assert_value(out, ref)


def test_get_template() -> None:
    """Test :func:`CAT.utils.get_template`."""
    ref = {'CRSParameters': {'_1': 'HB_HNOF', '_2': 'HB_TEMP', '_3': 'FAST', '_4': 'COMBI2005', 'rav': 0.4, 'aprime': 1550.0, 'fcorr': 2.802, 'chb': 0.0, 'sigmahbond': 0.00978, 'aeff': 5.96, 'Lambda': 0.135, 'omega': -0.212, 'Eta': -9.65, 'chortf': 0.816}, 'Dispersion': {'H': -0.034, 'C': -0.0356, 'N': -0.0224, 'O': -0.0333, 'F': -0.026, 'Si': -0.04, 'P': -0.045, 'S': -0.052, 'Cl': -0.0485, 'Br': -0.055, 'I': -0.062}, 'Technical': {'rsconv': '1e-7', 'maxiter': 10000, 'bpconv': '1e-6', 'bpmaxiter': 40, 'solconv': '1e-5', 'solmaxiter': 40, 'solxilarge': 0.99, 'ehdeltaT': 1.0}}  # noqa
    out = get_template('crs.yaml')['MOPAC PM6']['input'].as_dict()
    assert_value(out, ref)
