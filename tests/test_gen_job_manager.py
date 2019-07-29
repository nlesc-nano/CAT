"""Tests for the :mod:`CAT.gen_job_manager.GenJobManager` class in :mod:`CAT.gen_job_manager`."""

from os.path import (join, abspath)
from collections import abc

from scm.plams import Settings, JobManager, AMSJob, AMSResults

from CAT.gen_job_manager import GenJobManager
from CAT.assertion_functions import (
    assert_eq, assert_instance, assert_subclass, assert_isin, Invert
)

SETTINGS = Settings({'counter_len': 3, 'hashing': 'input', 'remove_empty_directories': True})
PATH = join('tests', 'test_files')
ABSPATH = abspath(PATH)
FOLDER = 'test_plams_workdir'


def test_init() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager.__init__`."""
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    assert_instance(manager, GenJobManager)
    assert_subclass(manager.__class__, JobManager)
    assert_eq(manager.settings, SETTINGS)
    assert_eq(manager.jobs, [])
    assert_eq(manager.names, {})
    assert_eq(manager.hashes, {})
    assert_eq(manager.path, ABSPATH)
    assert_eq(manager.foldername, FOLDER)
    assert_eq(manager.workdir, join(ABSPATH, FOLDER))
    assert_eq(manager.logfile, join(ABSPATH, FOLDER, 'logfile'))
    assert_eq(manager.input, join(ABSPATH, FOLDER, SETTINGS.hashing))


def test_load_job() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager.load_job`."""
    filename = join(PATH, FOLDER, 'QD_opt_part1', 'QD_opt_part1.hash')
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    manager.load_job(filename)

    k, v = next(iter(manager.hashes.items()))
    assert_eq(k, '0da9b13507022986d26bbc57b4c366cf1ead1fe70ff750e071e79e393b14dfb5')
    assert_instance(v, abc.Callable)
    assert_eq(v.__name__, 'unpickle_job')

    job = v()
    assert_instance(job, AMSJob)
    assert_eq(job.status, 'successful')
    assert_instance(job.results, AMSResults)
    assert_eq(job.name, 'QD_opt_part1')
    assert_eq(job.path, '/Users/basvanbeek/Documents/CdSe/Week_5/qd/QD_optimize/QD_opt_part1')
    assert_instance(job.settings, Settings)
    assert_eq(job.depend, [])
    assert_eq(job._dont_pickle, [])
    assert_eq(job.molecule.get_formula(), 'C78Cd68H182O26Se55')


def test_check_hash() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager._check_hash`."""
    filename = join(PATH, FOLDER, 'QD_opt_part1', 'QD_opt_part1.hash')
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    manager.load_job(filename)

    v = next(iter(manager.hashes.values()))
    j1 = v()
    j2 = manager._check_hash(j1)

    with Invert(assert_eq) as func:
        func(j2, j1)
    for k in vars(j1):
        if k in ('results', 'molecule', '_hash'):
            continue
        attr1, attr2 = getattr(j1, k), getattr(j2, k)
        assert_eq(attr1, attr2)


def test_remove_job() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager.test_remove_job`."""
    filename = join(PATH, FOLDER, 'QD_opt_part1', 'QD_opt_part1.hash')
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    manager.load_job(filename)

    k, v = next(iter(manager.hashes.items()))
    job = v()

    assert_isin(k, manager.hashes)
    manager.remove_job(job)
    assert_eq(manager.hashes, {})
