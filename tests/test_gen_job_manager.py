"""Tests for the :mod:`CAT.gen_job_manager.GenJobManager` class in :mod:`CAT.gen_job_manager`."""

from os.path import (join, abspath)
from collections import abc

from scm.plams import Settings, JobManager, AMSJob, AMSResults
from assertionlib import assertion

from CAT.gen_job_manager import GenJobManager

SETTINGS = Settings({'counter_len': 3, 'hashing': 'input', 'remove_empty_directories': True})
PATH = join('tests', 'test_files')
ABSPATH = abspath(PATH)
FOLDER = 'test_plams_workdir'


def test_init() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager.__init__`."""
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    assertion.isinstance(manager, GenJobManager)
    assertion.issubclass(manager.__class__, JobManager)
    assertion.eq(manager.settings, SETTINGS)
    assertion.eq(manager.jobs, [])
    assertion.eq(manager.names, {})
    assertion.eq(manager.hashes, {})
    assertion.eq(manager.path, ABSPATH)
    assertion.eq(manager.foldername, FOLDER)
    assertion.eq(manager.workdir, join(ABSPATH, FOLDER))
    assertion.eq(manager.logfile, join(ABSPATH, FOLDER, 'logfile'))
    assertion.eq(manager.input, join(ABSPATH, FOLDER, SETTINGS.hashing))


def test_load_job() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager.load_job`."""
    filename = join(PATH, FOLDER, 'QD_opt_part1', 'QD_opt_part1.hash')
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    manager.load_job(filename)

    k, v = next(iter(manager.hashes.items()))
    assertion.eq(k, '0da9b13507022986d26bbc57b4c366cf1ead1fe70ff750e071e79e393b14dfb5')
    assertion.isinstance(v, abc.Callable)
    assertion.eq(v.__name__, 'unpickle_job')

    job = v()
    assertion.isinstance(job, AMSJob)
    assertion.eq(job.status, 'successful')
    assertion.isinstance(job.results, AMSResults)
    assertion.eq(job.name, 'QD_opt_part1')
    assertion.eq(job.path, '/Users/basvanbeek/Documents/CdSe/Week_5/qd/QD_optimize/QD_opt_part1')
    assertion.isinstance(job.settings, Settings)
    assertion.eq(job.depend, [])
    assertion.eq(job._dont_pickle, [])
    assertion.eq(job.molecule.get_formula(), 'C78Cd68H182O26Se55')


def _test_check_hash() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager._check_hash`."""
    filename = join(PATH, FOLDER, 'QD_opt_part1', 'QD_opt_part1.hash')
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    manager.load_job(filename)

    v = next(iter(manager.hashes.values()))
    j1 = v()
    j2 = manager._check_hash(j1)

    assertion.ne(j2, j1)
    for k in vars(j1):
        if k in ('results', 'molecule', '_hash'):
            continue
        attr1, attr2 = getattr(j1, k), getattr(j2, k)
        assertion.eq(attr1, attr2)


def _test_remove_job() -> None:
    """Test :meth:`CAT.gen_job_manager.GenJobManager.test_remove_job`."""
    filename = join(PATH, FOLDER, 'QD_opt_part1', 'QD_opt_part1.hash')
    manager = GenJobManager(SETTINGS, PATH, FOLDER)
    manager.load_job(filename)

    k, v = next(iter(manager.hashes.items()))
    job = v()

    assertion.contains(manager.hashes, k)
    manager.remove_job(job)
    assertion.eq(manager.hashes, {})
