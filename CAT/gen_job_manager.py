"""
CAT.gen_job_manager
===================

A module to hold the :class:`GenJobManager` class.

Index
-----
.. currentmodule:: CAT.gen_job_manager
.. autosummary::
    GenJobManager

API
---
.. autoclass:: GenJobManager
   :members:
   :private-members:
   :special-members:

"""

import threading
from typing import (Optional, Sequence, Callable)
from os.path import (join, isfile, abspath, dirname)

try:
    import dill as pickle
except ModuleNotFoundError:
    import pickle

from scm.plams import (JobManager, log, config, FileError)
from scm.plams.core.basejob import (Job, MultiJob)


__all__ = ['GenJobManager']


class GenJobManager(JobManager):
    """A modified version of the PLAMS :class:`JobManager` class.

    In order to minimize its memory footprint, :meth:`GenJobManager.load_job` will now populate
    the values of :attr:`GenJobManager.hashes` with a callable the creates and returns a
    :class:`Job` instance (rather than permantly sotring the :class:`Job` instance).

    Attributes
    ----------
    foldername : str
        The working folder name.

    workdir : str
        The absolute path to the working folder.

    logfile : str
        The absolute path to the logfile.

    input : str
        The absolute path to the copy of the input file in the working folder.

    settings : |plams.Settings|_
        A :class:`Settings` instance for this job manager.

    jobs : |list|_ [|plams.Jobs|_]
        A list of all jobs managed with this instance (in order of run calls).

    names : |dict|_ [|str|_, |int|_]
        A dictionary with names of jobs.
        For each name an integer value is stored indicating
        how many jobs with that basename have already been run.

    hashes : |dict|_ [|str|_, |Callable|_]
        A dictionary working as a hash-table for jobs.
        Values created by :meth:`GenJobManager.load_job` are stored as callables which create
        :class:`Job` instances.

    """

    def __init__(self, *args: Sequence, **kwargs: dict) -> None:
        """Initialize the :class:`GenJobManager` instance."""
        super().__init__(*args, **kwargs)

    @staticmethod
    def _unpickle(filename: str) -> Optional[Job]:
        """Attempt to unpickle and return a :class:`Job` containing file."""
        with open(filename, 'rb') as f:
            try:
                return pickle.load(f)
            except Exception as e:
                log(f"Unpickling of {filename} failed. Caught the following Exception:\n{e}", 1)
                return None

    @staticmethod
    def _get_job(filename: str) -> Callable:
        """Return a callable which converts **filename** into a :class:`Job` instance."""
        def unpickle_job() -> Optional[Job]:
            return GenJobManager._unpickle(filename)
        return unpickle_job

    def load_job(self, filename: str) -> None:
        """Load previously saved job from **filename**.

        Parameters
        ----------
        filename : str
            A path to a .dill file in some job folder.
            A |Job| instance stored there is loaded and returned.
            All attributes of this instance removed before pickling are restored.

        """
        def setstate(job, path, parent=None):
            job.parent = parent
            job.jobmanager = self
            job.default_settings = [config.job]
            job.path = path
            if isinstance(job, MultiJob):
                job._lock = threading.Lock()
                for child in job:
                    setstate(child, join(path, child.name), job)
                for otherjob in job.other_jobs():
                    setstate(otherjob, join(path, otherjob.name), job)

            job.results.refresh()
            h = job.hash()
            if h is not None:
                filename = join(path, job.name + '.dill')
                self.hashes[h] = self._get_job(filename)

        # Raise an error if **filename** cannot be found
        if isfile(filename):
            filename = abspath(filename)
        else:
            raise FileError('File {} not present'.format(filename))

        path = dirname(filename)
        job = self._unpickle(filename)
        setstate(job, path)

    def remove_job(self, job: Job) -> None:
        """Remove **job** from the job manager; forget its hash."""
        if job in self.jobs:
            self.jobs.remove(job)
            job.jobmanager = None
        h = job.hash()
        if h in self.hashes and self.hashes[h]() == job:
            del self.hashes[h]

    def _check_hash(self, job: Job) -> Optional[Job]:
        """Calculate and check the hash of **job**.

        If the hash is not ``None``, previous jobs are searched for the same hash.
        Matches are returned if applicable; returns ``None`` otherwise.

        Parameter
        ---------
        job : |plams.Job|_
            A :class:`Job` instance.

        Returns
        -------
        |plams.Job|_
            Optional: A matching :class:`Job` instance stored in :attr:`GenJobManager.hashes`.
            If no matching hash is found return ``None``.

        """
        h = job.hash()
        if h is None:  # No hash available for **job**, move along
            return None

        if h in self.hashes:
            prev = self.hashes[h]()
            log(f'Job {job.name} previously run as {prev.name}, using old results', 1)
            return prev
        else:
            filename = join(job.path, job.name + '.dill')
            func = self._get_job(filename)
            self.hashes[h] = func  # Set a callable that returns a Job instance
            return None
