"""A module with miscellaneous functions.

Index
-----
.. currentmodule:: CAT.utils
.. autosummary::
    type_to_string
    get_time
    check_sys_var
    dict_concatenate
    get_template
    KindEnum
    AnchorTup

API
---
.. autofunction:: type_to_string
.. autofunction:: get_time
.. autofunction:: check_sys_var
.. autofunction:: dict_concatenate
.. autofunction:: get_template
.. autoclass:: KindEnum
.. autoclass:: AnchorTup

"""

import os
import enum
import yaml
import inspect
import textwrap
import operator
import threading
import pkg_resources as pkg
from types import MappingProxyType, TracebackType
from shutil import rmtree
from pprint import pformat
from logging import Logger
from os.path import join, isdir, isfile, exists
from itertools import cycle, chain, repeat
from contextlib import redirect_stdout
from collections import abc, Counter
from typing import (
    Iterable, Union, TypeVar, Mapping, Type, Generator, Iterator, Optional,
    Any, NoReturn, Dict, overload, Callable, NamedTuple, Tuple,
)

import numpy as np
from scipy.spatial import cKDTree

from nanoutils import as_nd_array
from rdkit.Chem import Mol
from scm.plams.core.basejob import Job
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams import (
    config, Settings, Molecule, MoleculeError, PeriodicTable, init,
    from_smiles, AMSJob, ADFJob, Cp2kJob, DiracJob, GamessJob, CRSJob
)

from .logger import logger
from .mol_utils import to_atnum
from .gen_job_manager import GenJobManager
from ._mol_str_parser import FormatEnum

__all__ = [
    'JOB_MAP', 'check_sys_var', 'dict_concatenate', 'get_template',
    'cycle_accumulate', 'iter_repeat', 'get_nearest_neighbors',
    'log_traceback_locals', 'KindEnum', 'AnchorTup', 'AllignmentEnum',
    'AllignmentTup', 'FormatEnum',
]

JOB_MAP: Mapping[Type[Job], str] = MappingProxyType({
    ADFJob: 'adf',
    AMSJob: 'ams',
    DiracJob: 'dirac',
    Cp2kJob: 'cp2k',
    GamessJob: 'gamess',
    ORCAJob: 'orca',
    CRSJob: 'cosmo-rs'
})

T1 = TypeVar('T1')
T2 = TypeVar('T2')
KT = TypeVar('KT')
VT = TypeVar('VT')
P = TypeVar('P', str, bytes, os.PathLike)
Dtype = Union[type, str, np.dtype]


def type_to_string(job: Type[Job]) -> str:
    """Turn a :class:`type` instance into a :class:`str`."""
    try:
        return JOB_MAP[job]
    except KeyError:
        logger.error(f"No default settings available for type: '{job.__class__.__name__}'")
        return ''


def check_sys_var() -> None:
    """Validate all ADF environment variables.

    Raises
    ------
    EnvironmentError
        Raised if one or more of the following environment variables are absent:
        * ``'AMSBIN'``
        * ``'AMSHOME'``
        * ``'AMSRESOURCES'``
        * ``'SCMLICENSE'``

    """
    sys_var = ['AMSBIN', 'AMSHOME', 'AMSRESOURCES', 'SCMLICENSE']
    sys_var_exists = [item in os.environ for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            logger.error(f"The environment variable {sys_var[i]!r} has not been set")

    if not all(sys_var_exists):
        raise EnvironmentError('One or more AMS environment variables have not been set, '
                               'aborting AMS job')


def dict_concatenate(dict_list: Iterable[Mapping[KT, VT]]) -> Dict[KT, VT]:
    """Concatenates a list of dictionaries."""
    ret: Dict[KT, VT] = {}
    for item in dict_list:
        ret.update(item)
    return ret


def get_template(template_name: str, from_cat_data: bool = True) -> Settings:
    """Grab a yaml template and return it as Settings object."""
    if from_cat_data:
        path = join('data/templates', template_name)
        xs = pkg.resource_string('CAT', path)
        return Settings(yaml.load(xs.decode(), Loader=yaml.FullLoader))
    with open(template_name, 'r') as file:
        return Settings(yaml.load(file, Loader=yaml.FullLoader))


@overload
def validate_path(path: None) -> str:
    ...
@overload   # noqa: E302
def validate_path(path: P) -> P:
    ...
def validate_path(path):  # noqa: E302
    """Validate a provided directory path.

    Parameters
    ----------
    path : str
        Optional: A path to a directory.
        Will default to the current working directory if ``None``.

    Results
    -------
    |str|_
        Returns either **path** or the current working directory.

    Raises
    ------
    FileNotFoundError
        Raised if **path** cannot be found.

    NotADirectoryError
        Raised if **path** is not a directory.

    """
    if path in {None, '.', ''}:
        return os.getcwd()
    elif isdir(path):
        return path

    if not exists(path):
        raise FileNotFoundError(f"{path!r} not found")
    elif isfile(path):
        raise NotADirectoryError(f"{path!r} is not a directory")


@overload
def validate_core_atom(atom: int) -> int:
    ...
@overload  # noqa: E302
def validate_core_atom(atom: str) -> Union[Molecule, int]:
    ...
def validate_core_atom(atom):  # noqa: E302
    """Parse and validate the ``["optional"]["qd"]["dissociate"]["core_atom"]`` argument."""
    # Potential atomic number or symbol
    if isinstance(atom, int) or atom in PeriodicTable.symtonum:
        return to_atnum(atom)

    # Potential SMILES string
    try:
        mol = from_smiles(atom)
    except Exception as ex:
        raise ValueError(f'Failed to recognize {atom!r} as a valid atomic number, '
                         'atomic symbol or SMILES string') from ex

    # Double check the SMILES string:
    charge_dict = {}
    for at in mol:
        charge = at.properties.charge
        try:
            charge_dict[charge] += 1
        except KeyError:
            charge_dict[charge] = 1
    if 0 in charge_dict:
        del charge_dict[0]

    # Only a single charged atom is allowed
    if len(charge_dict) > 1:
        charge_count = sum([v for v in charge_dict.values()])
        raise MoleculeError(f'The SMILES string {repr(atom)} contains more than one charged atom: '
                            f'charged atom count: {charge_count}')
    return mol


def parallel_init(path: Union[str, 'os.PathLike[str]'],
                  folder: Union[str, 'os.PathLike[str]'],
                  hashing: str = 'input') -> None:
    """Construct a workdir with a thread-safe name.

    Serves as a wrapper around :func:`plams.init<scm.plams.core.functions.init>`.

    Paramaters
    ----------
    path : str
        The path to the PLAMS workdir.
    folder : str
        The name of the PLAMS workdir.
    hashing : str
        The type of hashing used by the PLAMS :class:`JobManager`.
        Only ``"input"`` is considered an acceptable value here.

    """
    if hashing != 'input':
        raise ValueError(f"Invalid value: {hashing!r}")

    folder_ = f'{folder}.{threading.get_ident()}'
    path_abs = join(path, folder_)
    if isdir(path_abs):
        rmtree(path_abs)
    init(path=path, folder=folder_)

    # Not quite thread-safe, but as long as all threads assign a single value
    # there should be no problem
    config.log.file = 3
    config.log.stdout = 0
    return None


def restart_init(path: Union[str, 'os.PathLike[str]'],
                 folder: Union[str, 'os.PathLike[str]'],
                 hashing: str = 'input') -> None:
    """Wrapper around the plams.init_ function; used for importing one or more previous jobs.

    All pickled .dill files in **path**/**folder**/ will be loaded into the
    :class:`GenJobManager` instance initiated by :func:`init()<scm.plams.core.functions.init>`.

    .. _plams.init: https://www.scm.com/doc/plams/components/functions.html#scm.plams.core.functions.init

    Note
    ----
    Previous jobs are stored in a more generator-esque manner in :attr:`GenJobManager.hashes` and
    :class:`Job` instances are thus created on demand rather than
    permanently storing them in memory.

    Paramaters
    ----------
    path : str
        The path to the PLAMS workdir.

    folder : str
        The name of the PLAMS workdir.

    hashing : str
        Optional: The type of hashing used by the PLAMS :class:`JobManager`.
        Accepted values are: ``"input"``, ``"runscript"``, ``"input+runscript"`` and ``None``.

    """  # noqa
    # Create a job manager
    settings = Settings({'counter_len': 3, 'hashing': hashing, 'remove_empty_directories': True})
    manager = GenJobManager(settings, path, folder, hashing)

    # Change the default job manager
    with open(os.devnull, 'w') as f_, redirect_stdout(f_):
        init()
    if os.path.isdir(config.default_jobmanager.workdir):
        rmtree(config.default_jobmanager.workdir)
    config.default_jobmanager = manager
    config.log.file = 3
    config.log.stdout = 0

    workdir = manager.workdir
    if not isdir(workdir):  # workdir does not exist, dont bother trying to load previous jobs
        os.mkdir(workdir)
        return None

    # Update the default job manager with previous Jobs
    for f in os.listdir(workdir):
        job_dir = join(workdir, f)
        if not isdir(job_dir):  # Not a directory; move along
            continue

        hash_file = join(job_dir, f + '.hash')
        if isfile(hash_file):  # Update JobManager.hashes
            manager.load_job(hash_file)

        # Grab the job name
        try:
            name, _num = f.rsplit('.', 1)
            num = int(_num)
        except ValueError:  # Jobname is not appended with a number
            name = f
            num = 1

        # Update JobManager.names
        try:
            manager.names[name] = max(manager.names[name], num)
        except KeyError:
            manager.names[name] = num
    return None


@overload
def cycle_accumulate(iterable: Iterable[T1]) -> Generator[T1, None, None]:
    ...
@overload  # noqa: E302
def cycle_accumulate(iterable: Iterable[T1], start: T1) -> Generator[T1, None, None]:
    ...
def cycle_accumulate(iterable, start=0):  # noqa: E302
    """Accumulate and return elements from **iterable** until it is exhausted.

    Then repeat (and keep accumulating) the sequence indefinitely.
    The elements of **iterable** must have access to the :func:`__iadd__` method,
    *e.g.* :class:`float`, :class:`int` or :class:`str`.

    """
    ret = start
    for i in cycle(iterable):
        ret += i
        yield ret


def iter_repeat(iterable: Iterable[T1], times: int) -> Iterator[T1]:
    """Iterate over an iterable and apply :func:`itertools.repeat` to each element.

    Examples
    --------
    .. code:: python

        >>> from CAT.utils import iter_repeat

        >>> iterable = range(3)
        >>> n = 2
        >>> iterator = iter_repeat(iterable, n)
        >>> for i in iterator:
        ...     print(i)
        0
        0
        1
        1
        2
        2

    Parameters
    ----------
    iterable : :class:`Iterable<collections.abc.Iterable>`
        An iterable.

    times : :class:`int`
        The number of times each element should be repeated.

    Returns
    -------
    :class:`Iterator<collections.abc.Iterator>`
        An iterator that yields each element from **iterable** multiple **times**.

    """
    return chain.from_iterable(repeat(i, times) for i in iterable)


def get_nearest_neighbors(center: Union[Molecule, np.ndarray],
                          neighbor: Union[Molecule, np.ndarray],
                          k: Union[int, Iterable[int]],
                          distance_upper_bound: float = np.inf,
                          return_dist: bool = False, **kwargs: Any) -> np.ndarray:
    r"""Return the :math:`k` nearest-neighbors (from **neighbor**) for all user-specified atoms in **center**.

    The Euclidean distance is herein used for identifying nearest-neighbors.

    Warning
    -------
    Missing neighbors are denoted with :code:`len(center)` if
    an insufficient amount of neighbors can be identified.

    Parameters
    ----------
    center : array-like [:class:`float`], shape :math:`(n, 3)`
        A 2D array-like object with the Cartesian coordinates of all central atoms.

    neighbor : array-like [:class:`float`], shape :math:`(m, 3)`
        A 2D array-like object with the Cartesian coordinates of all
        (potentially) neighboring atoms.

    k : :class:`int` or :data:`Iterable<collections.abc.Iterable>` [:class:`int`]
        The number of to-be returned nearest neighbors.
        If :math:`k` is an iterable than it will return all nearest-neighbors specified in there.
        For example, :code:`[1, 3, 5]` will returned the first, third and fifth nearest-neighbors.
        Note that counting starts from 1.

    distance_upper_bound : :class:`float`
        Limit the nearest-neighbor search to neighbors within a certain radius with respect to
        each atom in **center**.

    return_dist : :class:`bool`
        If ``True``, return both the indices of the :math:`k` nearest neighbors and the respective
        :math:`(n, k)` distance matrix.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for SciPy's :class:`cKDTree<scipy.spatial.cKDTree>` class.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`int`], shape :math:`(n, k)`
        A 2D array with indices of the  :math:`k` nearest neighbors.

    See Also
    --------
    :class:`cKDTree<scipy.spatial.cKDTree>`
        kd-tree for quick nearest-neighbor lookup.

    """  # noqa
    if center is neighbor:
        xyz1 = xyz2 = np.asarray(center)
    else:
        xyz1 = np.array(center, ndmin=2, copy=False)
        xyz2 = np.array(neighbor, ndmin=2, copy=False)

    if isinstance(k, abc.Iterable):
        k = as_nd_array(k, dtype=int)

    tree = cKDTree(xyz2, **kwargs)
    try:
        dist, idx = tree.query(xyz1, k=k, distance_upper_bound=distance_upper_bound)
    except ValueError as ex:
        _parse_value_error(ex, k)

    if idx.ndim == 1:  # Always return the indices as 2D array
        idx.shape += (1,)

    if return_dist:
        return dist, idx
    return idx


def _parse_value_error(ex: Exception, k: Any) -> NoReturn:
    """Post-process a :exc:`ValueError` raised by :func:`get_nearest_neighbors`."""
    if isinstance(k, abc.Iterable) and min(k) < 1:
        raise ValueError("All elements of 'k' must be larger than or equal to 1; "
                         f"observed minimum: {min(k)!r}") from ex
    elif hasattr(k, '__int__') and k < 1:
        raise ValueError("'k' must be larger than or equal to 1; "
                         f"observed value: {k!r}") from ex
    raise ex


class SetEnviron:
    """A reentrant, re-usable context manager for temporarily setting environment variables."""

    __slots__ = ('__weakref__', '_kwargs', '_kwargs_old')
    _kwargs: Mapping[str, str]
    _kwargs_old: Mapping[str, Optional[str]]

    def __init__(self, **kwargs: str) -> None:
        r"""Initialize the context manager.

        Parameters
        ----------
        \**kwargs : :class:`str`
            The to-be updated parameters.

        """
        self._kwargs = MappingProxyType(kwargs)
        self._kwargs_old = MappingProxyType({
            k: os.environ.get(k) for k in self._kwargs
        })

    def __enter__(self) -> None:
        """Enter the context manager."""
        os.environ.update(self._kwargs)

    def __exit__(
        self,
        __exc_type: Optional[Type[BaseException]],
        __exc_value: Optional[BaseException],
        __traceback: Optional[TracebackType],
    ) -> None:
        """Exit the context manager."""
        for k, v in self._kwargs_old.items():
            if v is None:
                # Use `pop` instead of `del` to ensure thread-safety
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def log_traceback_locals(logger: Logger, level: int = -1,
                         str_func: Callable[[object], str] = pformat) -> None:
    """Log all local variables at the specified traceback level.

    Parameters
    ----------
    logger : :class:`~logging.Logger`
        A logger for writing the local variables.
    level : :class:`int`
        The traceback level.
    str_func : :data:`Callable[[object], str]<typing.Callable>`
        The callable for creating the variables string representation.

    """
    try:
        local_dct = inspect.trace()[level].frame.f_locals
    except IndexError:
        i = operator.index(level)
        raise RuntimeError(f"No traceback was found at level {i}") from None

    for name, _value in local_dct.items():
        prefix = f"    {name}: {_value.__class__.__name__} = "
        n = len(prefix)
        value_str = textwrap.indent(str_func(_value), n * ' ')[n:].split('\n')
        value_str[0] = prefix + value_str[0]
        for v in value_str:
            logger.debug(v)


class KindEnum(enum.Enum):
    """An enum with different anchoring operation kinds (see :class:`AnchorTup`)."""

    FIRST = 0
    MEAN = 1
    MEAN_TRANSLATE = 2


class AllignmentEnum(enum.Enum):
    """An enum with different core vector orientations (see :class:`AllignmentTup`)."""

    SPHERE = 0
    SURFACE = 1


class MultiAnchorEnum(enum.Enum):
    """An enum with different actions for when ligands with multiple anchors are found."""

    ALL = 0
    FIRST = 1
    RAISE = 2


class AnchorTup(NamedTuple):
    """A named tuple with anchoring operation instructions."""

    mol: "None | Mol"
    group: "None | str" = None
    group_idx: Tuple[int, ...] = (0,)
    anchor_idx: Tuple[int, ...] = ()
    remove: "None | Tuple[int, ...]" = None
    kind: KindEnum = KindEnum.FIRST
    angle_offset: "None | float" = None
    dihedral: "None | float" = None
    group_format: FormatEnum = FormatEnum.SMILES
    multi_anchor_filter: MultiAnchorEnum = MultiAnchorEnum.ALL


class AllignmentTup(NamedTuple):
    """A named tuple with core vector orientation options."""

    kind: AllignmentEnum
    invert: bool


class MoleculeWarning(Warning, MoleculeError):
    """A :exc:`Warning` subclass for :class:`~scm.plams.mol.molecule.Molecule` related errors."""


def get_formula(mol: Molecule) -> str:
    """Backport of the PLAMS <= 1.5.1 ``Molecule.get_formula`` method.

    The resulting atoms are reported in alphabetical order,
    contrary to the Hill system (that prioritizes ``CH`` pairs) utilized after 1.5.1.
    """
    dct = Counter(at.symbol for at in mol)
    return "".join(f"{at}{i}" for at, i in sorted(dct.items()))
