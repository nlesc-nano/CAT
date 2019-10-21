import operator
from shutil import rmtree
from pathlib import Path
from contextlib import AbstractContextManager
from typing import (
    Optional, Union, Dict, Hashable, MutableMapping, TypeVar, Iterable, Container, Tuple, Callable,
    Any, List, Type
)

import numpy as np
import pandas as pd

import rdkit
import qmflows
from rdkit.Chem.AllChem import UFFGetMoleculeForceField as UFF
from scm.plams import finish, Settings
from scm.plams.core.basejob import Job
from assertionlib.dataclass import AbstractDataClass

from ..utils import restart_init, _job_dict
from ..logger import logger
from ..settings_dataframe import SettingsDataFrame
from ..workflows.workflow_dicts import finalize_templates as load_templates
from ..frozen_settings import FrozenSettings

T = TypeVar('T')

ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')
MOL = ('mol', '')
OPT = ('opt', '')



def _return_True(value: Any) -> bool: return True
def _lt_0(value) -> int: return value < 0


def pop_and_concatenate(mapping: MutableMapping[Hashable, T], base_key: Hashable,
                        filter_func: Callable[[Any], bool] = _return_True) -> Tuple[T, ...]:
    """Take a key and :meth:`pop<dict.pop>` all values from **mapping**.

    The popping will continue as long as :code:`base_key + str(i)` is available in the mapping,
    where ``i`` is an :class:`int` larger than 1.
    The value if ``i`` will start from 1 and increase by `+ 1` every iteration.

    Examples
    --------
    .. code:: python

        >>> mapping: dict = {
        ...     'job1': 1,
        ...     'job2': 2,
        ...     'job3': 3,
        ...     'final_key': True
        ... }

        >>> base_key: str = 'job'
        >>> value_tuple: tuple = concatenate_values(mapping, base_key)
        >>> print(value_tuple)
        (1, 2, 3)

        >>> print(mapping)
        {'final_key': True}

    Parameters
    ----------
    mapping : :data:`MutableMapping`
        A dictionary or other mutable mapping.

    base_key : :data:`Hashable<typing.Hashable>`
        The base key which will be appended with successively increasing integers.

    filter_func : :data:`Callable<typing.Callable>`
        A callable for truth-testing each extracted **mapping** value.
        Values returning `False` will not be added to the to-be returned :class:`tuple`.

    Returns
    -------
    :class:`tuple`
        A tuple with values popped from **mapping**.

    """
    i = 1
    ret = []
    while True:
        key = f'{base_key}{i}'
        try:
            value = mapping.pop(key)
        except KeyError:
            return tuple(ret)
        else:
            if filter_func(value):
                ret.append(value)
            i += 1


class WorkFlow(AbstractDataClass):
    """A workflow manager.

    Examples
    --------
    Typical usage example:

    .. code:: python

        >>> import pandas as pd

        >>> # Prepare workflow parameters
        >>> df = pd.DataFrame(...)
        >>> settings = Settings(...)
        >>> def fancy_df_func(df, **kwargs):
        ...     pass

        >>> # Create the workflow
        >>> workflow = WorkFlow.from_template(settings, name='asa')
        >>> print(workflow)
        WorkFlow(
            name       = 'asa',
            db         = None,
            read       = False,
            write      = False,
            overwrite  = False,
            path       = '.',
            keep_files = True,
            jobs       = None,
            settings   = None
        )

        # Run the workflow
        >>> idx = workflow.from_db(df)
        >>> workflow(fancy_df_func, df, index=idx)
        >>> workflow.to_db(df)

    """

    #: Map a name to a workflow template.
    _WORKFLOW_TEMPLATES: FrozenSettings = load_templates()

    #: A context manager for supressing Pandas :exc:`SettingwithCopyWarning`.
    _SUPRESS_SETTINGWITHCOPYWARNING: AbstractContextManager = pd.option_context(
        'mode.chained_assignment', None
    )

    # Get-only properties

    @property
    def template(self) -> Dict[str, Tuple[str, ...]]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"template"``]."""
        return self._WORKFLOW_TEMPLATES[self.name].template

    @property
    def mol_type(self) -> str:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"mol_type"``]."""
        return self._WORKFLOW_TEMPLATES[self.name].mol_type

    @property
    def description(self) -> str:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"description"``]."""
        return self._WORKFLOW_TEMPLATES[self.name].description

    @property
    def import_columns(self) -> Dict[str, Tuple[str, str]]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"import_columns"``]."""  # noqa
        return self._WORKFLOW_TEMPLATES[self.name].import_columns

    @property
    def export_columns(self) -> Tuple[Tuple[str, str], ...]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"export_columns"``]."""  # noqa
        return self._WORKFLOW_TEMPLATES[self.name].export_columns

    # Getter and setter properties

    @property
    def read(self) -> bool: return self._read

    @read.setter
    def read(self, value: Union[bool, Container]) -> None:
        try:
            self._read = bool(self.db) and self.mol_type in value
        except TypeError:  # value is not a container
            self._read = bool(value)

    @property
    def write(self) -> bool: return self._write

    @write.setter
    def write(self, value: Union[bool, Container]) -> None:
        try:
            self._write = bool(self.db) and self.mol_type in value
        except TypeError:  # value is not a container
            self._write = bool(value)

    @property
    def overwrite(self) -> bool: return self._overwrite

    @overwrite.setter
    def overwrite(self, value: Union[bool, Container]) -> None:
        try:
            self._overwrite = bool(self.db) and self.mol_type in value
        except TypeError:  # value is not a container
            self._overwrite = bool(value)

    @property
    def jobs(self) -> Optional[Tuple[Job, ...]]: return self._jobs

    @jobs.setter
    def jobs(self, value: Optional[Iterable[Job]]) -> None:
        self._jobs = (None,) if value is None else tuple(value)

    @property
    def settings(self) -> Optional[Tuple[Settings, ...]]: return self._settings

    @settings.setter
    def settings(self, value: Optional[Iterable[Settings]]) -> None:
        self._settings = (None,) if value is None else tuple(value)

    # Methods and magic methods

    def __init__(self, name: str,
                 db: Optional['Database'] = None,
                 read: bool = False,
                 write: bool = False,
                 overwrite: bool = False,
                 path: str = '.',
                 keep_files: bool = True,
                 read_template: bool = True,
                 jobs: Optional[Iterable[Job]] = None,
                 settings: Optional[Iterable[Settings]] = None,
                 **kwargs: Any) -> None:
        if name not in self._WORKFLOW_TEMPLATES:
            err = (f"Invalid value for the 'name' parameter: {repr(name)}\n"
                   f"Allowed values: {', '.join(repr(k) for k in self._WORKFLOW_TEMPLATES)}")
            raise ValueError(err)

        self.name: str = name
        self.db = db

        self.read: bool = read
        self.write: bool = write
        self.overwrite: bool = overwrite

        self.path: str = path
        self.keep_files: bool = keep_files
        self.read_template: bool = read_template
        self.jobs: Iterable[Job] = jobs
        self.settings: Iterable[Settings] = settings

        for k, v in kwargs.items():
            setattr(self, k, v)

    @AbstractDataClass.inherit_annotations()
    def _str_iterator(self):
        iterator = super()._str_iterator()
        return ((k.strip('_'), v) for k, v in iterator)

    def __call__(self, func: Callable, df: pd.DataFrame,
                 idx_slice: Union[slice, pd.Series] = slice(None),
                 columns: Optional[List[Hashable]] = None, **kwargs) -> None:
        """Initialize the workflow.

        Parameters
        ----------
        func : :data:`Callable<typing.Callable>`
            A callable object which will recieve **df**, all :class:`WorkFlow` instance
            attributes and ***kwargs** as arguments.
            The callable is expected to conduct some operation and export the results **dfs'**
            :attr:`WorkFlow.import_columns` columns.

        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        idx_slice : :class:`slice` or :class:`pandas.Series` [:class:`bool`]
            An object for slicing the rows of **df** (*i.e.* :attr:`pandas.DataFrame.index`).

        See Also
        --------
        :meth:`Workflow.from_db`:
            Returns a value for the **index** parameter.

        """
        # Prepare slices
        slice1 = idx_slice, MOL
        slice2 = idx_slice, list(self.import_columns.keys()) if columns is None else columns

        # Run the workflow
        logger.info(f"Starting {self.description}")
        with PlamsInit(path=self.path, folder=self.name), self._SUPRESS_SETTINGWITHCOPYWARNING:
            self_vars = {k.strip('_'): v for k, v in vars(self).items()}
            value = func(df.loc[slice1], **self_vars, **kwargs)
            df.loc[slice2] = value
        logger.info(f"Finishing {self.description}\n")

    def from_db(self, df: pd.DataFrame) -> Union[slice, pd.Series]:
        """Ensure that all required keys are present in **df** and import from the database.

        Returns a :class:`pandas.index` with all to-be updated rows, as based on how many
        previous results were imported from :attr:`WorkFlow.db`.
        If no results were pulled from :attr:`WorkFlow.db` (or :attr:`WorkFlow.overwrite` is
        ``True``), then return :code:`slice(None)`.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        Returns
        -------
        :class:`pandas.Series` [:class:`bool`] or :class:`slice`
            A Series for slicing a part of **df** or  a :class:`slice` object for
            slicing the entirety of **df** (*i.e.* :code:`slice(0, None`).

        """
        # Add all necasary keys to the DataFrame
        for key, value in self.import_columns.items():
            if key not in df:
                df[key] = value

        if not self.read:  # Nothing to see here, move along
            return slice(None)

        # Import from the database
        with self._SUPRESS_SETTINGWITHCOPYWARNING:
            self.db.from_csv(df, database=self.mol_type)

        # Return a new DataFrame slice based on previously calculated results
        if self.overwrite:
            return slice(None)
        else:
            keys = list(self.import_columns.keys())
            return self._isnull(df, keys).any(axis=1)

    def to_db(self, df: pd.DataFrame, status: Optional[str] = None,
              job_recipe: Optional[dict] = None,
              idx_slice: Union[slice, pd.Series] = slice(None)) -> None:
        """Export results to the database.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        status : :class:`str`, optional
            Whether or not **df** contains structures resulting from a geometry optimization.

        idx_slice : :class:`slice` or :class:`pandas.Series` [:class:`bool`]
            An object for slicing the rows of **df** (*i.e.* :attr:`pandas.DataFrame.index`).

        job_recipe : :class:`dict`
            A (nested) dictionary with the used job settings.

        """
        # Dont export any settings columns if job_recipe is None
        # No job recipe == no settings to export anyway
        if job_recipe is None:
            export_columns = [i for i in self.export_columns if i[0] != 'settings']
        else:
            export_columns = list(self.export_columns)

        # Set the optimization status of the molecules to True
        if status == 'optimized':
            df.loc[idx_slice, OPT] = True

        # Write results to the database
        if self.write:
            with self._SUPRESS_SETTINGWITHCOPYWARNING:
                self.db.update_csv(
                    df.loc[idx_slice],
                    database=self.mol_type,
                    columns=export_columns,
                    overwrite=self.overwrite,
                    job_recipe=job_recipe,
                    status=status,
                )

        # Remove the PLAMS results directories
        if not self.keep_files:
            rmtree(Path(self.path) / self.name)

    @classmethod
    def from_template(cls, settings: Union[Settings, SettingsDataFrame], name: str) -> 'WorkFlow':
        """Construct a :class:`WorkFlow` instance from a |plams.Settings| object.

        Parameters
        ----------
        settings : |plams.Settings|
            A Settings instance with all CAT settings.
            Certain values are extracted from **settings** based on the supplied template
            (see **name**).

        name : :class:`str`
            The name of the settings template.

        See Also
        --------
        :attr:`WorkFlow._WORKFLOW_TEMPLATES`
            A dictionary with all available template names (*i.e.* its keys).

        """
        # Extract the settings object from the SettingsDataFrame
        if isinstance(settings, SettingsDataFrame):
            settings = settings.settings

        kwargs = {}
        kwargs['name'] = name

        # Raise a KeyError if a key cannot be found
        with Settings.supress_missing():
            try:  # Extract the correct template
                template: Dict[str, Tuple[str, ...]] = cls._WORKFLOW_TEMPLATES[name].template
            except KeyError as ex:
                err = (f"Invalid value for the 'name' parameter: {repr(name)}\n"
                       f"Allowed values: {', '.join(repr(k) for k in cls._WORKFLOW_TEMPLATES)}")
                raise ValueError(err).with_traceback(ex.__traceback__)

            # Create a dictionary with keyword arguments
            for k, v in template.items():
                kwargs[k] = settings.get_nested(v)

        # Post process all jobs and job settings
        kwargs['jobs'] = pop_and_concatenate(kwargs, 'job')
        kwargs['settings'] = pop_and_concatenate(kwargs, 's')
        return cls.from_dict(kwargs)

    def get_recipe(self) -> Settings:
        """Create a recipe for :meth:`WorkFlow.to_db`."""
        settings_names = [i[1:] for i in self.export_columns if i[0] == 'settings']
        uff_fallback = {
            'key': f'RDKit_{rdkit.__version__}', 'value': f'{UFF.__module__}.{UFF.__name__}'
        }

        ret = Settings()
        for name, job, settings in zip(settings_names, self.jobs, self.settings):
            # job is None, *i.e.* it's an RDKit UFF optimziation
            if job is None:
                ret[name].update(uff_fallback)
                continue

            settings = Settings(settings)
            if self.read_template:  # Update the settings using a QMFlows template
                template = qmflows.geometry['specific'][self.type_to_string(job)].copy()
                settings.soft_update(template)
            ret[name].key = job
            ret[name].value = settings
        return ret

    @staticmethod
    def _isnull(df: pd.DataFrame, columns: List[Hashable]) -> pd.DataFrame:
        """A more expansive version of the :func:`pandas.isnull` function.

        :class:`int` series now also return ``False`` if smaller than ``0`` and :class:`bool`
        series are simply inverted.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame.

        columns : :class:`list`
            A list of column keys from **df**.

        """
        dtype_dict = {
            np.dtype(bool): operator.invert,
            np.dtype(int): _lt_0,
            np.dtype(float): pd.isnull,
            np.dtype(object): pd.isnull
        }

        ret = pd.DataFrame(index=df.index)
        for key, series in df[columns].items():
            try:
                ret[key] = dtype_dict[series.dtype](series)
            except KeyError:  # Plan b
                ret[key] = series.isnull()
        return df

    @staticmethod
    def type_to_string(job: Union[Job, Type[Job]]) -> Optional[None]:
        """Turn a :class:`type` instance into a :class:`str`."""
        if not isinstance(job, type):
            job = type(job)
        try:
            return _job_dict[job]
        except KeyError:
            logger.error(f"No default settings available for type: '{job.__class__.__name__}'")
            return None


class PlamsInit(AbstractContextManager):
    """A context manager for calling :func:`.restart_init` and |plams.finish|."""

    def __init__(self, path: str, folder: str, hashing: str = 'input'):
        self.path = path
        self.folder = folder
        self.hashing = hashing

    def __enter__(self) -> None:
        """Enter the context manager; call :func:`.restart_init`."""
        restart_init(self.path, self.folder, self.hashing)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; call |plams.finish|."""
        finish()