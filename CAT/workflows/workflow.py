from shutil import rmtree
from pathlib import Path
from contextlib import AbstractContextManager
from typing import (
    Optional, Union, Dict, Hashable, MutableMapping, TypeVar, Iterable, Container, Tuple, Callable,
    Any
)

import pandas as pd

from scm.plams import finish, Settings
from scm.plams.core.basejob import Job
from assertionlib.dataclass import AbstractDataClass

from CAT.utils import restart_init
from CAT.logger import logger
from CAT.settings_dataframe import SettingsDataFrame
from CAT.workflows.workflow_dicts import finalize_templates as load_templates
from CAT.frozen_settings import FrozenSettings

T = TypeVar('T')

ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')
MOL = ('mol', '')


def _return_True(value: Any) -> bool: return True


def concatenate_values(mapping: MutableMapping[Hashable, T], base_key: Hashable,
                       filter_func: Callable[[Any], bool] = _return_True) -> Tuple[T, ...]:
    """Take a key and :meth:`pop<dict.pop>` all values from **mapping**.

    The popping will continue as long as :code:`base_key + str(i)` is available in the mapping,
    where ``i`` is an :class:`int` larger than 1.
    The value if ``i`` will start from 1 and increase by `+1` every iteration.

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
        A callable for truth-testing each extracted value.
        Values returning `False` will not be added to the to-be returned :class:`tuple`.

    Returns
    -------
    :class:`tuple`
        A tuple with all values popped from **mapping**.

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
    """A workflow manager."""

    #: Map a name to a workflow template.
    _WORKFLOW_TEMPLATES: FrozenSettings = load_templates()

    #: A context manager for supressing Pandas :exc:`SettingwithCopyWarning`.
    _SUPRESS_SETTINGWITHCOPYWARNING: AbstractContextManager = pd.option_context(
        'mode.chained_assignment', None
    )

    # Get-only properties

    @property
    def template(self) -> Dict[str, Tuple[str, ...]]:
        return self._WORKFLOW_TEMPLATES[self.name].template

    @property
    def mol_type(self) -> str:
        return self._WORKFLOW_TEMPLATES[self.name].mol_type

    @property
    def description(self) -> str:
        return self._WORKFLOW_TEMPLATES[self.name].description

    @property
    def import_columns(self) -> Dict[str, Tuple[str, str]]:
        return self._WORKFLOW_TEMPLATES[self.name].import_columns

    @property
    def export_columns(self) -> Tuple[Tuple[str, str], ...]:
        return self._WORKFLOW_TEMPLATES[self.name].export_columns

    # Getter and setter properties

    @property
    def read(self) -> bool: return self._read

    @read.setter
    def read(self, value: Container) -> None:
        self._read = self.db and self.mol_type in value

    @property
    def write(self) -> bool: return self._write

    @write.setter
    def write(self, value: Container) -> None:
        self._write = self.db and self.mol_type in value

    @property
    def overwrite(self) -> bool: return self._overwrite

    @overwrite.setter
    def overwrite(self, value: Container) -> None:
        self._overwrite = self.db and self.mol_type in value

    @property
    def jobs(self) -> Optional[Tuple[Job, ...]]: return self._jobs

    @jobs.setter
    def jobs(self, value: Optional[Tuple[Job, ...]]) -> None:
        self._jobs = None if not value else value

    @property
    def settings(self) -> Optional[Tuple[Settings, ...]]: return self._settings

    @settings.setter
    def settings(self, value: Optional[Tuple[Settings, ...]]) -> None:
        self._settings = None if not value else value

    # Methods and magic methods

    def __init__(self, name: str, db=None, read=False, write=False, overwrite=False, path='.',
                 keep_files=True, jobs=None, settings=None, **kwargs) -> None:
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
        self.jobs: Iterable[Job] = jobs
        self.settings: Iterable[Settings] = settings

        for k, v in kwargs.items():
            setattr(self, k, v)

    @AbstractDataClass.inherit_annotations()
    def _str_iterator(self):
        iterator = super()._str_iterator()
        return ((k.strip('_'), v) for k, v in iterator)

    def __call__(self, func: Callable, df: pd.DataFrame,
                 index: Optional[pd.Series] = slice(None), **kwargs) -> None:
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

        index : :class:`slice` or :class:`pandas.index` [:class:`bool`]
            An object for slicing the index of **df**.

        See Also
        --------
        :meth:`Workflow.from_db`:
            Returns a value for the **index** parameter.

        """
        slice1 = index, MOL
        slice2 = index, list(self.import_columns.keys())

        logger.info(f"Starting {self.description}")
        with PlamsInit(path=self.path, folder=self.name), self._SUPRESS_SETTINGWITHCOPYWARNING:
            value = func(df.loc[slice1], **vars(self), **kwargs)
            df.loc[slice2] = value
        logger.info(f"Finishing {self.description}\n")

    def from_db(self, df: pd.DataFrame) -> Union[slice, pd.Series]:
        """Ensure that all required keys are present in **df** and import results from the database.

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
        for key, value in self.import_columns.items():
            if key not in df:
                df[key] = value
        if not self.read:  # Nothing to see here, move along
            return slice(None)

        # Add all necasary keys to the DataFrame and import from the database
        self.db.from_csv(df, database=self.mol_type)

        # Return a new DataFrame slice based on previously calculated results
        if self.overwrite:
            return slice(None)
        else:
            ret = df[self.import_columns].isnull().any(axis=1)
            if ret.all():
                return slice(None)
            return ret

    def to_db(self, df: pd.DataFrame, job_recipe: Optional[dict] = None, **kwargs) -> None:
        """Export results to the database.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        job_recipe : :class:`dict`
            A (nested) dictionary with the used job settings.

        """
        # Write results to the database
        if self.write:
            self.db.update_csv(
                df,
                columns=list(self.export_columns),
                database=self.mol_type,
                overwrite=self.overwrite,
                job_recipe=job_recipe,
                **kwargs
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
        kwargs['jobs'] = concatenate_values(kwargs, 'job', filter_func=bool)
        kwargs['settings'] = concatenate_values(kwargs, 's', filter_func=bool)
        return cls.from_dict(kwargs)


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
