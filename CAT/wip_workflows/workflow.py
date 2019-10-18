from shutil import rmtree
from pathlib import Path
from contextlib import AbstractContextManager
from typing import (
    Optional, Union, Dict, Hashable, MutableMapping, TypeVar, Iterable, Container, Tuple, Callable
)

import pandas as pd

from scm.plams import finish, Settings
from scm.plams.core.basejob import Job
from assertionlib.dataclass import AbstractDataClass

from CAT.utils import restart_init
from CAT.logger import logger
from CAT.settings_dataframe import SettingsDataFrame
from CAT.wip_workflows.workflow_dicts import finilize_templats as load_templates
from CAT.frozen_settings import FrozenSettings

T = TypeVar('T')

ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')


def concatenate_values(mapping: MutableMapping[Hashable, T], base_key: Hashable) -> Tuple[T, ...]:
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
            ret.append(value)
            i += 1


class WorkFlow(AbstractDataClass):

    #: Map a name to a workflow template
    _WORKFLOW_TEMPLATES: FrozenSettings = load_templates()

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
    def read(self, value: Container):
        self._read = self.db and self.mol_type in value

    @property
    def write(self) -> bool: return self._write

    @write.setter
    def write(self, value: Container):
        self._write = self.db and self.mol_type in value

    @property
    def overwrite(self) -> bool: return self._overwrite

    @overwrite.setter
    def overwrite(self, value: Container):
        self._overwrite = self.db and self.mol_type in value

    # Methods and magic methods

    def __init__(self, name: str, db=None, read=False, write=False, overwrite=False, path='.',
                 keep_files=True, jobs=None, settings=None, **kwargs) -> None:
        if name not in self._MOL_TYPE_MAPPING:
            err = (f"Invalid value for the 'name' parameter: {repr(name)}\n"
                   f"Allowed values: {', '.join(repr(k) for k in self._MOL_TYPE_MAPPING)}")
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

    def __call__(self, func: Callable, df: pd.DataFrame) -> None:
        """Initialize the workflow."""
        name = self.name

        self.from_db(df)
        logger.info(f"Starting {self.job_description}")

        with PlamsInit(path=self.path, folder=name):
            func(df, **vars(self))

        self.to_db(df)
        logger.info(f"Finishing {self.job_description}")

        if not self.keep_files:
            rmtree(Path(self.path) / name)

    def from_db(self, df: pd.DataFrame) -> None:
        """Import results from the database."""
        if not self.read:
            return

        for key, value in self.column_dict.items():
            if key not in df:
                df[key] = value
        self.db.from_csv(df, database=self.mol_type)

    def to_db(self, df: pd.DataFrame) -> None:
        """Export results to the database."""
        if not self.write:
            return

        # Update the database
        self.db.update_csv(
            df, columns=self.columns, database=self.mol_type, overwrite=self.overwrite
        )

    @classmethod
    def from_template(cls, settings: Union[Settings, SettingsDataFrame], name: str) -> 'WorkFlow':
        if isinstance(settings, SettingsDataFrame):
            settings = settings.settings

        kwargs = {}

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

        kwargs['jobs'] = concatenate_values(kwargs, 'job')
        kwargs['settings'] = concatenate_values(kwargs, 's')
        kwargs['name'] = name

        return cls.from_dict(kwargs)


class PlamsInit(AbstractContextManager):
    def __init__(self, path: str, folder: str, hashing: Optional[str] = 'input'):
        self.path = path
        self.folder = folder
        self.hashing = hashing

    def __enter__(self) -> None:
        restart_init(self.path, self.folder, self.hashing)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        finish()
