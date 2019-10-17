from contextlib import AbstractContextManager
from typing import (
    Optional, Union, Dict, List, Hashable, MutableMapping, TypeVar, Iterable, Any, Collection
)

import numpy as np
import pandas as pd

from scm.plams import finish, Settings
from scm.plams.core.basejob import Job
from assertionlib.dataclass import AbstractDataClass

from CAT.logger import logger
from CAT.utils import restart_init
from CAT.settings_dataframe import SettingsDataFrame
from CAT.wip_workflows.workflow_dicts import ASA
from CAT.frozen_settings import FrozenSettings

T = TypeVar('T')

ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')


def concatenate_values(mapping: MutableMapping[Hashable, T], base_key: Any) -> List[T]:
    """

    Parameters
    ----------
    settings

    """
    i = 1
    ret = []
    while True:
        key = f'{base_key}{i}'
        try:
            value = mapping.pop(key)
        except KeyError:
            return ret
        else:
            ret.append(value)
            i += 1


class WorkFlow(AbstractDataClass):

    #: Map a name to a workflow template
    _TEMPLATE_MAPPING: FrozenSettings = FrozenSettings({
        'asa': ASA
    })

    #: Map a name to a workflow mol type
    _MOL_TYPE_MAPPING: FrozenSettings = FrozenSettings({
        'asa': 'qd'
    })

    #: Map a name to a workflow description
    _JOB_DESCRIPTION_MAPPING: FrozenSettings = FrozenSettings({
        'asa': 'ligand activation strain analyses'
    })

    #: Map a name to a new columns and their placeholder values
    _PULL_COLUMN_MAPPING: FrozenSettings = FrozenSettings({
        'asa': {ASA_INT: np.nan, ASA_STRAIN: np.nan, ASA_E: np.nan}
    })

    #: Map a name to a new columns and their placeholder values
    _PUSH_COLUMN_MAPPING: FrozenSettings = FrozenSettings({
        'asa': (ASA_INT, ASA_STRAIN, ASA_E)
    })

    # Get-only properties

    @property
    def template(self) -> Dict[str, List[str]]: return self._TEMPLATE_MAPPING[self.name]

    @property
    def mol_type(self) -> str: return self._MOL_TYPE_MAPPING[self.name]

    @property
    def job_description(self) -> str: return self._JOB_DESCRIPTION_MAPPING[self.name]

    @property
    def column_dict(self) -> str: return self._COLUMN_MAPPING[self.name]

    @property
    def columns(self) -> str: return self._PUSH_COLUMN_MAPPING[self.name]

    # Getter and setter properties

    @property
    def read(self) -> bool: return self._read

    @read.setter
    def read(self, value: Collection):
        self._read = self.db and self._MOL_TYPE_MAPPING[self.name] in value

    @property
    def write(self) -> bool: return self._write

    @write.setter
    def write(self, value: Collection):
        self._write = self.db and self._MOL_TYPE_MAPPING[self.name] in value

    @property
    def overwrite(self) -> bool: return self._overwrite

    @overwrite.setter
    def overwrite(self, value: Collection):
        self._overwrite = self.db and self._MOL_TYPE_MAPPING[self.name] in value

    # Methods and magic methods

    def __init__(self, name: str, db=None, read=None, write=None, overwrite=None, path=None,
                 folder=None, keep_files=None, jobs=None, settings=None, **kwargs) -> None:
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
        logger.info(f"Finishing {self.job_description}")
        self.to_db(df)

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
    def from_dict(cls, settings: Union[Settings, SettingsDataFrame], name: str) -> 'WorkFlow':
        if isinstance(settings, SettingsDataFrame):
            settings = settings.settings

        kwargs = {}

        # Raise a KeyError if a key cannot be found
        with Settings.supress_missing():
            # Extract the correct template
            try:
                template: Dict[str, List[str]] = cls._TEMPLATE_MAPPING[name]
            except KeyError as ex:
                err = (f"Invalid value for the 'name' parameter: {repr(name)}\n"
                       f"Allowed values: {', '.join(repr(k) for k in cls._MOL_TYPE_MAPPING)}")
                raise ValueError(err).with_traceback(ex.__traceback__)

            # Create a dictionary with keyword arguments
            for k, v in template.items():
                kwargs[k] = settings.get_nested(v)

        kwargs['jobs'] = concatenate_values(kwargs, 'job')
        kwargs['settings'] = concatenate_values(kwargs, 's')
        kwargs['name'] = name

        return super().from_dict(kwargs)


class PlamsInit(AbstractContextManager):
    def __init__(self, path: str, folder: str, hashing: Optional[str] = 'input'):
        self.path = path
        self.folder = folder
        self.hashing = hashing

    def __enter__(self) -> None:
        restart_init(self.path, self.folder, self.hashing)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        finish()
