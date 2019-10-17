from contextlib import AbstractContextManager
from typing import Optional

from scm.plams import finish
from assertionlib.dataclass import AbstractDataClass

from CAT.utils import restart_init


class WorkFlow(AbstractDataClass):
    def __init__(self):
        self.db = None
        self.read: bool = False
        self.write: bool = False
        self.overwrite: bool = False

        self.path: str = ''
        self.folder: str = ''
        self.keep_files: bool = True
        self.jobs = ()
        self.settings = ()

    def __call__(self, df):
        self.from_db(df)
        with PlamsInit(self.path, self.folder):
            pass
        self.to_db(df)

    def from_db(self, df):
        pass

    def to_db(self, df):
        pass

    @classmethod
    def from_dict(cls, settings, template):
        dct = settings
        return super().from_dict(dct)


class PlamsInit(AbstractContextManager):
    def __init__(self, path: str, folder: str, hashing: Optional[str] = 'input'):
        self.path = path
        self.folder = folder
        self.hashing = hashing

    def __enter__(self) -> None:
        restart_init(self.path, self.folder, self.hashing)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        finish()
