from pathlib import Path

import yaml

from nanoutils import delete_finally
from scm.plams import Settings
from CAT.base import prep

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_thread_safe() -> None:
    """Tests for the ``optional.database.thread_safe`` option."""
    yaml_path = PATH / 'input4.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

    arg.path = PATH
    prep(arg)
