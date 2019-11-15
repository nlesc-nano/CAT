"""An example input file."""

from os.path import (dirname, join)

import yaml
from scm.plams import Settings

from CAT import base
from CAT.logger import logger


yaml_path = join(dirname(__file__), 'input_settings.yaml')
# yaml_path = '/Users/basvanbeek/Downloads/juliette.yaml'
with open(yaml_path, 'r') as file:
    arg = Settings(yaml.load(file, Loader=yaml.FullLoader))

try:
    qd_df, core_df, ligand_df = base.prep(arg)
except Exception as ex:
    logger.critical(f'{ex.__class__.__name__}: {ex}', exc_info=True)
    raise ex
