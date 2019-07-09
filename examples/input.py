"""An example input file."""

from os.path import (dirname, join)

import yaml
from scm.plams import Settings

from CAT import base


yaml_path = join(dirname(__file__), 'input_settings.yaml')
with open(yaml_path, 'r') as file:
    arg = Settings(yaml.load(file, Loader=yaml.FullLoader))

qd_df, core_df, ligand_df = base.prep(arg)
