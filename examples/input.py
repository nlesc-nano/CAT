""" An example input file. """

from os.path import (dirname, join)
import yaml
import CAT
from scm.plams.core.settings import Settings


yaml_path = join(dirname(__file__), 'input_settings.yaml')
with open(yaml_path, 'r') as file:
    arg = Settings(yaml.load(file))

qd_df, core_df, ligand_df = CAT.base.prep(arg)
