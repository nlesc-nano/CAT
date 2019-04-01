""" An example input file. """

from os.path import (dirname, join)
import yaml
import CAT
from scm.plams.core.settings import Settings


yaml_path = join(dirname(__file__), 'input_settings.yaml')
with open(yaml_path, 'r') as file:
    arg = Settings(yaml.load(file, Loader=yaml.FullLoader))

qd_list, core_list, ligand_list = CAT.base.prep(arg)
