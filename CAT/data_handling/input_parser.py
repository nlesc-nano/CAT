#!/usr/bin/env python

from os import getcwd
from os.path import (join, exists)
import argparse
import yaml

from scm.plams.core.settings import Settings
import CAT


def extract_args(args):
    """ Extract and return all arguments. """
    input_file = args.settings[0]
    if exists(input_file):
        pass
    elif exists(join(getcwd(), input_file)):
        input_file = join(getcwd(), input_file)
    else:
        error = 'No file found at ' + input_file + ' or ' + join(getcwd(), input_file)
        raise FileNotFoundError(error)

    with open(input_file, 'r') as file:
        return Settings(yaml.load(file))


parser = argparse.ArgumentParser(
        prog='CAT',
        usage='%(prog)s settings',
        description='Start CAT: Compound Attachment/Analysis Tool.'
)

parser.add_argument(
        'settings',
        nargs='+',
        type=str,
        help='A YAML file with the settings for CAT'
)

# args = parser.parse_args()
# args = extract_args(args)

yml = '/Users/basvanbeek/Documents/GitHub/CAT/examples/input.yml'
with open(yml, 'r') as file:
    arg = Settings(yaml.load(file))

qd_list, core_list, ligand_list = CAT.base.prep(arg)
