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

def main(args=None):
    parser = argparse.ArgumentParser(
            prog='CAT',
            usage='init_cat my_settings_file.yaml',
            description='Description: This script initalizes \
            the Compound Attachment/Analysis Tool.'
    )

    parser.add_argument(
            'YAML',
            nargs='+',
            type=str,
            help='A .yaml file with the settings for CAT'
    )

    args = parser.parse_args(args)
    args = extract_args(args)

    CAT.base.prep(args)
