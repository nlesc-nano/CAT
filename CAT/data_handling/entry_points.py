#!/usr/bin/env python
"""
CAT.data_handling.input_parser
==============================

Entry point for CAT.

Index
-----
.. currentmodule:: CAT.data_handling.input_parser
.. autosummary::
    extract_args
    main

API
---
.. autofunction:: extract_args
.. autofunction:: main

"""

import argparse
from os import getcwd
from os.path import (join, exists)
from typing import (Optional, List)

import yaml

from scm.plams.core.settings import Settings
import CAT


def extract_args(args: Optional[List[str]] = None) -> Settings:
    """Extract and return all arguments."""
    input_file = args.YAML[0]
    if exists(input_file):
        pass
    elif exists(join(getcwd(), input_file)):
        input_file = join(getcwd(), input_file)
    else:
        input_file2 = join(getcwd(), input_file)
        err = f'No file found at {input_file} or {input_file2}'
        CAT.logger.logger('FileNotFoundError: ' + err)
        raise FileNotFoundError(err)

    with open(input_file, 'r') as file:
        return Settings(yaml.load(file, Loader=yaml.FullLoader))


def main(args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        prog='CAT',
        usage='init_cat my_settings_file.yaml',
        description=('Description: This script initalizes '
                     'the Compound Attachment Tool (CAT).')
    )

    parser.add_argument(
        'YAML', nargs=1, type=str, metavar='input.yaml',
        help='Required: A .yaml file with the settings for CAT'
    )

    args = parser.parse_args(args)
    CAT.base.prep(extract_args(args), return_mol=False)
