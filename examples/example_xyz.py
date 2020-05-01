""" A module for grabbing the example multi-xyz file. """

__all__ = []

from os.path import (join, dirname)


def get_data_dir(path=None, name='data'):
    """ Return the path + name of the CAT/data directory. """
    path = path or join(dirname(dirname(__file__)), name)
    return join(path, name)
