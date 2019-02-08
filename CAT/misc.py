""" A module with miscellaneous functions. """

__all__ = ['check_sys_var', 'dict_concatenate', 'create_dir', 'get_time']

import os
import time


def get_time():
    """ Returns the current time as string. """
    return '[' + time.strftime('%H:%M:%S') + '] '


def check_sys_var():
    """
    Check if all ADF environment variables are set and if the 2018 version of ADF is installed.
    """
    sys_var = ['ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE']
    sys_var_exists = [item in os.environ for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            print(get_time() +
                  'WARNING: The environment variable ' + sys_var[i] + ' has not been set')
    if False in sys_var_exists:
        raise EnvironmentError(get_time() + 'One or more ADF environment variables have '
                               'not been set, aborting ADF job.')
    if '2018' not in os.environ['ADFHOME']:
        error = get_time() + 'No ADF version 2018 detected in ' + os.environ['ADFHOME']
        error += ', aborting ADF job.'
        raise ImportError(error)


def dict_concatenate(dic):
    """
    Concatenates a list of dictionaries.
    """
    concact_dic = {}
    for item in dic:
        concact_dic.update(item)
    return concact_dic


def create_dir(dir_name, path=os.getcwd()):
    """
    Creates a new directory if this directory does not yet exist.
    """
    if not os.path.exists(path):
        error = path + ' not found, aborting run'
        raise FileNotFoundError(error)
    dir_path = os.path.join(path, str(dir_name))
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path
