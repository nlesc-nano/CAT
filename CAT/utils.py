""" A module with miscellaneous functions. """

__all__ = ['check_sys_var', 'dict_concatenate', 'get_time', 'get_template']

import os
import time
import json
import yaml
import pkg_resources as pkg
from os.path import join

from scm.plams.core.settings import Settings


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


def get_template(template_name):
    """
    Grab a template and return it as Settings object.
    """
    path = join('data/templates', template_name)
    xs = pkg.resource_string('CAT', path)
    if 'json' in template_name:
        s = json.loads(xs.decode())
    elif 'yaml' in template_name or 'yml' in template_name:
        s = yaml.load(xs.decode())
    return Settings(s)
